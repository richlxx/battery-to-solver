#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cuda.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nvrtc.h>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "dfn/fit/spm_batch.hpp"
#include "dfn/sim/config.hpp"

namespace {

constexpr int kMaxCells = 32;
constexpr int kMaxSegments = 16;

struct DevElectrode {
  double thickness;
  double radius;
  double porosity;
  double c_max;
  double c0;
  double d_s;
  double k;
  double alpha_a;
  double alpha_c;
};

struct DevParams {
  DevElectrode neg;
  DevElectrode pos;
  double c_e;
  double temperature;
  double faraday;
  double gas;
};

struct DevSegment {
  double duration;
  double current;
};

const char* kKernelSource = R"NVRTC(
#define MAX_CELLS 32

struct DevElectrode {
  double thickness;
  double radius;
  double porosity;
  double c_max;
  double c0;
  double d_s;
  double k;
  double alpha_a;
  double alpha_c;
};

struct DevParams {
  DevElectrode neg;
  DevElectrode pos;
  double c_e;
  double temperature;
  double faraday;
  double gas;
};

struct DevSegment {
  double duration;
  double current;
};

__device__ double clamp_stoich(double cs, double cmax) {
  return fmin(fmax(cs, 1e-6 * cmax), cmax * (1.0 - 1e-6));
}

__device__ double neg_ocp(double x) {
  x = fmin(fmax(x, 0.0), 1.0);
  const double one_minus = 1.0 - x;
  return 0.05 + 0.8 * exp(-10.0 * x) + 0.1 * one_minus * one_minus;
}

__device__ double pos_ocp(double y) {
  y = fmin(fmax(y, 0.0), 1.0);
  return 3.9 + 0.2 * (1.0 - y);
}

__device__ double overpotential(const DevElectrode& e, double j, double cs,
                                double c_e, double T, double F, double R) {
  if (j == 0.0) {
    return 0.0;
  }
  cs = clamp_stoich(cs, e.c_max);
  const double i0 = F * e.k * pow(c_e, e.alpha_a) *
                    pow(e.c_max - cs, e.alpha_a) * pow(cs, e.alpha_c);
  if (i0 <= 0.0) {
    return 0.0;
  }
  const double a = e.alpha_a * F / (R * T);
  const double c = e.alpha_c * F / (R * T);
  double eta = 0.0;
  for (int iteration = 0; iteration < 60; ++iteration) {
    const double ea = exp(fmin(fmax(a * eta, -60.0), 60.0));
    const double ec = exp(fmin(fmax(-c * eta, -60.0), 60.0));
    const double residual = i0 / F * (ea - ec) - j;
    const double derivative = i0 / F * (a * ea + c * ec);
    const double step = residual / derivative;
    eta -= step;
    if (fabs(step) < 1e-12) {
      break;
    }
  }
  return eta;
}

__device__ void prepare(const DevElectrode& e, int nr, double dt,
                        double* lower, double* diag, double* upper,
                        double* alpha, double* surf_coeff, double* j_coeff,
                        double F) {
  const double dr = e.radius / static_cast<double>(nr);
  double area[MAX_CELLS + 1];
  double volume[MAX_CELLS];
  for (int i = 0; i <= nr; ++i) {
    const double r = i * dr;
    area[i] = 4.0 * 3.14159265358979323846 * r * r;
  }
  for (int i = 0; i < nr; ++i) {
    const double r0 = i * dr;
    const double r1 = (i + 1) * dr;
    volume[i] = 4.0 * 3.14159265358979323846 / 3.0 *
                (r1 * r1 * r1 - r0 * r0 * r0);
  }
  for (int i = 0; i < nr; ++i) {
    const double west = i > 0 ? e.d_s * area[i] / dr / volume[i] : 0.0;
    const double east =
        i + 1 < nr ? e.d_s * area[i + 1] / dr / volume[i] : 0.0;
    lower[i] = -dt * west;
    diag[i] = 1.0 - dt * (-(west + east));
    upper[i] = -dt * east;
  }
  *alpha = dr / (2.0 * e.d_s);
  *surf_coeff = -area[nr] / volume[nr - 1];
  const double a_s = 3.0 * (1.0 - e.porosity) / e.radius;
  *j_coeff = 1.0 / (a_s * F * e.thickness);
}

__device__ void advance(double* c, int nr, const double* lower,
                        const double* diag, const double* upper, double source) {
  double rhs[MAX_CELLS];
  double cp[MAX_CELLS];
  double dp[MAX_CELLS];
  for (int i = 0; i < nr; ++i) {
    rhs[i] = c[i];
  }
  rhs[nr - 1] += source;
  cp[0] = upper[0] / diag[0];
  dp[0] = rhs[0] / diag[0];
  for (int i = 1; i < nr; ++i) {
    const double denom = diag[i] - lower[i] * cp[i - 1];
    cp[i] = upper[i] / denom;
    dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denom;
  }
  c[nr - 1] = dp[nr - 1];
  for (int i = nr - 2; i >= 0; --i) {
    c[i] = dp[i] - cp[i] * c[i + 1];
  }
}

extern "C" __global__ void spm_batch_kernel(
    int n, const double* kinetic, const double* d_neg, const double* d_pos,
    DevParams base, int nr, double dt, const DevSegment* segments,
    int n_segments, double* probes) {
  const int sample = blockIdx.x * blockDim.x + threadIdx.x;
  if (sample >= n) {
    return;
  }
  DevParams params = base;
  const double scale = exp(kinetic[sample]);
  params.neg.k *= scale;
  params.pos.k *= scale;
  params.neg.d_s = d_neg[sample];
  params.pos.d_s = d_pos[sample];

  double lower_n[MAX_CELLS], diag_n[MAX_CELLS], upper_n[MAX_CELLS];
  double lower_p[MAX_CELLS], diag_p[MAX_CELLS], upper_p[MAX_CELLS];
  double alpha_n, alpha_p, surf_n, surf_p, jcoeff_n, jcoeff_p;
  prepare(params.neg, nr, dt, lower_n, diag_n, upper_n, &alpha_n, &surf_n,
          &jcoeff_n, params.faraday);
  prepare(params.pos, nr, dt, lower_p, diag_p, upper_p, &alpha_p, &surf_p,
          &jcoeff_p, params.faraday);

  double c_neg[MAX_CELLS];
  double c_pos[MAX_CELLS];
  for (int i = 0; i < nr; ++i) {
    c_neg[i] = params.neg.c0;
    c_pos[i] = params.pos.c0;
  }

  for (int s = 0; s < n_segments; ++s) {
    const double current = segments[s].current;
    const int steps =
        static_cast<int>(segments[s].duration / dt + 0.5);
    for (int k = 0; k < steps; ++k) {
      advance(c_neg, nr, lower_n, diag_n, upper_n,
              dt * surf_n * (current * jcoeff_n));
      advance(c_pos, nr, lower_p, diag_p, upper_p,
              dt * surf_p * (-current * jcoeff_p));
    }
    const double next =
        s + 1 < n_segments ? segments[s + 1].current : current;
    const double j_neg = next * jcoeff_n;
    const double j_pos = -next * jcoeff_p;
    const double cs_neg =
        clamp_stoich(c_neg[nr - 1] - alpha_n * j_neg, params.neg.c_max);
    const double cs_pos =
        clamp_stoich(c_pos[nr - 1] - alpha_p * j_pos, params.pos.c_max);
    const double eta_neg = overpotential(params.neg, j_neg, cs_neg, params.c_e,
                                         params.temperature, params.faraday,
                                         params.gas);
    const double eta_pos = overpotential(params.pos, j_pos, cs_pos, params.c_e,
                                         params.temperature, params.faraday,
                                         params.gas);
    const double voltage =
        pos_ocp(cs_pos / params.pos.c_max) -
        neg_ocp(cs_neg / params.neg.c_max) + eta_pos - eta_neg;
    probes[static_cast<long long>(s) * n + sample] = voltage;
  }
}
)NVRTC";

void check_cuda(CUresult status, const char* what) {
  if (status != CUDA_SUCCESS) {
    const char* message = nullptr;
    cuGetErrorString(status, &message);
    throw std::runtime_error(std::string(what) + ": " +
                             (message ? message : "CUDA error"));
  }
}

void check_nvrtc(nvrtcResult status, const char* what) {
  if (status != NVRTC_SUCCESS) {
    throw std::runtime_error(std::string(what) + ": " +
                             nvrtcGetErrorString(status));
  }
}

DevParams make_device(const dfn::Params& p) {
  DevParams d;
  d.neg = {p.neg.thickness, p.neg.radius, p.neg.porosity, p.neg.c_max,
           p.neg.c0,        p.neg.d_s,  p.neg.k,        p.neg.alpha_a,
           p.neg.alpha_c};
  d.pos = {p.pos.thickness, p.pos.radius, p.pos.porosity, p.pos.c_max,
           p.pos.c0,        p.pos.d_s,  p.pos.k,        p.pos.alpha_a,
           p.pos.alpha_c};
  d.c_e = p.elyte.c0;
  d.temperature = p.temperature;
  d.faraday = dfn::fit::kFaraday;
  d.gas = dfn::fit::kGasConstant;
  return d;
}

std::vector<dfn::fit::BatchSegment> protocol_c() {
  return {{600.0, 0.0},   {1200.0, 18.0}, {1800.0, 0.0},
          {1200.0, -18.0}, {1800.0, 0.0},  {600.0, 36.0},
          {1800.0, 0.0},  {300.0, 72.0},   {1800.0, 0.0}};
}

struct Samples {
  std::vector<double> kinetic;
  std::vector<double> d_neg;
  std::vector<double> d_pos;
};

Samples make_samples(int n, const dfn::Params& p, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  Samples samples;
  samples.kinetic.resize(static_cast<std::size_t>(n));
  samples.d_neg.resize(static_cast<std::size_t>(n));
  samples.d_pos.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    samples.kinetic[static_cast<std::size_t>(i)] = u(rng);
    samples.d_neg[static_cast<std::size_t>(i)] =
        p.neg.d_s * std::exp(u(rng));
    samples.d_pos[static_cast<std::size_t>(i)] =
        p.pos.d_s * std::exp(u(rng));
  }
  return samples;
}

std::vector<double> cpu_probes(const dfn::Params& base, const Samples& samples,
                               int n, int nr, double dt,
                               const std::vector<dfn::fit::BatchSegment>& segs,
                               int threads) {
  const std::vector<dfn::fit::BatchRow> rows =
      dfn::fit::rows_from_protocol(base, segs, nr, dt);
  std::vector<std::size_t> indices;
  std::size_t index = 0;
  for (const auto& segment : segs) {
    index += static_cast<std::size_t>(
        std::lround(segment.duration / dt));
    indices.push_back(index);
  }
  std::vector<double> out(static_cast<std::size_t>(n) * segs.size(), 0.0);
  auto run = [&](int i) {
    dfn::Params p = base;
    const double scale = std::exp(samples.kinetic[static_cast<std::size_t>(i)]);
    p.neg.k *= scale;
    p.pos.k *= scale;
    p.neg.d_s = samples.d_neg[static_cast<std::size_t>(i)];
    p.pos.d_s = samples.d_pos[static_cast<std::size_t>(i)];
    const std::vector<double> voltage =
        dfn::fit::simulate_rows(p, rows, nr, dt);
    for (std::size_t s = 0; s < indices.size(); ++s) {
      out[s * static_cast<std::size_t>(n) + static_cast<std::size_t>(i)] =
          voltage[indices[s]];
    }
  };
  if (threads <= 1) {
    for (int i = 0; i < n; ++i) {
      run(i);
    }
  } else {
    std::atomic<int> next{0};
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));
    for (int t = 0; t < threads; ++t) {
      workers.emplace_back([&]() {
        for (;;) {
          const int i = next.fetch_add(1);
          if (i >= n) {
            break;
          }
          run(i);
        }
      });
    }
    for (auto& worker : workers) {
      worker.join();
    }
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  std::string config_path = "configs/fit_demo.yaml";
  int n = 65536;
  int nr = 20;
  double dt = 5.0;
  int repeat = 3;
  int threads = 1;
  bool check = false;
  bool cpu_only = false;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    const auto value = [&](const char* name) -> const char* {
      if (i + 1 >= argc) {
        std::cerr << "error: " << name << " requires a value\n";
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "-h" || arg == "--help") {
      std::cout
          << "spm_batch: CUDA (NVRTC) batch SPM experiment\n\n"
          << "  --config PATH   parameters (default configs/fit_demo.yaml)\n"
          << "  --samples N     number of parameter samples (default 65536)\n"
          << "  --nr N          radial cells (default 20, max 32)\n"
          << "  --dt S          time step (default 5 s)\n"
          << "  --repeat N      timed repetitions (default 3)\n"
          << "  --threads N     CPU reference threads (default 1)\n"
          << "  --check         compare GPU probes with the CPU reference\n"
          << "  --cpu           run the CPU reference only\n";
      return 0;
    }
    if (arg == "--config") {
      config_path = value("--config");
      continue;
    }
    if (arg == "--samples") {
      n = std::stoi(value("--samples"));
      continue;
    }
    if (arg == "--nr") {
      nr = std::stoi(value("--nr"));
      continue;
    }
    if (arg == "--dt") {
      dt = std::stod(value("--dt"));
      continue;
    }
    if (arg == "--repeat") {
      repeat = std::stoi(value("--repeat"));
      continue;
    }
    if (arg == "--threads") {
      threads = std::max(1, std::stoi(value("--threads")));
      continue;
    }
    if (arg == "--check") {
      check = true;
      continue;
    }
    if (arg == "--cpu") {
      cpu_only = true;
      continue;
    }
    std::cerr << "error: unknown option " << arg << "\n";
    return 2;
  }

  try {
    if (nr < 1 || nr > kMaxCells) {
      throw std::runtime_error("nr must be in [1, 32]");
    }
    const dfn::Params params = dfn::load_params_yaml(config_path);
    const std::vector<dfn::fit::BatchSegment> segments = protocol_c();
    const Samples samples = make_samples(n, params, 11);

    if (cpu_only) {
      const auto start = std::chrono::steady_clock::now();
      const std::vector<double> reference =
          cpu_probes(params, samples, n, nr, dt, segments, threads);
      const double ms = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - start)
                            .count();
      std::cout << "cpu n=" << n << " probes=" << reference.size()
                << " threads=" << threads << " ms=" << ms
                << " per_sample_ms="
                << ms / static_cast<double>(n) << '\n';
      return 0;
    }

    check_cuda(cuInit(0), "cuInit");
    CUdevice device;
    check_cuda(cuDeviceGet(&device, 0), "cuDeviceGet");
    int major = 0;
    int minor = 0;
    check_cuda(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR,
                                    device),
               "compute capability");
    check_cuda(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR,
                                    device),
               "compute capability");
    CUcontext context;
    check_cuda(cuCtxCreate(&context, 0, device), "cuCtxCreate");

    nvrtcProgram program;
    check_nvrtc(nvrtcCreateProgram(&program, kKernelSource, "spm_batch.cu", 0,
                                   nullptr, nullptr),
                "nvrtcCreateProgram");
    const std::string arch_option = "--gpu-architecture=compute_" +
                                    std::to_string(major) +
                                    std::to_string(minor);
    const char* options[] = {"--std=c++17", arch_option.c_str()};
    const nvrtcResult compiled =
        nvrtcCompileProgram(program, 2, options);
    if (compiled != NVRTC_SUCCESS) {
      std::size_t log_size = 0;
      nvrtcGetProgramLogSize(program, &log_size);
      std::string log(log_size, '\0');
      nvrtcGetProgramLog(program, log.data());
      std::cerr << "NVRTC compile failed:\n" << log << "\n";
      throw std::runtime_error("NVRTC compile failed");
    }
    std::size_t ptx_size = 0;
    check_nvrtc(nvrtcGetPTXSize(program, &ptx_size), "nvrtcGetPTXSize");
    std::vector<char> ptx(ptx_size);
    check_nvrtc(nvrtcGetPTX(program, ptx.data()), "nvrtcGetPTX");
    nvrtcDestroyProgram(&program);

    CUmodule module;
    check_cuda(cuModuleLoadData(&module, ptx.data()), "cuModuleLoadData");
    CUfunction kernel;
    check_cuda(cuModuleGetFunction(&kernel, module, "spm_batch_kernel"),
               "cuModuleGetFunction");

    DevParams device_params = make_device(params);
    DevSegment host_segments[kMaxSegments];
    if (segments.size() > static_cast<std::size_t>(kMaxSegments)) {
      throw std::runtime_error("too many segments");
    }
    for (std::size_t s = 0; s < segments.size(); ++s) {
      host_segments[s] = {segments[s].duration, segments[s].current};
    }

    CUdeviceptr d_kinetic = 0;
    CUdeviceptr d_dneg = 0;
    CUdeviceptr d_dpos = 0;
    CUdeviceptr d_segments = 0;
    CUdeviceptr d_probes = 0;
    const std::size_t probe_count =
        static_cast<std::size_t>(n) * segments.size();
    check_cuda(cuMemAlloc(&d_kinetic, n * sizeof(double)), "cuMemAlloc");
    check_cuda(cuMemAlloc(&d_dneg, n * sizeof(double)), "cuMemAlloc");
    check_cuda(cuMemAlloc(&d_dpos, n * sizeof(double)), "cuMemAlloc");
    check_cuda(cuMemAlloc(&d_segments, sizeof(DevSegment) * segments.size()),
               "cuMemAlloc");
    check_cuda(cuMemAlloc(&d_probes, probe_count * sizeof(double)),
               "cuMemAlloc");
    check_cuda(cuMemcpyHtoD(d_kinetic, samples.kinetic.data(),
                            n * sizeof(double)),
               "cuMemcpyHtoD");
    check_cuda(cuMemcpyHtoD(d_dneg, samples.d_neg.data(), n * sizeof(double)),
               "cuMemcpyHtoD");
    check_cuda(cuMemcpyHtoD(d_dpos, samples.d_pos.data(), n * sizeof(double)),
               "cuMemcpyHtoD");
    check_cuda(cuMemcpyHtoD(d_segments, host_segments,
                            sizeof(DevSegment) * segments.size()),
               "cuMemcpyHtoD");

    const int block = 128;
    const int grid = (n + block - 1) / block;
    int n_segments = static_cast<int>(segments.size());
    void* args[] = {&n,        &d_kinetic, &d_dneg,      &d_dpos,
                    &device_params, &nr,  &dt,          &d_segments,
                    &n_segments, &d_probes};

    auto launch = [&]() {
      check_cuda(cuLaunchKernel(kernel, grid, 1, 1, block, 1, 1, 0, nullptr,
                                args, nullptr),
                 "cuLaunchKernel");
    };
    launch();
    check_cuda(cuCtxSynchronize(), "warm-up");

    CUevent start_event;
    CUevent stop_event;
    check_cuda(cuEventCreate(&start_event, CU_EVENT_DEFAULT), "cuEventCreate");
    check_cuda(cuEventCreate(&stop_event, CU_EVENT_DEFAULT), "cuEventCreate");
    check_cuda(cuEventRecord(start_event, nullptr), "cuEventRecord");
    for (int r = 0; r < repeat; ++r) {
      launch();
    }
    check_cuda(cuEventRecord(stop_event, nullptr), "cuEventRecord");
    check_cuda(cuEventSynchronize(stop_event), "cuEventSynchronize");
    float kernel_ms = 0.0f;
    check_cuda(cuEventElapsedTime(&kernel_ms, start_event, stop_event),
               "cuEventElapsedTime");
    kernel_ms /= static_cast<float>(repeat);

    std::vector<double> probes(probe_count);
    check_cuda(cuMemcpyDtoH(probes.data(), d_probes,
                            probe_count * sizeof(double)),
               "cuMemcpyDtoH");

    double max_diff = 0.0;
    if (check) {
      const int check_n = std::min(n, 16);
      const std::vector<double> reference =
          cpu_probes(params, samples, check_n, nr, dt, segments, 1);
      for (int i = 0; i < check_n; ++i) {
        for (std::size_t s = 0; s < segments.size(); ++s) {
          const double diff =
              std::abs(probes[s * static_cast<std::size_t>(n) +
                              static_cast<std::size_t>(i)] -
                       reference[s * static_cast<std::size_t>(check_n) +
                                 static_cast<std::size_t>(i)]);
          max_diff = std::max(max_diff, diff);
        }
      }
    }

    std::cout << std::setprecision(6);
    std::cout << "gpu n=" << n << " nr=" << nr << " dt=" << dt
              << " cc=" << major << minor
              << " kernel_ms=" << kernel_ms
              << " per_sample_us="
              << 1000.0 * kernel_ms / static_cast<double>(n);
    if (check) {
      std::cout << " max_diff=" << max_diff;
    }
    std::cout << '\n';

    cuMemFree(d_kinetic);
    cuMemFree(d_dneg);
    cuMemFree(d_dpos);
    cuMemFree(d_segments);
    cuMemFree(d_probes);
    cuEventDestroy(start_event);
    cuEventDestroy(stop_event);
    cuModuleUnload(module);
    cuCtxDestroy(context);
    return check && max_diff > 1e-9 ? 1 : 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }
}
