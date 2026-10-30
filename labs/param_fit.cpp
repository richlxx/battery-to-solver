#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "dfn/fit/levenberg_marquardt.hpp"
#include "dfn/fit/spm_batch.hpp"
#include "dfn/linalg/dense.hpp"
#include "dfn/sim/config.hpp"

namespace {

using dfn::fit::BatchRow;
using dfn::fit::BatchSegment;

constexpr double kDefaultDt = 5.0;
constexpr int kDefaultNr = 20;
constexpr double kDefaultSigma = 1e-3;

struct Protocol {
  std::string name;
  std::vector<BatchSegment> segments;
};

Protocol make_protocol(std::string_view name) {
  if (name == "a") {
    return {"a", {{1200.0, 18.0}}};
  }
  if (name == "b") {
    return {"b", {{1200.0, 18.0}, {1800.0, 0.0}, {1200.0, -18.0}}};
  }
  if (name == "c") {
    return {"c",
            {{600.0, 0.0},
             {1200.0, 18.0},
             {1800.0, 0.0},
             {1200.0, -18.0},
             {1800.0, 0.0},
             {600.0, 36.0},
             {1800.0, 0.0},
             {300.0, 72.0},
             {1800.0, 0.0}}};
  }
  throw std::runtime_error("protocol must be a, b or c");
}

void write_rows(const std::vector<BatchRow>& rows, const std::string& path) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("cannot open " + path + " for writing");
  }
  file << "t_s,current_a_m2,voltage_v\n";
  file << std::setprecision(10);
  for (const BatchRow& row : rows) {
    file << row.t << ',' << row.current << ',' << row.voltage << '\n';
  }
}

std::vector<BatchRow> read_rows(const std::string& path) {
  std::ifstream file(path);
  if (!file) {
    throw std::runtime_error("cannot open " + path + " for reading");
  }
  std::string line;
  std::getline(file, line);
  std::vector<BatchRow> rows;
  while (std::getline(file, line)) {
    if (line.empty()) {
      continue;
    }
    std::stringstream stream(line);
    std::string field;
    BatchRow row;
    std::getline(stream, field, ',');
    row.t = std::stod(field);
    std::getline(stream, field, ',');
    row.current = std::stod(field);
    std::getline(stream, field, ',');
    row.voltage = std::stod(field);
    rows.push_back(row);
  }
  if (rows.size() < 3) {
    throw std::runtime_error("data file needs at least three samples");
  }
  return rows;
}

// theta = (ln kinetic scale, ln D_neg, ln D_pos). A single kinetic scale is
// used because the full-cell voltage cannot separate the two electrodes'
// rate constants; chapter 11 discusses this at length.
std::vector<double> log_parameters(const dfn::Params& params) {
  return {0.0, std::log(params.neg.d_s), std::log(params.pos.d_s)};
}

dfn::Params apply_parameters(dfn::Params params,
                             const std::vector<double>& theta) {
  const double kinetic_scale = std::exp(theta[0]);
  params.neg.k *= kinetic_scale;
  params.pos.k *= kinetic_scale;
  params.neg.d_s = std::exp(theta[1]);
  params.pos.d_s = std::exp(theta[2]);
  return params;
}

std::vector<std::vector<double>> sensitivity(const dfn::Params& params,
                                             const std::vector<BatchRow>& rows,
                                             int nr, double dt,
                                             double step = 1e-3) {
  const std::vector<double> base = log_parameters(params);
  const std::size_t m = rows.size();
  std::vector<std::vector<double>> columns(
      m, std::vector<double>(base.size(), 0.0));
  for (std::size_t j = 0; j < base.size(); ++j) {
    std::vector<double> up = base;
    std::vector<double> down = base;
    up[j] += step;
    down[j] -= step;
    const std::vector<double> v_up =
        dfn::fit::simulate_rows(apply_parameters(params, up), rows, nr, dt);
    const std::vector<double> v_down =
        dfn::fit::simulate_rows(apply_parameters(params, down), rows, nr, dt);
    for (std::size_t i = 0; i < m; ++i) {
      columns[i][j] = (v_up[i] - v_down[i]) / (2.0 * step);
    }
  }
  return columns;
}

double condition_number(const std::vector<std::vector<double>>& columns) {
  const std::size_t n = columns.front().size();
  std::vector<double> norm(n, 0.0);
  for (const auto& row : columns) {
    for (std::size_t j = 0; j < n; ++j) {
      norm[j] += row[j] * row[j];
    }
  }
  for (double& value : norm) {
    value = std::sqrt(value);
  }
  std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      if (norm[i] == 0.0 || norm[j] == 0.0) {
        return std::numeric_limits<double>::infinity();
      }
      double sum = 0.0;
      for (const auto& row : columns) {
        sum += (row[i] / norm[i]) * (row[j] / norm[j]);
      }
      matrix[i][j] = sum;
    }
  }

  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    std::size_t p = 0;
    std::size_t q = 1;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i + 1; j < n; ++j) {
        if (std::abs(matrix[i][j]) > off) {
          off = std::abs(matrix[i][j]);
          p = i;
          q = j;
        }
      }
    }
    if (off < 1e-12) {
      break;
    }
    const double angle = 0.5 * std::atan2(2.0 * matrix[p][q],
                                         matrix[q][q] - matrix[p][p]);
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    for (std::size_t k = 0; k < n; ++k) {
      const double mkp = matrix[k][p];
      const double mkq = matrix[k][q];
      matrix[k][p] = c * mkp - s * mkq;
      matrix[k][q] = s * mkp + c * mkq;
    }
    for (std::size_t k = 0; k < n; ++k) {
      const double mpk = matrix[p][k];
      const double mqk = matrix[q][k];
      matrix[p][k] = c * mpk - s * mqk;
      matrix[q][k] = s * mpk + c * mqk;
    }
  }
  double largest = 0.0;
  double smallest = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < n; ++i) {
    largest = std::max(largest, matrix[i][i]);
    smallest = std::min(smallest, std::abs(matrix[i][i]));
  }
  return smallest > 0.0 ? largest / smallest
                        : std::numeric_limits<double>::infinity();
}

void print_correlations(const std::vector<std::vector<double>>& columns) {
  const std::size_t n = columns.front().size();
  std::vector<double> norm(n, 0.0);
  for (const auto& row : columns) {
    for (std::size_t j = 0; j < n; ++j) {
      norm[j] += row[j] * row[j];
    }
  }
  for (double& value : norm) {
    value = std::sqrt(value);
  }
  std::cout << "C";
  for (std::size_t j = 0; j < n; ++j) {
    std::cout << ",p" << j;
  }
  std::cout << '\n';
  for (std::size_t i = 0; i < n; ++i) {
    std::cout << "p" << i;
    for (std::size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (const auto& row : columns) {
        sum += row[i] * row[j];
      }
      std::cout << ',' << (norm[i] * norm[j] > 0.0
                               ? sum / (norm[i] * norm[j])
                               : 0.0);
    }
    std::cout << '\n';
  }
}

std::vector<std::vector<double>> normalized_jacobian(
    const dfn::Params& params, const std::vector<BatchRow>& rows, int nr,
    double dt, double sigma) {
  const std::vector<std::vector<double>> columns =
      sensitivity(params, rows, nr, dt);
  const std::size_t n = columns.front().size();
  std::vector<std::vector<double>> matrix(rows.size(),
                                          std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < rows.size(); ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      matrix[i][j] = columns[i][j] / sigma;
    }
  }
  return matrix;
}

void print_confidence(const std::vector<std::vector<double>>& jacobian) {
  const std::size_t n = jacobian.front().size();
  std::vector<std::vector<double>> normal(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (const auto& row : jacobian) {
        normal[i][j] += row[i] * row[j];
      }
    }
  }
  std::vector<std::vector<double>> inverse(n, std::vector<double>(n, 0.0));
  try {
    for (std::size_t j = 0; j < n; ++j) {
      std::vector<double> unit(n, 0.0);
      unit[j] = 1.0;
      const std::vector<double> column = dfn::solve_dense(normal, unit);
      for (std::size_t i = 0; i < n; ++i) {
        inverse[i][j] = column[i];
      }
    }
  } catch (const std::exception&) {
    std::cout << "CI: unavailable (singular normal matrix)\n";
    return;
  }
  std::cout << "sigma_theta";
  for (std::size_t i = 0; i < n; ++i) {
    std::cout << ',' << std::sqrt(std::max(0.0, inverse[i][i]));
  }
  std::cout << '\n';
}

double gaussian(std::mt19937& rng, double sigma) {
  std::normal_distribution<double> distribution(0.0, sigma);
  return distribution(rng);
}

struct FitSummary {
  std::vector<double> theta;
  double cost = 0.0;
  int iterations = 0;
};

}  // namespace

int main(int argc, char** argv) {
  std::string config_path = "configs/fit_demo.yaml";
  std::string protocol_name = "c";
  std::string data_path;
  std::string out_path;
  std::string curve_path;
  std::string sensitivity_path;
  std::string theta_text;
  std::string action;
  int nr = kDefaultNr;
  double dt = kDefaultDt;
  double noise = kDefaultSigma;
  double sigma = kDefaultSigma;
  int starts = 8;
  unsigned seed = 7;

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
          << "param_fit: SPM parameter fitting lab\n\n"
          << "  --config PATH       truth/fixed parameters (default "
             "configs/fit_demo.yaml)\n"
          << "  --protocol a|b|c    experiment design (default c)\n"
          << "  --generate          write synthetic data (needs --out)\n"
          << "  --fit               fit the data (needs --data)\n"
          << "  --sensitivity       print S(t), correlations and condition\n"
          << "  --data PATH         input CSV\n"
          << "  --out PATH          output CSV for --generate\n"
          << "  --curve PATH        fitted voltage curve CSV\n"
          << "  --sensitivity-csv PATH  write S(t) columns to CSV\n"
          << "  --evaluate          score --theta on --data\n"
          << "  --theta S,DN,DP     kinetic scale and two D_s values\n"
          << "  --noise V           synthetic noise std (default 1e-3)\n"
          << "  --sigma V           fit noise std (default 1e-3)\n"
          << "  --seed N            random seed (default 7)\n"
          << "  --starts N          multi-start count (default 8)\n"
          << "  --dt S              sample interval (default 5 s)\n"
          << "  --nr N              radial cells (default 20)\n";
      return 0;
    }
    if (arg == "--config") {
      config_path = value("--config");
      continue;
    }
    if (arg == "--protocol") {
      protocol_name = value("--protocol");
      continue;
    }
    if (arg == "--generate") {
      action = "generate";
      continue;
    }
    if (arg == "--fit") {
      action = "fit";
      continue;
    }
    if (arg == "--sensitivity") {
      action = "sensitivity";
      continue;
    }
    if (arg == "--evaluate") {
      action = "evaluate";
      continue;
    }
    if (arg == "--theta") {
      theta_text = value("--theta");
      continue;
    }
    if (arg == "--data") {
      data_path = value("--data");
      continue;
    }
    if (arg == "--out") {
      out_path = value("--out");
      continue;
    }
    if (arg == "--curve") {
      curve_path = value("--curve");
      continue;
    }
    if (arg == "--sensitivity-csv") {
      sensitivity_path = value("--sensitivity-csv");
      continue;
    }
    if (arg == "--noise") {
      noise = std::stod(value("--noise"));
      continue;
    }
    if (arg == "--sigma") {
      sigma = std::stod(value("--sigma"));
      continue;
    }
    if (arg == "--seed") {
      seed = static_cast<unsigned>(std::stoul(value("--seed")));
      continue;
    }
    if (arg == "--starts") {
      starts = std::stoi(value("--starts"));
      continue;
    }
    if (arg == "--dt") {
      dt = std::stod(value("--dt"));
      continue;
    }
    if (arg == "--nr") {
      nr = std::stoi(value("--nr"));
      continue;
    }
    std::cerr << "error: unknown option " << arg << "\n";
    return 2;
  }

  try {
    const dfn::Params params = dfn::load_params_yaml(config_path);
    const Protocol protocol = make_protocol(protocol_name);

    if (action == "generate") {
      if (out_path.empty()) {
        throw std::runtime_error("--generate needs --out");
      }
      std::vector<BatchRow> rows =
          dfn::fit::rows_from_protocol(params, protocol.segments, nr, dt);
      std::mt19937 rng(seed);
      for (BatchRow& row : rows) {
        row.voltage += gaussian(rng, noise);
      }
      write_rows(rows, out_path);
      std::cout << "wrote " << rows.size() << " samples to " << out_path
                << " (protocol " << protocol.name << ", dt " << dt
                << " s, noise " << noise << " V)\n";
      return 0;
    }

    if (data_path.empty()) {
      throw std::runtime_error("--fit/--sensitivity need --data");
    }
    const std::vector<BatchRow> rows = read_rows(data_path);
    const double sample_dt = rows[1].t - rows[0].t;

    if (action == "sensitivity") {
      const auto columns = sensitivity(params, rows, nr, sample_dt);
      std::cout << std::setprecision(6);
      std::cout << "protocol " << protocol.name << " samples " << rows.size()
                << " kappa " << condition_number(columns) << '\n';
      print_correlations(columns);
      if (!sensitivity_path.empty()) {
        std::ofstream file(sensitivity_path);
        if (!file) {
          throw std::runtime_error("cannot open " + sensitivity_path +
                                   " for writing");
        }
        file << "t_s,s_k,s_dneg,s_dpos\n";
        file << std::setprecision(10);
        for (std::size_t i = 0; i < rows.size(); ++i) {
          file << rows[i].t << ',' << columns[i][0] << ',' << columns[i][1]
               << ',' << columns[i][2] << '\n';
        }
      }
      return 0;
    }

    if (action == "evaluate") {
      if (theta_text.empty()) {
        throw std::runtime_error("--evaluate needs --theta");
      }
      std::vector<double> theta;
      std::stringstream stream(theta_text);
      std::string field;
      while (std::getline(stream, field, ',')) {
        theta.push_back(std::log(std::stod(field)));
      }
      if (theta.size() != 3) {
        throw std::runtime_error("--theta needs scale,d_neg,d_pos");
      }
      const std::vector<double> predicted =
          dfn::fit::simulate_rows(apply_parameters(params, theta), rows, nr,
                                  sample_dt);
      double sum = 0.0;
      double worst = 0.0;
      for (std::size_t i = 0; i < rows.size(); ++i) {
        const double diff = predicted[i] - rows[i].voltage;
        sum += diff * diff;
        worst = std::max(worst, std::abs(diff));
      }
      const double rmse = std::sqrt(sum / static_cast<double>(rows.size()));
      std::cout << std::setprecision(8) << "evaluate protocol " << protocol.name
                << " samples " << rows.size() << " rmse_mV " << 1000.0 * rmse
                << " max_mV " << 1000.0 * worst << '\n';
      return 0;
    }

    if (action != "fit") {
      throw std::runtime_error("choose --generate, --fit, --sensitivity or --evaluate");
    }

    const std::vector<double> truth = log_parameters(params);
    std::mt19937 rng(seed + 1);
    std::uniform_real_distribution<double> spread(-1.2, 1.2);
    const auto simulate = [&](const std::vector<double>& theta) {
      return dfn::fit::simulate_rows(apply_parameters(params, theta), rows, nr,
                                     sample_dt);
    };
    const auto residual = [&](const std::vector<double>& theta,
                              std::vector<double>& values) {
      try {
        const std::vector<double> predicted = simulate(theta);
        values.resize(predicted.size());
        for (std::size_t i = 0; i < predicted.size(); ++i) {
          if (!std::isfinite(predicted[i])) {
            return false;
          }
          values[i] = (predicted[i] - rows[i].voltage) / sigma;
        }
        return true;
      } catch (const std::exception&) {
        return false;
      }
    };

    FitSummary best;
    best.cost = std::numeric_limits<double>::infinity();
    const auto clock_start = std::chrono::steady_clock::now();
    for (int start = 0; start < std::max(1, starts); ++start) {
      std::vector<double> theta = truth;
      if (start > 0) {
        for (double& value : theta) {
          value += spread(rng);
        }
      }
      dfn::LmOptions options;
      options.max_iterations = 60;
      const dfn::LmResult result =
          dfn::levenberg_marquardt(theta, residual, options);
      if (result.ok && result.cost < best.cost) {
        best.cost = result.cost;
        best.theta = result.theta;
        best.iterations = result.iterations;
      }
    }
    const double wall_ms = std::chrono::duration<double, std::milli>(
                               std::chrono::steady_clock::now() - clock_start)
                               .count();
    if (best.theta.empty()) {
      throw std::runtime_error("all LM starts failed");
    }

    const std::vector<double> predicted = simulate(best.theta);
    std::cout << std::setprecision(7);
    std::cout << "protocol=" << protocol.name << " samples=" << rows.size()
              << " starts=" << starts << " wall_ms=" << wall_ms << '\n';
    std::cout << "k_scale_fit=" << std::exp(best.theta[0])
              << " (k_neg=" << params.neg.k * std::exp(best.theta[0])
              << ", k_pos=" << params.pos.k * std::exp(best.theta[0]) << ")\n";
    std::cout << "d_neg_true=" << params.neg.d_s
              << " d_neg_fit=" << std::exp(best.theta[1])
              << " d_pos_true=" << params.pos.d_s
              << " d_pos_fit=" << std::exp(best.theta[2]) << '\n';
    std::cout << "relerr=";
    for (std::size_t j = 0; j < truth.size(); ++j) {
      const double reference = j == 0 ? 1.0 : std::exp(truth[j]);
      std::cout << (j ? "," : "")
                << std::abs(std::exp(best.theta[j]) / reference - 1.0);
    }
    std::cout << " chi2_red="
              << 2.0 * best.cost /
                     static_cast<double>(rows.size() - best.theta.size())
              << " iterations=" << best.iterations << '\n';
    const auto jacobian = normalized_jacobian(
        apply_parameters(params, best.theta), rows, nr, sample_dt, sigma);
    print_confidence(jacobian);
    std::cout << "kappa="
              << condition_number(sensitivity(
                     apply_parameters(params, best.theta), rows, nr, sample_dt))
              << '\n';
    if (!curve_path.empty()) {
      std::ofstream file(curve_path);
      if (!file) {
        throw std::runtime_error("cannot open " + curve_path + " for writing");
      }
      file << "t_s,current_a_m2,data_v,fit_v\n";
      file << std::setprecision(10);
      for (std::size_t i = 0; i < rows.size(); ++i) {
        file << rows[i].t << ',' << rows[i].current << ',' << rows[i].voltage
             << ',' << predicted[i] << '\n';
      }
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }
}
