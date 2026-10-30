#include "dfn/solvers/idaklu.hpp"

#include <ida/ida.h>
#include <ida/ida_ls.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <algorithm>
#include <exception>
#include <string>

namespace dfn {

namespace {

struct CallbackData {
  const DaeModel* model = nullptr;
  const std::vector<EventFunction>* events = nullptr;
  CscPattern pattern;
  bool dense = false;
  std::string error;
};

void fill_sparse_pattern(SUNMatrix a, const CscPattern& pattern) {
  sunindextype* col_ptr = SM_INDEXPTRS_S(a);
  sunindextype* rows = SM_INDEXVALS_S(a);
  for (int j = 0; j <= pattern.n; ++j) {
    col_ptr[j] = static_cast<sunindextype>(pattern.col_ptr[static_cast<std::size_t>(j)]);
  }
  for (int k = 0; k < pattern.nnz(); ++k) {
    rows[k] = static_cast<sunindextype>(pattern.row_index[static_cast<std::size_t>(k)]);
  }
}

int residual_callback(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rr,
                      void* user_data) {
  auto* data = static_cast<CallbackData*>(user_data);
  try {
    data->model->residual(t, N_VGetArrayPointer(y), N_VGetArrayPointer(yp),
                          N_VGetArrayPointer(rr));
  } catch (const std::exception& e) {
    data->error = e.what();
    return -1;
  }
  return 0;
}

int jacobian_callback(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp,
                      N_Vector /*r*/, SUNMatrix jac, void* user_data,
                      N_Vector /*tmp1*/, N_Vector /*tmp2*/,
                      N_Vector /*tmp3*/) {
  auto* data = static_cast<CallbackData*>(user_data);
  try {
    CscBuilder builder(data->pattern);
    data->model->jacobian(t, N_VGetArrayPointer(y), N_VGetArrayPointer(yp), cj,
                          builder);
    const std::vector<double>& values = builder.values();
    if (data->dense) {
      std::copy(values.begin(), values.end(), SM_DATA_D(jac));
    } else {
      // IDA calls SUNMatZero before every Jacobian, and for sparse matrices
      // that also clears the index arrays. Rewrite the structure each time.
      fill_sparse_pattern(jac, data->pattern);
      std::copy(values.begin(), values.end(), SM_DATA_S(jac));
    }
  } catch (const std::exception& e) {
    data->error = e.what();
    return -1;
  }
  return 0;
}

int root_callback(sunrealtype t, N_Vector y, N_Vector yp, sunrealtype* gout,
                  void* user_data) {
  auto* data = static_cast<CallbackData*>(user_data);
  try {
    const auto& events = *data->events;
    for (std::size_t k = 0; k < events.size(); ++k) {
      gout[k] = events[k](t, N_VGetArrayPointer(y), N_VGetArrayPointer(yp));
    }
  } catch (const std::exception& e) {
    data->error = e.what();
    return -1;
  }
  return 0;
}

void destroy(N_Vector* y, N_Vector* yp, SUNMatrix* a, SUNLinearSolver* ls,
             void** ida_mem, SUNContext* ctx) {
  if (ida_mem != nullptr && *ida_mem != nullptr) {
    IDAFree(ida_mem);
  }
  if (ls != nullptr && *ls != nullptr) {
    SUNLinSolFree(*ls);
  }
  if (a != nullptr && *a != nullptr) {
    SUNMatDestroy(*a);
  }
  if (y != nullptr && *y != nullptr) {
    N_VDestroy(*y);
  }
  if (yp != nullptr && *yp != nullptr) {
    N_VDestroy(*yp);
  }
  if (ctx != nullptr && *ctx != nullptr) {
    SUNContext_Free(ctx);
  }
}

}  // namespace

std::string ida_error_message(int flag) {
  switch (flag) {
    case IDA_SUCCESS:
      return "IDA_SUCCESS: integration completed";
    case IDA_TSTOP_RETURN:
      return "IDA_TSTOP_RETURN: stopped at the requested final time";
    case IDA_ROOT_RETURN:
      return "IDA_ROOT_RETURN: stopped at an event";
    case IDA_WARNING:
      return "IDA_WARNING: completed with a warning";
    case IDA_TOO_MUCH_WORK:
      return "IDA_TOO_MUCH_WORK: step count limit reached (stiff or wrong tolerances)";
    case IDA_TOO_MUCH_ACC:
      return "IDA_TOO_MUCH_ACC: accuracy could not be reached (tolerances too tight)";
    case IDA_ERR_FAIL:
      return "IDA_ERR_FAIL: error test failures (problem may be ill-posed)";
    case IDA_CONV_FAIL:
      return "IDA_CONV_FAIL: nonlinear solver failed to converge";
    case IDA_LINIT_FAIL:
      return "IDA_LINIT_FAIL: linear solver initialization failed";
    case IDA_LSETUP_FAIL:
      return "IDA_LSETUP_FAIL: linear setup failed, often a singular or incomplete Jacobian";
    case IDA_LSOLVE_FAIL:
      return "IDA_LSOLVE_FAIL: linear solve failed, often a singular or incomplete Jacobian";
    case IDA_RES_FAIL:
      return "IDA_RES_FAIL: the residual function returned an error";
    case IDA_REP_RES_ERR:
      return "IDA_REP_RES_ERR: repeated residual errors on the same step";
    case IDA_RTFUNC_FAIL:
      return "IDA_RTFUNC_FAIL: the root function failed";
    case IDA_CONSTR_FAIL:
      return "IDA_CONSTR_FAIL: solution violates the constraints";
    case IDA_FIRST_RES_FAIL:
      return "IDA_FIRST_RES_FAIL: the residual failed on the first step";
    case IDA_LINESEARCH_FAIL:
      return "IDA_LINESEARCH_FAIL: nonlinear linesearch failed (step too small)";
    case IDA_NO_RECOVERY:
      return "IDA_NO_RECOVERY: an error occurred and recovery is disabled";
    case IDA_NLS_INIT_FAIL:
      return "IDA_NLS_INIT_FAIL: nonlinear solver initialization failed";
    case IDA_NLS_SETUP_FAIL:
      return "IDA_NLS_SETUP_FAIL: nonlinear solver setup failed";
    case IDA_NLS_FAIL:
      return "IDA_NLS_FAIL: nonlinear solver failed to converge";
    case IDA_BAD_EWT:
      return "IDA_BAD_EWT: a weight in the error vector is non-positive";
    case IDA_TOO_CLOSE:
      return "IDA_TOO_CLOSE: t_end is too close to t0";
    case IDA_MEM_NULL:
      return "IDA_MEM_NULL: null solver memory";
    case IDA_MEM_FAIL:
      return "IDA_MEM_FAIL: memory allocation failed";
    case IDA_ILL_INPUT:
      return "IDA_ILL_INPUT: invalid input (check tolerances, matrices, callbacks)";
    case IDA_NO_MALLOC:
      return "IDA_NO_MALLOC: IDAInit was not called";
    default:
      return "unknown IDA flag " + std::to_string(flag);
  }
}

SolveReport ida_integrate(const DaeModel& model, const std::vector<double>& y0,
                          const std::vector<double>& yp0, double t0,
                          double t_end, const IdaOptions& options) {
  SolveReport report;
  const int n = model.size();
  if (n <= 0) {
    report.message = "ida_integrate: empty model";
    return report;
  }
  if (static_cast<int>(y0.size()) != n || static_cast<int>(yp0.size()) != n) {
    report.message = "ida_integrate: y0/yp0 size mismatch";
    return report;
  }
  if (t_end <= t0) {
    report.message = "ida_integrate: t_end must be greater than t0";
    return report;
  }
  if (!options.supply_jacobian && options.linear_solver == LinearSolverKind::klu) {
    report.message =
        "ida_integrate: KLU needs an explicit Jacobian; set supply_jacobian or "
        "use the dense linear solver";
    return report;
  }

  if (!options.atol.empty() && options.atol.size() != 1 &&
      static_cast<int>(options.atol.size()) != n) {
    report.message = "ida_integrate: atol must have size 1 or n";
    return report;
  }

  SUNContext ctx = nullptr;
  if (SUNContext_Create(SUN_COMM_NULL, &ctx) != 0) {
    report.message = "ida_integrate: SUNContext_Create failed";
    return report;
  }

  N_Vector y = N_VNew_Serial(n, ctx);
  N_Vector yp = N_VNew_Serial(n, ctx);
  if (y == nullptr || yp == nullptr) {
    destroy(&y, &yp, nullptr, nullptr, nullptr, &ctx);
    report.message = "ida_integrate: N_VNew_Serial failed";
    return report;
  }
  std::copy(y0.begin(), y0.end(), N_VGetArrayPointer(y));
  std::copy(yp0.begin(), yp0.end(), N_VGetArrayPointer(yp));

  void* ida_mem = IDACreate(ctx);
  if (ida_mem == nullptr) {
    destroy(&y, &yp, nullptr, nullptr, &ida_mem, &ctx);
    report.message = "ida_integrate: IDACreate failed";
    return report;
  }

  CallbackData data;
  data.model = &model;
  data.dense = options.linear_solver == LinearSolverKind::dense;
  data.pattern = data.dense ? dense_pattern(n) : model.jacobian_pattern();

  int flag = IDAInit(ida_mem, residual_callback, t0, y, yp);
  if (flag < 0) {
    destroy(&y, &yp, nullptr, nullptr, &ida_mem, &ctx);
    report.error_code = flag;
    report.message = "IDAInit: " + ida_error_message(flag);
    return report;
  }
  IDASetUserData(ida_mem, &data);
  IDASetStopTime(ida_mem, t_end);
  data.events = &options.events;
  if (!options.events.empty()) {
    flag = IDARootInit(ida_mem, static_cast<int>(options.events.size()),
                       root_callback);
    if (flag < 0) {
      destroy(&y, &yp, nullptr, nullptr, &ida_mem, &ctx);
      report.error_code = flag;
      report.message = "IDARootInit: " + ida_error_message(flag);
      return report;
    }
  }

  if (options.atol.size() <= 1) {
    const double atol = options.atol.empty() ? 1e-8 : options.atol[0];
    IDASStolerances(ida_mem, options.rtol, atol);
  } else {
    N_Vector atol = N_VNew_Serial(n, ctx);
    if (atol == nullptr) {
      destroy(&y, &yp, nullptr, nullptr, &ida_mem, &ctx);
      report.message = "ida_integrate: N_VNew_Serial failed for atol";
      return report;
    }
    std::copy(options.atol.begin(), options.atol.end(),
              N_VGetArrayPointer(atol));
    IDASVtolerances(ida_mem, options.rtol, atol);
    N_VDestroy(atol);
  }

  if (options.max_step > 0.0) {
    IDASetMaxStep(ida_mem, options.max_step);
  }

  SUNMatrix a = nullptr;
  SUNLinearSolver ls = nullptr;
  if (data.dense) {
    a = SUNDenseMatrix(n, n, ctx);
    if (a != nullptr) {
      ls = SUNLinSol_Dense(y, a, ctx);
    }
  } else {
    a = SUNSparseMatrix(n, n, data.pattern.nnz(), CSC_MAT, ctx);
    if (a != nullptr) {
      fill_sparse_pattern(a, data.pattern);
      ls = SUNLinSol_KLU(y, a, ctx);
    }
  }
  if (a == nullptr || ls == nullptr) {
    destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
    report.message = "ida_integrate: failed to create the matrix or linear solver";
    return report;
  }

  flag = IDASetLinearSolver(ida_mem, ls, a);
  if (flag < 0) {
    destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
    report.error_code = flag;
    report.message = "IDASetLinearSolver: " + ida_error_message(flag);
    return report;
  }
  if (options.supply_jacobian) {
    flag = IDASetJacFn(ida_mem, jacobian_callback);
    if (flag < 0) {
      destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
      report.error_code = flag;
      report.message = "IDASetJacFn: " + ida_error_message(flag);
      return report;
    }
  }

  if (options.calculate_consistent_initial) {
    // IDACalcIC with IDA_YA_YDP_INIT requires the id vector even when every
    // variable is differential.
    const std::vector<int> types = model.variable_types();
    N_Vector id = N_VNew_Serial(n, ctx);
    if (id == nullptr) {
      destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
      report.message = "ida_integrate: N_VNew_Serial failed for id";
      return report;
    }
    for (int i = 0; i < n; ++i) {
      NV_Ith_S(id, i) = static_cast<sunrealtype>(types[static_cast<std::size_t>(i)]);
    }
    IDASetId(ida_mem, id);
    N_VDestroy(id);

    flag = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, t_end);
    if (flag < 0) {
      const std::string detail = data.error.empty() ? "" : " (" + data.error + ")";
      destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
      report.error_code = flag;
      report.message = "IDACalcIC: " + ida_error_message(flag) + detail;
      return report;
    }
  }

  std::vector<double> times;
  for (double value : options.output_times) {
    if (value > t0 && value <= t_end) {
      times.push_back(value);
    }
  }
  std::sort(times.begin(), times.end());
  times.erase(std::unique(times.begin(), times.end()), times.end());
  if (times.empty() || times.back() < t_end) {
    times.push_back(t_end);
  }

  sunrealtype t_final = t0;
  report.output_t.push_back(t0);
  report.output_y.emplace_back(N_VGetArrayPointer(y),
                               N_VGetArrayPointer(y) + n);

  sunrealtype* y_data = N_VGetArrayPointer(y);
  for (double tout : times) {
    flag = IDASolve(ida_mem, tout, &t_final, y, yp, IDA_NORMAL);
    if (flag < 0) {
      break;
    }
    if (flag == IDA_ROOT_RETURN) {
      std::vector<int> roots(options.events.size(), 0);
      if (IDAGetRootInfo(ida_mem, roots.data()) == IDA_SUCCESS) {
        for (std::size_t k = 0; k < roots.size(); ++k) {
          if (roots[k] != 0) {
            report.event_index = static_cast<int>(k);
            break;
          }
        }
      }
      report.t_event = t_final;
      report.output_t.push_back(t_final);
      report.output_y.emplace_back(y_data, y_data + n);
      break;
    }
    report.output_t.push_back(t_final);
    report.output_y.emplace_back(y_data, y_data + n);
    if (flag == IDA_TSTOP_RETURN) {
      break;
    }
  }

  long int steps = 0;
  long int residual_evals = 0;
  long int jacobian_evals = 0;
  long int linear_setups = 0;
  long int nonlinear_iterations = 0;
  IDAGetNumSteps(ida_mem, &steps);
  IDAGetNumResEvals(ida_mem, &residual_evals);
  IDAGetNumJacEvals(ida_mem, &jacobian_evals);
  IDAGetNumLinSolvSetups(ida_mem, &linear_setups);
  IDAGetNumNonlinSolvIters(ida_mem, &nonlinear_iterations);

  report.error_code = flag;
  report.steps = static_cast<int>(steps);
  report.residual_evals = static_cast<int>(residual_evals);
  report.jacobian_evals = static_cast<int>(jacobian_evals);
  report.linear_setups = static_cast<int>(linear_setups);
  report.nonlinear_iterations = static_cast<int>(nonlinear_iterations);
  report.t_final = t_final;
  report.y.assign(y_data, y_data + n);
  report.yp.assign(N_VGetArrayPointer(yp), N_VGetArrayPointer(yp) + n);
  report.ok = (flag == IDA_SUCCESS || flag == IDA_TSTOP_RETURN ||
               flag == IDA_ROOT_RETURN || flag == IDA_WARNING);
  report.message = ida_error_message(flag);
  if (!data.error.empty()) {
    report.message += " (" + data.error + ")";
  }

  destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
  return report;
}

}  // namespace dfn
