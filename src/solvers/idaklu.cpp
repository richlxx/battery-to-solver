// SUNDIALS IDA 封装：把 DaeModel 的残差 / Jacobian / 事件回调接到 IDA，
// 提供 KLU 稀疏或稠密线性求解器、一致初值计算与输出采样。
// 所有 IDA 返回标志 < 0 都视为错误，并转换成 SolveReport 中的可读信息。

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

// 传给 IDA 的用户数据：模型、事件函数、Jacobian 模式、后端类型与错误信息。
// 回调是 C 函数，不能向外抛异常，因此错误先暂存在 error 里再在主流程检查。
struct CallbackData {
  const DaeModel* model = nullptr;
  const std::vector<EventFunction>* events = nullptr;
  CscPattern pattern;
  bool dense = false;
  std::string error;
};

// 把本项目的 CscPattern 写进 SUNDIALS 稀疏矩阵的列指针与行下标数组。
// 第 j 列的非零区间是 colptr[j] .. colptr[j+1]，行号按 CSC 约定升序排列。
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

// IDA 残差回调：把 N_Vector 的数据指针直接交给 model->residual 写入 rr。
// 回调内不能抛异常，捕获后记到 data->error 并返回 -1，让 IDA 报 IDA_RES_FAIL。
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

// IDA Jacobian 回调：IDA 期望得到 dF/dy + cj·dF/dy'，由 model->jacobian 装配。
// cj 由 IDA 按当前步长给出，与 y、y' 一起传入 model->jacobian 完成组装。
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
      // 稠密后端：SUNDIALS 稠密矩阵按列优先存放，与本项目稠密模式的值顺序相同，
      // 因此可以整块拷贝。
      std::copy(values.begin(), values.end(), SM_DATA_D(jac));
    } else {
      // IDA 每次求 Jacobian 前都会调用 SUNMatZero，对稀疏矩阵它还会清空下标数组，
      // 所以这里每次都重新写入结构（列指针 + 行下标）。
      fill_sparse_pattern(jac, data->pattern);
      std::copy(values.begin(), values.end(), SM_DATA_S(jac));
    }
  } catch (const std::exception& e) {
    data->error = e.what();
    return -1;
  }
  return 0;
}

// 事件（root）回调：依次计算每个事件函数 g_k(t, y, y') 写入 gout，
// IDA 检测到某个 g_k 变号时就把该时刻作为事件返回（见 ida_integrate 主流程）。
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

// 统一清理：按与创建相反的顺序释放 IDA 内存、线性求解器、矩阵、
// y/yp 向量与 SUNContext；传入空指针的项会被跳过，便于在任意失败点调用。
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

// 把 IDA 返回标志翻译成可读文本。这些字符串会写进报告与日志，故保持英文，
// 以便和 SUNDIALS 文档及原始输出逐条对照。
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

// 用 IDA 积分 F(t, y, y') = 0：从 (t0, y0, yp0) 出发积分到 t_end。
// 参数：model 为 DAE 模型；y0/yp0 为状态及状态导数初值（长度均为 n）；
//       t0/t_end 为时间区间；options 见 IdaOptions（容差、线性求解器、事件等）。
// 返回：SolveReport，含 ok、t_final、末态 y/yp、统计计数与错误信息。
// 向量约定：所有 N_Vector 都是串行稠密向量，下标与 model 的未知量顺序一致；
//           F 的第 i 个方程对应 y 的第 i 个分量，Jacobian 亦按此排列。
// 输入检查：模型非空、初值长度、时间区间；KLU 只能配合显式 Jacobian 使用；
// 而 atol 允许三种写法（空 / 长度 1 / 长度 n）。
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

  // 创建 SUNDIALS 运行上下文，后续所有对象（向量、矩阵、求解器）都挂在它上面。
  SUNContext ctx = nullptr;
  if (SUNContext_Create(SUN_COMM_NULL, &ctx) != 0) {
    report.message = "ida_integrate: SUNContext_Create failed";
    return report;
  }

  // 串行 N_Vector：y 为状态、yp 为状态导数，长度 n；把初值拷贝进去。
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

  // 选择线性代数后端：dense 时用全稠密模式，否则用模型自带的稀疏模式交给 KLU。
  CallbackData data;
  data.model = &model;
  data.dense = options.linear_solver == LinearSolverKind::dense;
  data.pattern = data.dense ? dense_pattern(n) : model.jacobian_pattern();

  // 初始化 IDA：注册残差回调，并给出 t0 处的 (y, y')。
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
  // 可选：注册事件（root）函数，g_k 变号时 IDA 停在事件时刻并返回 IDA_ROOT_RETURN。
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

  // 容差：单值 atol 用 IDASStolerances，逐分量 atol 用 IDASVtolerances。
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

  // 可选的最大步长限制：<= 0 表示完全交给 IDA 自适应控制。
  if (options.max_step > 0.0) {
    IDASetMaxStep(ida_mem, options.max_step);
  }

  SUNMatrix a = nullptr;
  SUNLinearSolver ls = nullptr;
  // 创建矩阵与线性求解器：稠密后端用 SUNDenseMatrix + SUNLinSol_Dense，
  // 稀疏后端用 SUNSparseMatrix（CSC，预置模式）+ SUNLinSol_KLU。
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

  // 把线性求解器与它使用的矩阵挂到 IDA 上；之后 IDA 负责调用分解与回代。
  flag = IDASetLinearSolver(ida_mem, ls, a);
  if (flag < 0) {
    destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
    report.error_code = flag;
    report.message = "IDASetLinearSolver: " + ida_error_message(flag);
    return report;
  }
  // 仅在用户提供显式 Jacobian 时注册回调；否则由 IDA 内部做差分近似（精度较低）。
  if (options.supply_jacobian) {
    flag = IDASetJacFn(ida_mem, jacobian_callback);
    if (flag < 0) {
      destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
      report.error_code = flag;
      report.message = "IDASetJacFn: " + ida_error_message(flag);
      return report;
    }
  }

  // 一致初值：让 F(t0, y0, y0') = 0 与代数约束同时成立。
  if (options.calculate_consistent_initial) {
    // IDACalcIC 使用 IDA_YA_YDP_INIT 时即使所有变量都是微分变量也必须给 id 向量，
    // 其中 1 表示微分变量、0 表示代数变量（由 model.variable_types 提供）。
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

  // 整理输出时刻：只保留 (t0, t_end] 内的采样点，排序去重，并保证包含 t_end。
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

  // 输出序列先记录初始时刻 t0 的 y（一致初值修正后的值）。
  sunrealtype t_final = t0;
  report.output_t.push_back(t0);
  report.output_y.emplace_back(N_VGetArrayPointer(y),
                               N_VGetArrayPointer(y) + n);

  // 逐个输出时刻调用 IDASolve；IDA 内部按自适应步长走到 tout，
  // 步长由误差检验自动放大/缩小，这里只控制“在哪些时刻取快照”。
  sunrealtype* y_data = N_VGetArrayPointer(y);
  for (double tout : times) {
    flag = IDASolve(ida_mem, tout, &t_final, y, yp, IDA_NORMAL);
    if (flag < 0) {
      break;
    }
    if (flag == IDA_ROOT_RETURN) {
      // 事件发生：读出哪个 root 函数变号，记录事件下标、时刻与状态后停止。
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
      // 已经到达 t_end（IDASetStopTime 设定）：记录最后一点后结束。
      break;
    }
  }

  // 取统计量：步数、残差/Jacobian 求值次数、线性分解次数与非线性迭代次数，
  // 用于判断求解难度（例如 Jacobian 求值过多往往说明步长被误差检验反复压缩）。
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
  // 非负标志（成功/停在 t_end/停在事件/仅警告）都算本次积分成功。
  report.ok = (flag == IDA_SUCCESS || flag == IDA_TSTOP_RETURN ||
               flag == IDA_ROOT_RETURN || flag == IDA_WARNING);
  report.message = ida_error_message(flag);
  // 回调里捕获到的 C++ 异常信息（若有）附加在消息末尾，便于定位模型错误。
  if (!data.error.empty()) {
    report.message += " (" + data.error + ")";
  }

  destroy(&y, &yp, &a, &ls, &ida_mem, &ctx);
  return report;
}

}  // namespace dfn
