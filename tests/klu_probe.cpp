// 测试文件：KLU 稀疏线性求解器探针（对应 solvers/idaklu 的 CSC 装配路径）。
// 验证 SUNDIALS+KLU 可解小型稀疏系统，并能处理 ChainDae 的 2x2 Jacobian。

#include <cmath>
#include <iostream>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <algorithm>
#include <vector>

#include "dfn/models/chain_dae.hpp"

int main() {
  SUNContext context = nullptr;
  if (SUNContext_Create(SUN_COMM_NULL, &context) != 0) {
    std::cerr << "SUNContext_Create failed\n";
    return 1;
  }

  // 2x2 稀疏系统，A 按 CSC 存储，共 4 个非零元。
  const sunindextype n = 2;
  N_Vector x = N_VNew_Serial(n, context);
  N_Vector b = N_VNew_Serial(n, context);
  SUNMatrix A = SUNSparseMatrix(n, n, 4, CSC_MAT, context);

  // A = [[4, 1], [1, 3]]（CSC：按列依次存 row/val）。
  sunindextype* colptrs = SM_INDEXPTRS_S(A);
  sunindextype* rowvals = SM_INDEXVALS_S(A);
  sunrealtype* data = SM_DATA_S(A);
  colptrs[0] = 0;
  colptrs[1] = 2;
  colptrs[2] = 4;
  rowvals[0] = 0;
  data[0] = 4.0;
  rowvals[1] = 1;
  data[1] = 1.0;
  rowvals[2] = 0;
  data[2] = 1.0;
  rowvals[3] = 1;
  data[3] = 3.0;

  // 右端项 b = (7, 10)，对应解 x = (1, 3)。
  N_VConst(0.0, x);
  NV_Ith_S(b, 0) = 7.0;
  NV_Ith_S(b, 1) = 10.0;

  SUNLinearSolver solver = SUNLinSol_KLU(x, A, context);
  if (solver == nullptr) {
    std::cerr << "SUNLinSol_KLU failed\n";
    return 1;
  }

  if (SUNLinSolSetup(solver, A) != 0 ||
      SUNLinSolSolve(solver, A, x, b, 0.0) != 0) {
    std::cerr << "KLU solve failed\n";
    return 1;
  }

  // 直接法求解 2x2 系统，1e-12 的绝对容差只允许浮点舍入。
  const double x0 = NV_Ith_S(x, 0);
  const double x1 = NV_Ith_S(x, 1);
  const bool ok = std::abs(x0 - 1.0) < 1e-12 && std::abs(x1 - 3.0) < 1e-12;

  std::cout << "KLU solution: [" << x0 << ", " << x1 << "]\n";

  // 与 IDA Jacobian 回调相同的代码路径：用 dfn::CscPattern 构建 CSC 结构，
  // 其中一个微分行、一个代数行，验证 KLU 能分解该稀疏模式。
  {
    dfn::ChainDaeOptions options;
    options.cells = 1;
    options.lambda = 100.0;
    const dfn::ChainDae model(options);
    const dfn::CscPattern pattern = model.jacobian_pattern();
    SUNMatrix A2 = SUNSparseMatrix(2, 2, pattern.nnz(), CSC_MAT, context);
    // 把 dfn::CscPattern 的列指针与行下标原样复制到 SUNDIALS 稀疏矩阵。
    sunindextype* cp = SM_INDEXPTRS_S(A2);
    sunindextype* rv = SM_INDEXVALS_S(A2);
    for (int k = 0; k <= pattern.n; ++k) cp[k] = pattern.col_ptr[k];
    for (int k = 0; k < pattern.nnz(); ++k) rv[k] = pattern.row_index[k];

    // cj = 1e11 模拟极小时间步，检验主对角线量级相差很大时 KLU 仍可分解。
    dfn::CscBuilder builder(pattern);
    const std::vector<double> y = {0.0, 0.0};
    const std::vector<double> yp = {0.0, 1.0};
    model.jacobian(0.0, y.data(), yp.data(), 1.0e11, builder);
    std::copy(builder.values().begin(), builder.values().end(), SM_DATA_S(A2));

    N_Vector x2 = N_VNew_Serial(2, context);
    N_VConst(0.0, x2);
    SUNLinearSolver solver2 = SUNLinSol_KLU(x2, A2, context);
    const int setup_flag = SUNLinSolSetup(solver2, A2);
    std::cout << "assembled Jacobian setup flag: " << setup_flag << "\n";
    SUNLinSolFree(solver2);
    N_VDestroy(x2);
    SUNMatDestroy(A2);
    // setup_flag != 0 表示 KLU 符号/数值分解失败，例如结构奇异或内存不足。
    if (setup_flag != 0) {
      SUNLinSolFree(solver);
      SUNMatDestroy(A);
      N_VDestroy(x);
      N_VDestroy(b);
      SUNContext_Free(&context);
      return 1;
    }
  }

  SUNLinSolFree(solver);
  SUNMatDestroy(A);
  N_VDestroy(x);
  N_VDestroy(b);
  SUNContext_Free(&context);

  // 前面两块都成功后再确认数值解，避免提前销毁对象时遗漏错误。
  if (!ok) {
    std::cerr << "unexpected solution\n";
    return 1;
  }
  std::cout << "klu probe ok\n";
  return 0;
}
