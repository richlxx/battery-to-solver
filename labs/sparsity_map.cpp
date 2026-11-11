#include <iostream>

#include "dfn/models/dfn.hpp"

// 实验（第 4 章 稀疏结构）：以 8/4/8 个 x 节点、正负极各 6 个径向节点构造 DFN
// 的 Jacobian 稀疏模式，打印各未知量块的偏移区间和全部非零元坐标 (row, col)。
// 输出可被脚本画成稀疏图：先给 n 与 nnz，再给 block 区间，最后是 row,col 列表。
int main() {
  // 布局：负极 8 个 x 节点、隔膜 4 个、正极 8 个，每颗粒径向 6 个节点。
  const dfn::StateLayout layout = dfn::make_layout(8, 4, 8, 6, 6);
  dfn::DfnModel model(dfn::make_demo_params(), layout);
  // 稀疏模式只登记非零元的位置（CSC 结构），不涉及具体数值。
  const dfn::CscPattern pattern = model.jacobian_pattern();

  // n 为总未知量数，nnz 为 Jacobian 中非零元个数。
  std::cout << "n," << pattern.n << "\n";
  std::cout << "nnz," << pattern.nnz() << "\n";

  // 每个 block 行给出该未知量块的全局下标区间 [起始, 结束)：
  // cs_neg / cs_pos 为固相浓度（按 x、r 展平），ce 为电解液浓度，
  // phis_neg / phis_pos / phie 为固相与液相电位。
  std::cout << "block,cs_neg," << layout.cs_neg_offset() << ","
            << layout.cs_pos_offset() << "\n";
  std::cout << "block,cs_pos," << layout.cs_pos_offset() << ","
            << layout.ce_offset() << "\n";
  std::cout << "block,ce," << layout.ce_offset() << ","
            << layout.phis_neg_offset() << "\n";
  std::cout << "block,phis_neg," << layout.phis_neg_offset() << ","
            << layout.phis_pos_offset() << "\n";
  std::cout << "block,phis_pos," << layout.phis_pos_offset() << ","
            << layout.phie_offset() << "\n";
  std::cout << "block,phie," << layout.phie_offset() << ","
            << layout.size() << "\n";
  // cs_surf：每个颗粒表面浓度的代数变量，位于 layout.size() 之后。
  std::cout << "block,cs_surf," << layout.size() << "," << model.size() << "\n";

  // 按列压缩（CSC）遍历：col_ptr 给出第 col 列在 row_index 中的区间，
  // 输出 (row, col) 对，行号来自 row_index。
  std::cout << "row,col\n";
  for (int col = 0; col < pattern.n; ++col) {
    for (int k = pattern.col_ptr[static_cast<std::size_t>(col)];
         k < pattern.col_ptr[static_cast<std::size_t>(col) + 1]; ++k) {
      std::cout << pattern.row_index[static_cast<std::size_t>(k)] << "," << col
                << "\n";
    }
  }
  return 0;
}
