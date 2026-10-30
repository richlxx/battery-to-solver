#include <iostream>

#include "dfn/models/dfn.hpp"

int main() {
  const dfn::StateLayout layout = dfn::make_layout(8, 4, 8, 6, 6);
  dfn::DfnModel model(dfn::make_demo_params(), layout);
  const dfn::CscPattern pattern = model.jacobian_pattern();

  std::cout << "n," << pattern.n << "\n";
  std::cout << "nnz," << pattern.nnz() << "\n";

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
  std::cout << "block,cs_surf," << layout.size() << "," << model.size() << "\n";

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
