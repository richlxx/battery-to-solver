#include <string>

#include <pybind11/pybind11.h>

#include "dfn/version.hpp"

PYBIND11_MODULE(dfn_py, m) {
  m.attr("__version__") = std::string(dfn::kVersion);
}
