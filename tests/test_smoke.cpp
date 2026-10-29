#include <iostream>
#include <string>

#include "dfn/version.hpp"

int main() {
  const std::string version(dfn::kVersion);
  if (version != "0.1.0") {
    std::cerr << "unexpected version: " << version << "\n";
    return 1;
  }
  std::cout << "dfn smoke ok (" << version << ")\n";
  return 0;
}
