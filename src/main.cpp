#include <iostream>
#include <string_view>

#include "dfn/version.hpp"

namespace {

void print_help() {
  std::cout
      << "dfn-sim " << dfn::kVersion << "\n\n"
      << "Usage: dfn-sim [options]\n\n"
      << "Options:\n"
      << "  -h, --help     Show this message\n"
      << "      --version  Show the version\n"
      << "      --config   Load a YAML configuration (wired up in M2)\n";
}

}  // namespace

int main(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      print_help();
      return 0;
    }
    if (arg == "--version") {
      std::cout << dfn::kVersion << "\n";
      return 0;
    }
    if (arg == "--config") {
      if (i + 1 >= argc) {
        std::cerr << "error: --config requires a path\n";
        return 2;
      }
      std::cerr << "error: configuration loading is not implemented yet\n";
      return 2;
    }
    std::cerr << "error: unknown option " << arg << "\n";
    return 2;
  }

  print_help();
  return 0;
}
