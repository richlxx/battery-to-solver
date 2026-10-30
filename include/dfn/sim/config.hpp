#pragma once

#include <string>

#include "dfn/params.hpp"

namespace dfn {

// Load the flat YAML subset used by configs/demo.yaml: top-level sections
// (cell, negative, positive, separator, electrolyte) holding "key: number"
// lines. Unknown keys are rejected so that typos do not silently fall back to
// defaults.
Params load_params_yaml(const std::string& path);

}  // namespace dfn
