#include "dfn/sim/config.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace dfn {

namespace {

std::string trim(const std::string& text) {
  std::size_t begin = 0;
  while (begin < text.size() &&
         std::isspace(static_cast<unsigned char>(text[begin])) != 0) {
    ++begin;
  }
  std::size_t end = text.size();
  while (end > begin &&
         std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
    --end;
  }
  return text.substr(begin, end - begin);
}

double parse_number(const std::string& text, int line) {
  try {
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != text.size()) {
      throw std::invalid_argument("trailing characters");
    }
    return value;
  } catch (const std::exception&) {
    throw std::runtime_error("line " + std::to_string(line) +
                             ": expected a number, got '" + text + "'");
  }
}

void assign(Params& params, const std::string& section, const std::string& key,
            double value, int line) {
  ElectrodeParams* electrode = nullptr;
  if (section == "negative") {
    electrode = &params.neg;
  } else if (section == "positive") {
    electrode = &params.pos;
  }

  if (section == "cell") {
    if (key == "area") {
      params.area = value;
    } else if (key == "temperature") {
      params.temperature = value;
    } else {
      throw std::runtime_error("line " + std::to_string(line) +
                               ": unknown key 'cell." + key + "'");
    }
    return;
  }

  if (electrode != nullptr) {
    if (key == "thickness") {
      electrode->thickness = value;
    } else if (key == "radius") {
      electrode->radius = value;
    } else if (key == "porosity") {
      electrode->porosity = value;
    } else if (key == "bruggeman") {
      electrode->bruggeman = value;
    } else if (key == "c_max") {
      electrode->c_max = value;
    } else if (key == "c0") {
      electrode->c0 = value;
    } else if (key == "d_s") {
      electrode->d_s = value;
    } else if (key == "k") {
      electrode->k = value;
    } else if (key == "alpha_a") {
      electrode->alpha_a = value;
    } else if (key == "alpha_c") {
      electrode->alpha_c = value;
    } else if (key == "sigma") {
      electrode->sigma = value;
    } else {
      throw std::runtime_error("line " + std::to_string(line) +
                               ": unknown key '" + section + "." + key + "'");
    }
    return;
  }

  if (section == "separator") {
    if (key == "thickness") {
      params.sep.thickness = value;
    } else if (key == "porosity") {
      params.sep.porosity = value;
    } else if (key == "bruggeman") {
      params.sep.bruggeman = value;
    } else {
      throw std::runtime_error("line " + std::to_string(line) +
                               ": unknown key 'separator." + key + "'");
    }
    return;
  }

  if (section == "electrolyte") {
    if (key == "c0") {
      params.elyte.c0 = value;
    } else if (key == "d_e") {
      params.elyte.d_e = value;
    } else if (key == "t_plus") {
      params.elyte.t_plus = value;
    } else if (key == "kappa") {
      params.elyte.kappa = value;
    } else {
      throw std::runtime_error("line " + std::to_string(line) +
                               ": unknown key 'electrolyte." + key + "'");
    }
    return;
  }

  throw std::runtime_error("line " + std::to_string(line) + ": unknown section '" +
                           section + "'");
}

}  // namespace

Params load_params_yaml(const std::string& path) {
  std::ifstream file(path);
  if (!file) {
    throw std::runtime_error("cannot open parameter file " + path);
  }

  Params params = make_demo_params();
  std::string section;
  std::string line;
  int line_number = 0;

  while (std::getline(file, line)) {
    ++line_number;
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) {
      line.erase(comment);
    }
    line = trim(line);
    if (line.empty()) {
      continue;
    }

    const std::size_t colon = line.find(':');
    if (colon == std::string::npos) {
      throw std::runtime_error("line " + std::to_string(line_number) +
                               ": expected 'key: value' or 'section:'");
    }
    const std::string key = trim(line.substr(0, colon));
    const std::string rest = trim(line.substr(colon + 1));

    if (rest.empty()) {
      section = key;
      continue;
    }
    if (section.empty()) {
      throw std::runtime_error("line " + std::to_string(line_number) +
                               ": value outside of a section");
    }
    assign(params, section, key, parse_number(rest, line_number), line_number);
  }

  validate(params);
  return params;
}

}  // namespace dfn
