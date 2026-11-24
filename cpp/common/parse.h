#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <map>
#include <limits>
#include <filesystem>
#include <charconv>
#include <optional>
#include <string_view>
#include <system_error>
#include <set>

#include "json.hpp"
#include "polyEval_class.h"
#include "error_bound.h"

using json = nlohmann::json;
namespace fs = std::filesystem;

// Build the coefficient filename for one OPS criterion.
std::string remez_coeff_filename(
    int p_num,
    int e_num,
    const std::string& criteria = "depth"
);

std::vector<std::vector<double>> parse_remez_coeff(
    int p_num,
    int e_num,
    bool print,
    std::string dir,
    std::string filename
);

std::vector<std::vector<double>> parse_remez_coeff(
    int p_num,
    int e_num,
    const std::string& criteria = "depth",
    bool print = false,
    const std::string& dir = "cache/ops"
);

int calcuate_remez_depth(std::vector<std::vector<double>> coeffs);
nlohmann::json load_decomp_cache(const std::string& filename);
std::string get_poly_type_key(const std::vector<double>& coeff);

// 현재 polyEval JSON 실행계획으로부터 Decomp 복원
std::shared_ptr<Decomp> reconstruct_decomp_from_cache(
    const json& plan,
    Poly poly,
    int depth=0
);

std::map<std::string, fs::path> findFilesByPrefix(
    const fs::path& directory,
    std::string_view prefix
);

std::optional<std::string> extractLastString(
    const fs::path& path,
    std::string_view prefix
);
