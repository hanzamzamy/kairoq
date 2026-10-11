// Gerber RS-274X reader with X2 attributes (The Gerber Layer Format Specification,
// Ucamco, revision 2023.08).
#pragma once

#include "kairoq/feature.hpp"

#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace kairoq {

struct GerberImage {
    // File attributes (%TF): name without the leading dot -> fields.
    std::map<std::string, std::vector<std::string>> file_attributes;
    std::vector<Feature> features;
    std::vector<std::string> warnings;
    // True when at least one object carries a .N attribute.
    bool has_net_attributes = false;

    // %TF.FileFunction fields, or empty.
    std::vector<std::string> file_function() const;
    // %TF.FilePolarity: "Positive", "Negative" or empty.
    std::string file_polarity() const;
};

GerberImage parse_gerber(std::string_view text, const GeometryOptions& options = {});

// True if the text looks like a Gerber file (format or unit statement near the start).
bool looks_like_gerber(std::string_view text);

// Decode the \uXXXX escapes of an attribute field to UTF-8.
std::string decode_attribute_field(std::string_view field);

// Evaluate an aperture macro arithmetic expression ('x' is multiplication) with the
// given variables ($1 is vars[1]). Exposed for testing.
double evaluate_macro_expression(std::string_view expr, const std::vector<double>& vars);

}  // namespace kairoq
