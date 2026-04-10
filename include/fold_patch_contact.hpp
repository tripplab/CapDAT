#ifndef CAPDAT_FOLD_PATCH_CONTACT_HPP
#define CAPDAT_FOLD_PATCH_CONTACT_HPP

#include <vector>

#include "fold_patch_extraction.hpp"

struct FoldPatchAnalysisConfig;

struct PatchContactGrid {
    std::vector<double> node_x;
    std::vector<double> node_y;
    std::vector<double> raw_upper;
    std::vector<double> raw_lower;
    std::vector<bool> valid_mask;
    std::vector<int> diagnostic_flags;
};

[[nodiscard]] PatchContactGrid buildContactGrid(const std::vector<PatchAtom>& atoms,
                                                const FoldPatchAnalysisConfig& config);

#endif // CAPDAT_FOLD_PATCH_CONTACT_HPP
