#ifndef CAPDAT_FOLD_PATCH_SURFACE_FIT_HPP
#define CAPDAT_FOLD_PATCH_SURFACE_FIT_HPP

#include <vector>

#include "fold_patch_contact.hpp"

struct FoldPatchAnalysisConfig;

struct PatchReconstructedSurfaces {
    std::vector<double> z_upper;
    std::vector<double> z_lower;
    std::vector<double> z_mid;
    std::vector<double> dzdx_mid;
    std::vector<double> dzdy_mid;
    std::vector<double> d2zdx2_mid;
    std::vector<double> d2zdxdy_mid;
    std::vector<double> d2zdy2_mid;
    std::vector<bool> confidence_mask;
};

[[nodiscard]] PatchReconstructedSurfaces reconstructSurfaces(const PatchContactGrid& contact_grid,
                                                             const FoldPatchAnalysisConfig& config);

#endif // CAPDAT_FOLD_PATCH_SURFACE_FIT_HPP
