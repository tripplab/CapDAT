#ifndef CAPDAT_FOLD_PATCH_METRICS_HPP
#define CAPDAT_FOLD_PATCH_METRICS_HPP

#include <vector>

#include "fold_patch_surface_fit.hpp"

struct FoldPatchAnalysisConfig;

struct PatchMetrics {
    double t_f = 0.0;
    double H_f = 0.0;
    double K_f = 0.0;
    std::vector<double> t_field;
    std::vector<double> H_field;
    std::vector<double> K_field;
};

[[nodiscard]] PatchMetrics computeMetrics(const PatchReconstructedSurfaces& surfaces,
                                          const FoldPatchAnalysisConfig& config);

#endif // CAPDAT_FOLD_PATCH_METRICS_HPP
