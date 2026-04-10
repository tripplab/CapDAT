#ifndef CAPDAT_FOLD_PATCH_ANALYSIS_HPP
#define CAPDAT_FOLD_PATCH_ANALYSIS_HPP

#include <string>
#include <vector>

#include "capsid.hpp"
#include "fold_patch_contact.hpp"
#include "fold_patch_export.hpp"
#include "fold_patch_extraction.hpp"
#include "fold_patch_frame.hpp"
#include "fold_patch_metrics.hpp"
#include "fold_patch_reporter.hpp"
#include "fold_patch_surface_fit.hpp"

struct FoldPatchAnalysisConfig {
    int fold_type = 2;
    int fold_index = 0;
    std::string fold_name = "2_0";
    double R = 12.0;
    double h = 2.0;
    double z_margin = 11.0;
    double soft = 0.35;
    double R_eval = 10.0;
    double min_valid_fraction = 0.70;
    double min_central_valid_fraction = 0.80;
    int min_atoms_in_patch = 200;
    double max_invalid_central_fraction = 0.20;
    bool export_patch_pdb = true;
    bool export_grid_pdb = true;
    bool export_surface_obj = true;
    bool export_surface_ply = true;
    bool export_surface_stl = true;
    bool export_fields_csv = true;
    bool export_summary_json = true;
    bool export_summary_csv = true;
    std::string out_prefix = "./anal";
};

struct FoldPatchAnalysisResult {
    std::string fold_name;
    std::string frame_description;
    std::string orientation_state;
    int atom_count = 0;
    double z_min = 0.0;
    double z_max = 0.0;
    int grid_nx = 0;
    int grid_ny = 0;
    int valid_nodes = 0;
    int invalid_nodes = 0;
    std::vector<PatchAtom> patch_atoms;
    PatchContactGrid raw_grid;
    PatchReconstructedSurfaces surfaces;
    PatchMetrics metrics;
    std::vector<std::string> exported_files;
    bool success = false;
    std::string error_message;
};

[[nodiscard]] FoldPatchAnalysisResult runFoldPatchAnalysis(const Capsid& capsid,
                                                           const FoldPatchAnalysisConfig& config);

#endif // CAPDAT_FOLD_PATCH_ANALYSIS_HPP
