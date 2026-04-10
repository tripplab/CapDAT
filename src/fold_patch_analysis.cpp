#include "fold_patch_analysis.hpp"

#include <algorithm>

namespace {

std::string orientationStateToText(const Capsid::OrientationState& state) {
    if (state.in_original_parsed_frame) {
        return "original_parsed_frame";
    }
    if (state.reoriented_in_place) {
        return state.source_description.empty() ? "reoriented_in_place" : state.source_description;
    }
    return "unknown";
}

} // namespace

FoldPatchAnalysisResult runFoldPatchAnalysis(const Capsid& capsid,
                                             const FoldPatchAnalysisConfig& config) {
    FoldPatchAnalysisResult result;
    result.fold_name = config.fold_name;

    LocalFrame frame = buildFoldFrame(capsid, config);
    if (frame.fold_name.empty()) {
        result.success = false;
        result.error_message = "Could not resolve fold: " + config.fold_name;
        return result;
    }

    result.frame_description = frame.fold_name + " | origin=(" +
                               std::to_string(frame.origin.x) + "," +
                               std::to_string(frame.origin.y) + "," +
                               std::to_string(frame.origin.z) + ")";

    auto patch = extractPatch(capsid, frame, config);
    result.atom_count = static_cast<int>(patch.size());
    if (!patch.empty()) {
        result.z_min = std::min_element(patch.begin(), patch.end(),
                                        [](const PatchAtom& a, const PatchAtom& b) {
                                            return a.z < b.z;
                                        })
                           ->z;
        result.z_max = std::max_element(patch.begin(), patch.end(),
                                        [](const PatchAtom& a, const PatchAtom& b) {
                                            return a.z < b.z;
                                        })
                           ->z;
    }
    result.orientation_state = orientationStateToText(capsid.orientationState());

    if (result.atom_count < config.min_atoms_in_patch) {
        result.success = false;
        result.error_message = "Patch too small: " + std::to_string(result.atom_count) +
                               " atoms (min " + std::to_string(config.min_atoms_in_patch) + ")";
        return result;
    }

    result.patch_atoms = std::move(patch);

    exportFoldPatch(result, config);

    result.success = true;
    return result;
}
