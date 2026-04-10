#ifndef CAPDAT_FOLD_PATCH_EXTRACTION_HPP
#define CAPDAT_FOLD_PATCH_EXTRACTION_HPP

#include <string>
#include <vector>

#include "capsid.hpp"
#include "fold_patch_frame.hpp"

struct FoldPatchAnalysisConfig;

struct PatchAtom {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double vdw_radius = 0.0;
    std::string element;
    int serial = 0;
    std::string atom_name;
    std::string residue_id;
    std::string subunit_id;
};

[[nodiscard]] std::vector<PatchAtom> extractPatch(const Capsid& capsid,
                                                  const LocalFrame& frame,
                                                  const FoldPatchAnalysisConfig& config);

[[nodiscard]] double vdwRadius(const std::string& element);

#endif // CAPDAT_FOLD_PATCH_EXTRACTION_HPP
