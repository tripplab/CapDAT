#include "fold_patch_extraction.hpp"

#include "fold_patch_analysis.hpp"

std::vector<PatchAtom> extractPatch(const Capsid& capsid,
                                    const LocalFrame& frame,
                                    const FoldPatchAnalysisConfig& config) {
    (void)capsid;
    (void)frame;
    (void)config;

    return {};
}

double vdwRadius(const std::string& element) {
    (void)element;
    return 1.70;
}
