#ifndef CAPDAT_FOLD_PATCH_REPORTER_HPP
#define CAPDAT_FOLD_PATCH_REPORTER_HPP

struct FoldPatchAnalysisResult;
struct FoldPatchAnalysisConfig;

void reportFoldPatchAnalysis(const FoldPatchAnalysisResult& result,
                             const FoldPatchAnalysisConfig& config,
                             bool verbose);

#endif // CAPDAT_FOLD_PATCH_REPORTER_HPP
