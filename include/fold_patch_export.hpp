#ifndef CAPDAT_FOLD_PATCH_EXPORT_HPP
#define CAPDAT_FOLD_PATCH_EXPORT_HPP

struct FoldPatchAnalysisResult;
struct FoldPatchAnalysisConfig;

void exportFoldPatch(const FoldPatchAnalysisResult& result,
                     const FoldPatchAnalysisConfig& config);

#endif // CAPDAT_FOLD_PATCH_EXPORT_HPP
