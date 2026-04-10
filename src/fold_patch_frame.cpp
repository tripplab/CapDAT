#include "fold_patch_frame.hpp"

#include "fold_patch_analysis.hpp"

LocalFrame buildFoldFrame(const Capsid& capsid,
                          const FoldPatchAnalysisConfig& config) {
    (void)capsid;
    (void)config;

    LocalFrame frame;
    frame.origin = Eigen::Vector3d::Zero();
    frame.axes = Eigen::Matrix3d::Identity();
    return frame;
}
