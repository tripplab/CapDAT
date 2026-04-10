#ifndef CAPDAT_FOLD_PATCH_FRAME_HPP
#define CAPDAT_FOLD_PATCH_FRAME_HPP

#include <string>

#if __has_include(<Eigen/Core>)
#include <Eigen/Core>
#else
namespace Eigen {
struct Vector3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    static Vector3d Zero() {
        return {};
    }
};

struct Matrix3d {
    double m[3][3] = {{1.0, 0.0, 0.0},
                      {0.0, 1.0, 0.0},
                      {0.0, 0.0, 1.0}};

    static Matrix3d Identity() {
        return {};
    }
};
} // namespace Eigen
#endif

#include "capsid.hpp"

struct FoldPatchAnalysisConfig;

struct LocalFrame {
    Eigen::Vector3d origin;
    Eigen::Matrix3d axes;
    std::string fold_name;
};

[[nodiscard]] LocalFrame buildFoldFrame(const Capsid& capsid,
                                        const FoldPatchAnalysisConfig& config);

#endif // CAPDAT_FOLD_PATCH_FRAME_HPP
