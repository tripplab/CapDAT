#include "fold_patch_frame.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "fold_patch_analysis.hpp"
#include "geometry_symmetry.hpp"

namespace {

double dot(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Eigen::Vector3d cross(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

double norm(const Eigen::Vector3d& v) {
    return std::sqrt(dot(v, v));
}

Eigen::Vector3d normalize(const Eigen::Vector3d& v) {
    const double n = norm(v);
    if (n <= 1e-12) {
        return {0.0, 0.0, 0.0};
    }
    return {v.x / n, v.y / n, v.z / n};
}

Eigen::Vector3d subtract(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Eigen::Vector3d scale(const Eigen::Vector3d& v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}

Eigen::Vector3d add(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Eigen::Vector3d toEigen(const geometry_symmetry::Vector3& v) {
    return {v.x, v.y, v.z};
}

bool parseFoldName(const std::string& fold_name, int& fold_type, int& fold_index) {
    const std::size_t underscore = fold_name.find('_');
    if (underscore == std::string::npos || underscore == 0 || underscore == fold_name.size() - 1) {
        return false;
    }

    try {
        fold_type = std::stoi(fold_name.substr(0, underscore));
        fold_index = std::stoi(fold_name.substr(underscore + 1));
    } catch (const std::exception&) {
        return false;
    }
    return true;
}

double orthonormalError(const Eigen::Matrix3d& axes) {
    double err = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double dot_col = 0.0;
            for (int k = 0; k < 3; ++k) {
                dot_col += axes.m[k][i] * axes.m[k][j];
            }
            const double target = (i == j) ? 1.0 : 0.0;
            const double diff = dot_col - target;
            err += diff * diff;
        }
    }
    return std::sqrt(err);
}

Eigen::Vector3d centroidOfAtomsClosestToAxis(const Capsid& capsid,
                                             const Eigen::Vector3d& axis_unit) {
    std::vector<std::pair<double, Eigen::Vector3d>> candidates;
    candidates.reserve(capsid.atomCount());

    for (const Chain& chain : capsid.chains()) {
        for (const Residue& residue : chain.residues()) {
            for (const Atom& atom : residue.atoms()) {
                if (atom.isHetatm()) {
                    continue;
                }

                const Eigen::Vector3d p{atom.x(), atom.y(), atom.z()};
                const double axial = dot(p, axis_unit);
                if (axial <= 0.0) {
                    continue;
                }

                const Eigen::Vector3d perp = subtract(p, scale(axis_unit, axial));
                candidates.emplace_back(norm(perp), p);
            }
        }
    }

    if (candidates.empty()) {
        return Eigen::Vector3d::Zero();
    }

    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    const std::size_t near_axis_count = std::min<std::size_t>(400, candidates.size());
    std::vector<Eigen::Vector3d> near_axis_points;
    near_axis_points.reserve(near_axis_count);
    for (std::size_t i = 0; i < near_axis_count; ++i) {
        near_axis_points.push_back(candidates[i].second);
    }

    std::sort(near_axis_points.begin(), near_axis_points.end(), [&](const Eigen::Vector3d& a,
                                                                    const Eigen::Vector3d& b) {
        return dot(a, axis_unit) > dot(b, axis_unit);
    });

    const std::size_t n = std::min<std::size_t>(120, near_axis_points.size());
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    for (std::size_t i = 0; i < n; ++i) {
        sum = add(sum, near_axis_points[i]);
    }
    return scale(sum, 1.0 / static_cast<double>(n));
}

} // namespace

LocalFrame buildFoldFrame(const Capsid& capsid,
                          const FoldPatchAnalysisConfig& config) {
    LocalFrame frame;
    frame.origin = Eigen::Vector3d::Zero();
    frame.axes = Eigen::Matrix3d::Identity();
    frame.fold_name.clear();

    int fold_type = 0;
    int fold_index = 0;
    std::string fold_name = config.fold_name;
    if (!parseFoldName(config.fold_name, fold_type, fold_index)) {
        std::cerr << "[fold_patch] WARNING: unexpected fold_name format '" << config.fold_name
                  << "'; falling back to fold_type/fold_index config values\n";
        fold_type = config.fold_type;
        fold_index = config.fold_index;
        fold_name = std::to_string(fold_type) + "_" + std::to_string(fold_index);
    }

    geometry_symmetry::Vector3 fold_axis_raw{};
    try {
        const geometry_symmetry::FoldDefinition& fold = geometry_symmetry::foldByName(fold_name);
        fold_axis_raw = geometry_symmetry::foldUnitVector(fold.name);
    } catch (const std::exception&) {
        return frame;
    }

    const Eigen::Vector3d z_axis = normalize(toEigen(fold_axis_raw));

    Eigen::Vector3d reference = {1.0, 0.0, 0.0};
    if (std::fabs(dot(z_axis, reference)) > 0.99) {
        reference = {0.0, 1.0, 0.0};
    }

    Eigen::Vector3d x_axis = subtract(reference, scale(z_axis, dot(reference, z_axis)));
    x_axis = normalize(x_axis);
    const Eigen::Vector3d y_axis = normalize(cross(z_axis, x_axis));

    // NOTE: geometry_symmetry does not expose a capsid-fitted fold vertex point,
    // so we use the centroid of protein atoms closest to the positive fold axis.
    frame.origin = centroidOfAtomsClosestToAxis(capsid, z_axis);
    frame.axes.m[0][0] = x_axis.x;
    frame.axes.m[1][0] = x_axis.y;
    frame.axes.m[2][0] = x_axis.z;
    frame.axes.m[0][1] = y_axis.x;
    frame.axes.m[1][1] = y_axis.y;
    frame.axes.m[2][1] = y_axis.z;
    frame.axes.m[0][2] = z_axis.x;
    frame.axes.m[1][2] = z_axis.y;
    frame.axes.m[2][2] = z_axis.z;
    frame.fold_name = fold_name;

#ifndef NDEBUG
    assert(orthonormalError(frame.axes) < 1e-9);
#endif

    return frame;
}
