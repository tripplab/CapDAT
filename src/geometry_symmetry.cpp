#include "geometry_symmetry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace geometry_symmetry {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

bool nearlyEqual(double a, double b, double tolerance) {
    return std::fabs(a - b) <= tolerance;
}

double clampToUnit(double value) {
    return std::max(-1.0, std::min(1.0, value));
}

Vector3 subtract(const Vector3& a, const Vector3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector3 scale(const Vector3& v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}

Vector3 add(const Vector3& a, const Vector3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

/**
 * Canonical fold registry: this is the single source of truth for the VIPERdb
 * standard frame used by this module. These literal vectors must not be edited
 * casually because downstream geometry, IAU boundaries, and orientation logic
 * rely on this exact mapping from canonical fold name to direction.
 */
std::vector<FoldDefinition> buildCanonicalFolds() {
    const std::vector<FoldDefinition> raw = {
        {2, 0, "2_0", {0.000, 0.000, 125.100}, {}, 0.0, 2},
        {2, 1, "2_1", {23.892, 38.658, 125.100}, {}, 0.0, 2},
        {3, 0, "3_0", {47.784, 0.000, 125.100}, {}, 0.0, 3},
        {3, 1, "3_1", {-47.784, 0.000, 125.100}, {}, 0.0, 3},
        {5, 0, "5_0", {0.000, 77.316, 125.100}, {}, 0.0, 5},
    };

    std::vector<FoldDefinition> folds = raw;
    for (FoldDefinition& fold : folds) {
        fold.radius = norm(fold.reference_vector);
        fold.unit_vector = normalize(fold.reference_vector);
    }

    // Orientation anchor invariant: 2_0 is +Z in the canonical frame.
    const FoldDefinition& fold_20 = folds.front();
    if (!nearlyEqual(fold_20.unit_vector.x, 0.0, 1e-12) ||
        !nearlyEqual(fold_20.unit_vector.y, 0.0, 1e-12) ||
        !nearlyEqual(fold_20.unit_vector.z, 1.0, 1e-12)) {
        throw std::runtime_error("Canonical fold invariant violated: 2_0 must align with +Z");
    }

    return folds;
}

const std::vector<FoldDefinition>& foldRegistry() {
    static const std::vector<FoldDefinition> folds = buildCanonicalFolds();
    return folds;
}

const std::unordered_map<std::string, std::size_t>& foldNameIndex() {
    static const std::unordered_map<std::string, std::size_t> index = [] {
        std::unordered_map<std::string, std::size_t> lookup;
        const auto& folds = foldRegistry();
        for (std::size_t i = 0; i < folds.size(); ++i) {
            lookup.emplace(folds[i].name, i);
        }
        return lookup;
    }();
    return index;
}

Matrix3 identityMatrix() {
    return Matrix3{};
}

/**
 * Deterministic fallback axis for antiparallel direction alignment.
 *
 * Rule: pick the canonical basis vector least aligned with source direction,
 * then cross source with that basis. This guarantees a nonzero orthogonal axis
 * while yielding deterministic behavior for reproducible rotations.
 */
Vector3 antiparallelAxis(const Vector3& unit_source) {
    const std::array<Vector3, 3> basis = {{{1.0, 0.0, 0.0},
                                           {0.0, 1.0, 0.0},
                                           {0.0, 0.0, 1.0}}};

    std::size_t best_index = 0;
    double best_abs_dot = std::fabs(dot(unit_source, basis[0]));
    for (std::size_t i = 1; i < basis.size(); ++i) {
        const double candidate = std::fabs(dot(unit_source, basis[i]));
        if (candidate < best_abs_dot) {
            best_abs_dot = candidate;
            best_index = i;
        }
    }

    return normalize(cross(unit_source, basis[best_index]));
}

Matrix3 rotationMatrixFromAxisAngle(const Vector3& unit_axis, double angle_radians) {
    const double c = std::cos(angle_radians);
    const double s = std::sin(angle_radians);
    const double one_minus_c = 1.0 - c;

    const double x = unit_axis.x;
    const double y = unit_axis.y;
    const double z = unit_axis.z;

    Matrix3 rotation;
    rotation.m[0][0] = c + x * x * one_minus_c;
    rotation.m[0][1] = x * y * one_minus_c - z * s;
    rotation.m[0][2] = x * z * one_minus_c + y * s;

    rotation.m[1][0] = y * x * one_minus_c + z * s;
    rotation.m[1][1] = c + y * y * one_minus_c;
    rotation.m[1][2] = y * z * one_minus_c - x * s;

    rotation.m[2][0] = z * x * one_minus_c - y * s;
    rotation.m[2][1] = z * y * one_minus_c + x * s;
    rotation.m[2][2] = c + z * z * one_minus_c;

    return rotation;
}

IauDefinition buildCanonicalIau() {
    IauDefinition iau;
    iau.name = "IAU_2_0_3_0_5_0";
    iau.boundary_fold_names = {"2_0", "3_0", "5_0"};
    iau.boundary_folds = {foldByName("2_0"), foldByName("3_0"), foldByName("5_0")};

    // Edge sequence around the canonical IAU spherical triangle.
    const std::array<std::pair<int, int>, 3> edges = {{{0, 1}, {1, 2}, {2, 0}}};

    // Interior reference direction for sign disambiguation.
    const Vector3 centroid = normalize(add(add(iau.boundary_folds[0].unit_vector,
                                               iau.boundary_folds[1].unit_vector),
                                          iau.boundary_folds[2].unit_vector));

    for (std::size_t i = 0; i < edges.size(); ++i) {
        const Vector3 edge_a = iau.boundary_folds[static_cast<std::size_t>(edges[i].first)].unit_vector;
        const Vector3 edge_b = iau.boundary_folds[static_cast<std::size_t>(edges[i].second)].unit_vector;

        Vector3 normal = normalize(cross(edge_a, edge_b));
        double sign = 1.0;
        if (dot(normal, centroid) < 0.0) {
            normal = scale(normal, -1.0);
            sign = -1.0;
        }

        iau.boundary_normals[i] = normal;
        iau.boundary_signs[i] = sign;
    }

    return iau;
}

} // namespace

const std::vector<FoldDefinition>& canonicalFolds() {
    return foldRegistry();
}

const FoldDefinition& foldByName(const std::string& fold_name) {
    const auto& index = foldNameIndex();
    const auto it = index.find(fold_name);
    if (it == index.end()) {
        throw std::runtime_error("Invalid canonical fold name: '" + fold_name +
                                 "'. Valid names are 2_0, 2_1, 3_0, 3_1, 5_0.");
    }
    return foldRegistry()[it->second];
}

Vector3 foldReferenceVector(const std::string& fold_name) {
    return foldByName(fold_name).reference_vector;
}

Vector3 foldUnitVector(const std::string& fold_name) {
    return foldByName(fold_name).unit_vector;
}

double norm(const Vector3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vector3 normalize(const Vector3& v, double tolerance) {
    const double length = norm(v);
    if (length <= tolerance) {
        throw std::runtime_error("Cannot normalize zero-length vector");
    }
    return {v.x / length, v.y / length, v.z / length};
}

double dot(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 cross(const Vector3& a, const Vector3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double angleBetween(const Vector3& a, const Vector3& b, double tolerance) {
    const Vector3 a_unit = normalize(a, tolerance);
    const Vector3 b_unit = normalize(b, tolerance);
    const double cosine = clampToUnit(dot(a_unit, b_unit));
    return std::acos(cosine);
}

double angleToFoldAxis(const Vector3& position, const std::string& fold_name, double tolerance) {
    return angleBetween(position, foldUnitVector(fold_name), tolerance);
}

double cosineSimilarityToFoldAxis(const Vector3& position,
                                  const std::string& fold_name,
                                  double tolerance) {
    const Vector3 unit_position = normalize(position, tolerance);
    return clampToUnit(dot(unit_position, foldUnitVector(fold_name)));
}

std::string nearestFoldName(const Vector3& position,
                            const std::vector<std::string>& candidate_fold_names,
                            double tolerance) {
    if (candidate_fold_names.empty()) {
        throw std::runtime_error("nearestFoldName requires at least one candidate fold");
    }

    const Vector3 unit_position = normalize(position, tolerance);
    std::string best_name;
    double best_angle = std::numeric_limits<double>::infinity();

    for (const std::string& fold_name : candidate_fold_names) {
        const double angle = angleBetween(unit_position, foldUnitVector(fold_name), tolerance);
        if (angle < best_angle) {
            best_angle = angle;
            best_name = fold_name;
        }
    }

    return best_name;
}

double euclideanDistance(const Vector3& a, const Vector3& b) {
    return norm(subtract(a, b));
}

double pointToFoldAxisDistance(const Vector3& point,
                               const std::string& fold_name,
                               double tolerance) {
    const Vector3 axis = foldUnitVector(fold_name);
    const Vector3 point_cross_axis = cross(point, axis);

    // For a line through the origin with unit direction u, distance is |p x u|.
    (void)tolerance;
    return norm(point_cross_axis);
}

double radialDistance(const Vector3& point) {
    return norm(point);
}

RotationDefinition alignDirectionToDirection(const Vector3& source,
                                             const Vector3& target,
                                             double tolerance) {
    const Vector3 source_unit = normalize(source, tolerance);
    const Vector3 target_unit = normalize(target, tolerance);

    RotationDefinition result;
    result.source_direction = source_unit;
    result.target_direction = target_unit;

    const double cosine = clampToUnit(dot(source_unit, target_unit));

    if (cosine >= 1.0 - tolerance) {
        result.matrix = identityMatrix();
        result.axis = {0.0, 0.0, 1.0};
        result.angle_radians = 0.0;
        result.status = RotationStatus::identity;
        return result;
    }

    if (cosine <= -1.0 + tolerance) {
        const Vector3 axis = antiparallelAxis(source_unit);
        result.axis = axis;
        result.angle_radians = kPi;
        result.matrix = rotationMatrixFromAxisAngle(axis, result.angle_radians);
        result.status = RotationStatus::antiparallel;
        return result;
    }

    const Vector3 axis = normalize(cross(source_unit, target_unit), tolerance);
    const double angle = std::acos(cosine);

    result.axis = axis;
    result.angle_radians = angle;
    result.matrix = rotationMatrixFromAxisAngle(axis, angle);
    result.status = RotationStatus::valid;
    return result;
}

RotationDefinition alignFoldToPositiveZ(const std::string& fold_name, double tolerance) {
    return alignDirectionToDirection(foldUnitVector(fold_name), {0.0, 0.0, 1.0}, tolerance);
}

RotationDefinition alignDirectionToPositiveZ(const Vector3& direction, double tolerance) {
    return alignDirectionToDirection(direction, {0.0, 0.0, 1.0}, tolerance);
}

Vector3 rotateDirection(const Matrix3& rotation, const Vector3& direction) {
    return {rotation.m[0][0] * direction.x + rotation.m[0][1] * direction.y +
                rotation.m[0][2] * direction.z,
            rotation.m[1][0] * direction.x + rotation.m[1][1] * direction.y +
                rotation.m[1][2] * direction.z,
            rotation.m[2][0] * direction.x + rotation.m[2][1] * direction.y +
                rotation.m[2][2] * direction.z};
}

Vector3 rotatePoint(const Matrix3& rotation, const Vector3& point) {
    return rotateDirection(rotation, point);
}

std::vector<Vector3> rotatePoints(const Matrix3& rotation, const std::vector<Vector3>& points) {
    std::vector<Vector3> rotated;
    rotated.reserve(points.size());
    for (const Vector3& point : points) {
        rotated.push_back(rotatePoint(rotation, point));
    }
    return rotated;
}

const IauDefinition& canonicalIau() {
    static const IauDefinition iau = buildCanonicalIau();
    return iau;
}

std::array<double, 3> iauBoundaryMargins(const Vector3& direction, double tolerance) {
    const Vector3 unit_direction = normalize(direction, tolerance);
    const IauDefinition& iau = canonicalIau();

    std::array<double, 3> margins{};
    for (std::size_t i = 0; i < margins.size(); ++i) {
        // With canonicalized normal orientation, inside means margin >= 0.
        margins[i] = dot(iau.boundary_normals[i], unit_direction);
    }
    return margins;
}

IauClassification classifyDirectionInIau(const Vector3& direction, double tolerance) {
    const std::array<double, 3> margins = iauBoundaryMargins(direction, tolerance);

    bool any_boundary = false;
    for (double margin : margins) {
        if (margin < -tolerance) {
            return IauClassification::outside;
        }
        if (std::fabs(margin) <= tolerance) {
            any_boundary = true;
        }
    }

    return any_boundary ? IauClassification::boundary : IauClassification::inside;
}

bool isDirectionInsideIau(const Vector3& direction, double tolerance) {
    const IauClassification cls = classifyDirectionInIau(direction, tolerance);
    return cls == IauClassification::inside || cls == IauClassification::boundary;
}

bool isDirectionOnIauBoundary(const Vector3& direction, double tolerance) {
    return classifyDirectionInIau(direction, tolerance) == IauClassification::boundary;
}

bool isProperRotationMatrix(const Matrix3& matrix, double tolerance) {
    const Vector3 c0 = {matrix.m[0][0], matrix.m[1][0], matrix.m[2][0]};
    const Vector3 c1 = {matrix.m[0][1], matrix.m[1][1], matrix.m[2][1]};
    const Vector3 c2 = {matrix.m[0][2], matrix.m[1][2], matrix.m[2][2]};

    const bool orthonormal = nearlyEqual(norm(c0), 1.0, tolerance) &&
                             nearlyEqual(norm(c1), 1.0, tolerance) &&
                             nearlyEqual(norm(c2), 1.0, tolerance) &&
                             nearlyEqual(dot(c0, c1), 0.0, tolerance) &&
                             nearlyEqual(dot(c0, c2), 0.0, tolerance) &&
                             nearlyEqual(dot(c1, c2), 0.0, tolerance);

    const double determinant = matrix.m[0][0] * (matrix.m[1][1] * matrix.m[2][2] -
                                                  matrix.m[1][2] * matrix.m[2][1]) -
                               matrix.m[0][1] * (matrix.m[1][0] * matrix.m[2][2] -
                                                  matrix.m[1][2] * matrix.m[2][0]) +
                               matrix.m[0][2] * (matrix.m[1][0] * matrix.m[2][1] -
                                                  matrix.m[1][1] * matrix.m[2][0]);

    return orthonormal && nearlyEqual(determinant, 1.0, 1e-6);
}

} // namespace geometry_symmetry
