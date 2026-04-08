#ifndef CAPDAT_GEOMETRY_SYMMETRY_HPP
#define CAPDAT_GEOMETRY_SYMMETRY_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

/**
 * @file geometry_symmetry.hpp
 * @brief Post-parse icosahedral geometry and symmetry utilities for CapDAT.
 *
 * This module is intentionally independent from parsing logic. It provides
 * canonical fold definitions in the VIPERdb-like reference frame plus reusable
 * low-level geometry/rotation helpers that can be consumed by future analysis
 * or orientation workflows after a Capsid has already been reconstructed.
 */
namespace geometry_symmetry {

/**
 * @brief Lightweight 3D Cartesian vector used by geometry/symmetry utilities.
 *
 * Convention: vectors are interpreted as directions or positions from the
 * global origin unless a specific function comment states otherwise.
 */
struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/**
 * @brief Lightweight 3x3 rotation matrix in row-major order.
 *
 * Convention: matrices act on column vectors. Applying a matrix to a vector
 * returns a vector in the same coordinate frame after a proper rotation about
 * the origin.
 */
struct Matrix3 {
    std::array<std::array<double, 3>, 3> m = {{{1.0, 0.0, 0.0},
                                                {0.0, 1.0, 0.0},
                                                {0.0, 0.0, 1.0}}};
};

/**
 * @brief Canonical icosahedral fold definition.
 *
 * `reference_vector` stores the exact documented VIPERdb-frame coordinates.
 * `unit_vector` stores the same direction normalized once at initialization.
 */
struct FoldDefinition {
    int fold_type = 0;
    int fold_index = 0;
    std::string name;
    Vector3 reference_vector;
    Vector3 unit_vector;
    double radius = 0.0;
    int symmetry_order = 0;
};

/**
 * @brief Rotation quality/status classification for direction alignment.
 */
enum class RotationStatus {
    valid,
    identity,
    antiparallel,
    degenerate
};

/**
 * @brief Structured result for a direction-to-direction alignment rotation.
 */
struct RotationDefinition {
    Matrix3 matrix;
    Vector3 axis;
    double angle_radians = 0.0;
    Vector3 source_direction;
    Vector3 target_direction;
    RotationStatus status = RotationStatus::degenerate;
};

/**
 * @brief Direction-vs-IAU classification result.
 */
enum class IauClassification {
    inside,
    boundary,
    outside
};

/**
 * @brief Canonical triangular IAU geometry defined by three named folds.
 *
 * Boundary planes pass through the origin. Their normals are oriented so that
 * points inside the canonical IAU satisfy a consistent half-space inequality.
 */
struct IauDefinition {
    std::string name;
    std::array<std::string, 3> boundary_fold_names;
    std::array<FoldDefinition, 3> boundary_folds;
    std::array<Vector3, 3> boundary_normals;
    std::array<double, 3> boundary_signs;
};

/**
 * @brief Return all canonical folds in stable deterministic order.
 */
const std::vector<FoldDefinition>& canonicalFolds();

/**
 * @brief Look up a canonical fold by exact name (e.g. "2_0", "3_1", "5_0").
 *
 * @throws std::runtime_error when the name does not match a canonical fold.
 */
const FoldDefinition& foldByName(const std::string& fold_name);

/**
 * @brief Return the exact documented reference vector for a canonical fold.
 */
Vector3 foldReferenceVector(const std::string& fold_name);

/**
 * @brief Return the cached unit direction vector for a canonical fold.
 *
 * Callers should use this for angular logic and directional comparisons.
 */
Vector3 foldUnitVector(const std::string& fold_name);

/**
 * @brief Euclidean norm of a 3D vector.
 */
double norm(const Vector3& v);

/**
 * @brief Normalize a vector safely.
 * @param tolerance Zero-length threshold.
 * @throws std::runtime_error for zero/near-zero vectors.
 */
Vector3 normalize(const Vector3& v, double tolerance = 1e-12);

/**
 * @brief Dot product between vectors.
 */
double dot(const Vector3& a, const Vector3& b);

/**
 * @brief Cross product between vectors.
 */
Vector3 cross(const Vector3& a, const Vector3& b);

/**
 * @brief Angular distance in radians between two directions.
 *
 * Inputs are normalized internally, dot products are clamped to [-1, 1], and
 * zero vectors are rejected.
 */
double angleBetween(const Vector3& a, const Vector3& b, double tolerance = 1e-12);

/**
 * @brief Angle in radians between a position vector and a canonical fold axis.
 */
double angleToFoldAxis(const Vector3& position,
                       const std::string& fold_name,
                       double tolerance = 1e-12);

/**
 * @brief Cosine similarity between a position direction and fold axis.
 */
double cosineSimilarityToFoldAxis(const Vector3& position,
                                  const std::string& fold_name,
                                  double tolerance = 1e-12);

/**
 * @brief Return the nearest fold name by smallest angular distance.
 *
 * @param candidate_fold_names List of canonical fold names to search.
 */
std::string nearestFoldName(const Vector3& position,
                            const std::vector<std::string>& candidate_fold_names,
                            double tolerance = 1e-12);

/**
 * @brief Euclidean distance between two Cartesian points.
 */
double euclideanDistance(const Vector3& a, const Vector3& b);

/**
 * @brief Perpendicular distance from point to fold axis line through origin.
 */
double pointToFoldAxisDistance(const Vector3& point,
                               const std::string& fold_name,
                               double tolerance = 1e-12);

/**
 * @brief Radial distance from origin for a Cartesian position vector.
 */
double radialDistance(const Vector3& point);

/**
 * @brief Compute proper rotation mapping `source` direction to `target`.
 *
 * Special handling:
 * - identity status when already aligned within tolerance,
 * - antiparallel status with deterministic orthogonal axis and angle pi,
 * - degenerate status and exception for zero input vectors.
 */
RotationDefinition alignDirectionToDirection(const Vector3& source,
                                             const Vector3& target,
                                             double tolerance = 1e-12);

/**
 * @brief Convenience helper: align named canonical fold to +Z axis.
 */
RotationDefinition alignFoldToPositiveZ(const std::string& fold_name,
                                        double tolerance = 1e-12);

/**
 * @brief Convenience helper: align an arbitrary direction to +Z axis.
 */
RotationDefinition alignDirectionToPositiveZ(const Vector3& direction,
                                             double tolerance = 1e-12);

/**
 * @brief Apply rotation matrix to a direction vector.
 */
Vector3 rotateDirection(const Matrix3& rotation, const Vector3& direction);

/**
 * @brief Apply rotation matrix to a Cartesian point about origin.
 */
Vector3 rotatePoint(const Matrix3& rotation, const Vector3& point);

/**
 * @brief Apply rotation matrix to a collection of Cartesian points.
 */
std::vector<Vector3> rotatePoints(const Matrix3& rotation, const std::vector<Vector3>& points);

/**
 * @brief Return canonical IAU definition object.
 */
const IauDefinition& canonicalIau();

/**
 * @brief Return signed IAU boundary margins for a direction.
 *
 * Positive or near-zero margin values indicate inside/on-boundary for the
 * canonical half-space convention. Input must be nonzero and is normalized.
 */
std::array<double, 3> iauBoundaryMargins(const Vector3& direction, double tolerance = 1e-12);

/**
 * @brief Classify direction relative to canonical IAU as inside/boundary/outside.
 */
IauClassification classifyDirectionInIau(const Vector3& direction,
                                         double tolerance = 1e-9);

/**
 * @brief Check if a direction lies inside canonical IAU (boundary included).
 */
bool isDirectionInsideIau(const Vector3& direction, double tolerance = 1e-9);

/**
 * @brief Check if a direction lies on canonical IAU boundary within tolerance.
 */
bool isDirectionOnIauBoundary(const Vector3& direction, double tolerance = 1e-9);

/**
 * @brief Check rotation-matrix orthogonality and proper-rotation determinant.
 */
bool isProperRotationMatrix(const Matrix3& matrix, double tolerance = 1e-9);

} // namespace geometry_symmetry

#endif // CAPDAT_GEOMETRY_SYMMETRY_HPP
