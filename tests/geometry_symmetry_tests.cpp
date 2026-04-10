#include "geometry_symmetry.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool nearlyEqual(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void assertNear(double actual, double expected, const std::string& label, double eps = 1e-9) {
    if (!nearlyEqual(actual, expected, eps)) {
        throw std::runtime_error(label + " mismatch: expected " + std::to_string(expected) +
                                 ", got " + std::to_string(actual));
    }
}

void testFoldDefinitions() {
    const auto& folds = geometry_symmetry::canonicalFolds();
    assertTrue(folds.size() == 5, "Expected exactly five canonical folds");

    const std::string fold_20_name = "2_0";
    const auto& fold_20 = geometry_symmetry::foldByName(fold_20_name);
    assertTrue(fold_20.fold_type == 2, "2_0 fold_type mismatch");
    assertTrue(fold_20.fold_index == 0, "2_0 fold_index mismatch");
    assertNear(fold_20.reference_vector.x, 0.0, "2_0 ref x");
    assertNear(fold_20.reference_vector.y, 0.0, "2_0 ref y");
    assertNear(fold_20.reference_vector.z, 125.1, "2_0 ref z");
    assertNear(geometry_symmetry::norm(fold_20.unit_vector), 1.0, "2_0 unit norm", 1e-12);
    assertNear(fold_20.unit_vector.x, 0.0, "2_0 unit x", 1e-12);
    assertNear(fold_20.unit_vector.y, 0.0, "2_0 unit y", 1e-12);
    assertNear(fold_20.unit_vector.z, 1.0, "2_0 unit z", 1e-12);

    const std::string fold_21_name = "2_1";
    const auto& fold_21 = geometry_symmetry::foldByName(fold_21_name);
    assertNear(fold_21.reference_vector.x, 23.892, "2_1 ref x");
    assertNear(fold_21.reference_vector.y, 38.658, "2_1 ref y");
    assertNear(fold_21.reference_vector.z, 125.1, "2_1 ref z");
    assertNear(geometry_symmetry::norm(fold_21.unit_vector), 1.0, "2_1 unit norm", 1e-12);

    const std::string fold_30_name = "3_0";
    const auto& fold_30 = geometry_symmetry::foldByName(fold_30_name);
    assertNear(fold_30.reference_vector.x, 47.784, "3_0 ref x");
    assertNear(fold_30.reference_vector.y, 0.0, "3_0 ref y");
    assertNear(fold_30.reference_vector.z, 125.1, "3_0 ref z");

    const std::string fold_31_name = "3_1";
    const auto& fold_31 = geometry_symmetry::foldByName(fold_31_name);
    assertNear(fold_31.reference_vector.x, -47.784, "3_1 ref x");
    assertNear(fold_31.reference_vector.y, 0.0, "3_1 ref y");
    assertNear(fold_31.reference_vector.z, 125.1, "3_1 ref z");

    const std::string fold_50_name = "5_0";
    const auto& fold_50 = geometry_symmetry::foldByName(fold_50_name);
    assertNear(fold_50.reference_vector.x, 0.0, "5_0 ref x");
    assertNear(fold_50.reference_vector.y, 77.316, "5_0 ref y");
    assertNear(fold_50.reference_vector.z, 125.1, "5_0 ref z");
}

void testLookup() {
    const std::string fold_name = "3_1";
    const auto& fold_a = geometry_symmetry::foldByName(fold_name);
    const auto& fold_b = geometry_symmetry::foldByName(fold_name);
    assertTrue(&fold_a == &fold_b, "Fold lookup should be deterministic and stable");

    bool threw = false;
    try {
        (void)geometry_symmetry::foldByName("3_9");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assertTrue(threw, "Invalid fold lookup should throw");
}

void testAngleDistance() {
    const geometry_symmetry::Vector3 v = {1.0, 2.0, 3.0};
    const geometry_symmetry::Vector3 neg_v = {-1.0, -2.0, -3.0};

    assertNear(geometry_symmetry::angleBetween(v, v), 0.0, "angle self", 1e-12);
    assertNear(geometry_symmetry::angleBetween(v, neg_v), 3.14159265358979323846, "angle opposite", 1e-9);

    const geometry_symmetry::Vector3 on_axis = {0.0, 0.0, 20.0};
    assertNear(geometry_symmetry::pointToFoldAxisDistance(on_axis, "2_0"), 0.0, "point-to-axis on-axis", 1e-12);

    const geometry_symmetry::Vector3 p = {3.0, 4.0, 12.0};
    assertNear(geometry_symmetry::radialDistance(p), 13.0, "radial distance", 1e-12);

    const std::string nearest =
        geometry_symmetry::nearestFoldName({0.0, 0.0, 1.0}, {"2_1", "2_0", "5_0"});
    assertTrue(nearest == "2_0", "nearestFoldName should find 2_0 for +Z direction");
}

void testRotations() {
    const auto identity = geometry_symmetry::alignDirectionToDirection({0.0, 0.0, 1.0}, {0.0, 0.0, 1.0});
    assertTrue(identity.status == geometry_symmetry::RotationStatus::identity,
               "identity alignment should return identity status");
    assertTrue(geometry_symmetry::isProperRotationMatrix(identity.matrix),
               "identity matrix should be proper rotation");

    const auto fold_to_z = geometry_symmetry::alignFoldToPositiveZ("5_0");
    const auto rotated = geometry_symmetry::rotateDirection(
        fold_to_z.matrix, geometry_symmetry::foldUnitVector("5_0"));
    assertNear(rotated.x, 0.0, "fold->z x", 1e-8);
    assertNear(rotated.y, 0.0, "fold->z y", 1e-8);
    assertNear(rotated.z, 1.0, "fold->z z", 1e-8);
    assertTrue(geometry_symmetry::isProperRotationMatrix(fold_to_z.matrix),
               "fold->z matrix should be proper rotation");

    const auto antiparallel =
        geometry_symmetry::alignDirectionToDirection({0.0, 0.0, 1.0}, {0.0, 0.0, -1.0});
    assertTrue(antiparallel.status == geometry_symmetry::RotationStatus::antiparallel,
               "antiparallel status expected");
    assertNear(antiparallel.angle_radians, 3.14159265358979323846, "antiparallel angle", 1e-9);
    assertTrue(geometry_symmetry::isProperRotationMatrix(antiparallel.matrix),
               "antiparallel matrix should be proper rotation");

    bool threw = false;
    try {
        (void)geometry_symmetry::alignDirectionToDirection({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assertTrue(threw, "zero-vector source should throw");
}


void testFoldTypeIndexLookup() {
    assertTrue(geometry_symmetry::foldByTypeIndex(2, 0).name == "2_0", "(2,0) should resolve to 2_0");
    assertTrue(geometry_symmetry::foldByTypeIndex(2, 1).name == "2_1", "(2,1) should resolve to 2_1");
    assertTrue(geometry_symmetry::foldByTypeIndex(3, 0).name == "3_0", "(3,0) should resolve to 3_0");
    assertTrue(geometry_symmetry::foldByTypeIndex(3, 1).name == "3_1", "(3,1) should resolve to 3_1");
    assertTrue(geometry_symmetry::foldByTypeIndex(5, 0).name == "5_0", "(5,0) should resolve to 5_0");

    bool invalid_type = false;
    try {
        (void)geometry_symmetry::foldByTypeIndex(4, 0);
    } catch (const std::runtime_error&) {
        invalid_type = true;
    }
    assertTrue(invalid_type, "invalid fold type should throw");

    bool invalid_index = false;
    try {
        (void)geometry_symmetry::foldByTypeIndex(3, 9);
    } catch (const std::runtime_error&) {
        invalid_index = true;
    }
    assertTrue(invalid_index, "invalid fold index should throw");
}

void testIauHelpers() {
    const auto& iau = geometry_symmetry::canonicalIau();
    assertTrue(iau.boundary_fold_names[0] == "2_0", "IAU first fold should be 2_0");
    assertTrue(iau.boundary_fold_names[1] == "3_0", "IAU second fold should be 3_0");
    assertTrue(iau.boundary_fold_names[2] == "5_0", "IAU third fold should be 5_0");

    const auto cls_20 = geometry_symmetry::classifyDirectionInIau(geometry_symmetry::foldUnitVector("2_0"));
    assertTrue(cls_20 != geometry_symmetry::IauClassification::outside,
               "IAU defining fold 2_0 should not classify as outside");

    const geometry_symmetry::Vector3 interior = geometry_symmetry::normalize(
        {geometry_symmetry::foldUnitVector("2_0").x + geometry_symmetry::foldUnitVector("3_0").x +
             geometry_symmetry::foldUnitVector("5_0").x,
         geometry_symmetry::foldUnitVector("2_0").y + geometry_symmetry::foldUnitVector("3_0").y +
             geometry_symmetry::foldUnitVector("5_0").y,
         geometry_symmetry::foldUnitVector("2_0").z + geometry_symmetry::foldUnitVector("3_0").z +
             geometry_symmetry::foldUnitVector("5_0").z});

    assertTrue(geometry_symmetry::classifyDirectionInIau(interior) == geometry_symmetry::IauClassification::inside,
               "interior direction should classify as inside");

    assertTrue(geometry_symmetry::classifyDirectionInIau({-1.0, -1.0, -1.0}) ==
                   geometry_symmetry::IauClassification::outside,
               "opposite octant should classify as outside");

    const auto boundary_class =
        geometry_symmetry::classifyDirectionInIau(geometry_symmetry::foldUnitVector("3_0"), 1e-6);
    assertTrue(boundary_class == geometry_symmetry::IauClassification::boundary ||
                   boundary_class == geometry_symmetry::IauClassification::inside,
               "boundary tolerance should classify edge-adjacent direction consistently");
}

} // namespace

int main() {
    try {
        testFoldDefinitions();
        testLookup();
        testFoldTypeIndexLookup();
        testAngleDistance();
        testRotations();
        testIauHelpers();
        std::cout << "All geometry symmetry tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Geometry symmetry test failure: " << e.what() << '\n';
        return 1;
    }
}
