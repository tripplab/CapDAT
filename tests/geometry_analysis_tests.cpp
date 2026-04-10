#include "geometry_analysis.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

bool near(double a, double b, double eps = 1e-8) { return std::fabs(a - b) <= eps; }

void assertTrue(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

Capsid makeSimpleCapsid() {
    Capsid capsid("constructed");
    Chain chain(1, 'A');
    Residue residue("ALA", 1, ' ', 'A', 1);
    residue.addAtom(Atom(1, "CA", ' ', "ALA", 'A', 1, ' ', 1.0, 0.0, 0.0, 1.0, 20.0, "C", "", false));
    residue.addAtom(Atom(2, "CB", ' ', "ALA", 'A', 1, ' ', 0.0, 1.0, 0.0, 1.0, 20.0, "C", "", false));
    chain.addResidue(residue);
    capsid.addChain(chain);
    capsid.finalizeCounts();
    return capsid;
}

ParserConfig makeParserConfig() {
    ParserConfig config;
    config.include_hetatm = false;
    config.strict_mode = false;
    config.keep_altloc_all = true;
    config.ignore_blank_chain_id = false;
    config.verbose_warnings = false;
    config.protein_only = true;
    return config;
}

void testStage1IdentityFor2_0() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 12.0;

    const auto result = prepareGeometryAnalysisStage1(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1 should succeed for 2_0");
    assertTrue(result.resolved_fold_name == "2_0", "resolved fold should be 2_0");
    assertTrue(result.used_identity_rotation, "2_0 -> +Z should be identity");
    assertTrue(!result.coordinates_modified_in_place, "identity transform should not alter coordinates");

    const auto& state = capsid.orientationState();
    assertTrue(state.reoriented_in_place, "orientation state should mark reoriented workflow");
    assertTrue(state.requested_target_axis == 'z', "target axis should be z");
    assertTrue(state.already_aligned_identity, "identity should be recorded in orientation state");
}

void testStage1NonIdentityFor3_0() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 3;
    config.fold_index = 0;
    config.cylinder_radius = 12.0;

    const double x_before = capsid.chains()[0].residues()[0].atoms()[0].x();
    const double z_before = capsid.chains()[0].residues()[0].atoms()[0].z();

    const auto result = prepareGeometryAnalysisStage1(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1 should succeed for 3_0");
    assertTrue(result.resolved_fold_name == "3_0", "resolved fold should be 3_0");
    assertTrue(!result.used_identity_rotation, "3_0 -> +Z should not be identity");
    assertTrue(result.coordinates_modified_in_place, "non-identity transform should modify coordinates");

    const double x_after = capsid.chains()[0].residues()[0].atoms()[0].x();
    const double z_after = capsid.chains()[0].residues()[0].atoms()[0].z();
    assertTrue(!near(x_before, x_after) || !near(z_before, z_after), "atom coordinates should change");

    const auto& state = capsid.orientationState();
    assertTrue(state.reoriented_in_place, "orientation state should be set");
    assertTrue(!state.already_aligned_identity, "non-identity should be recorded");
}

void testInvalidFoldIndex() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 3;
    config.fold_index = 9;
    config.cylinder_radius = 12.0;

    bool threw = false;
    try {
        (void)prepareGeometryAnalysisStage1(capsid, config, makeParserConfig(), nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("Invalid geometry fold index for fold type 3") != std::string::npos;
    }
    assertTrue(threw, "invalid fold index should fail clearly");
}

void testInvalidCylinderRadius() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 0.0;

    bool threw = false;
    try {
        (void)prepareGeometryAnalysisStage1(capsid, config, makeParserConfig(), nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("Geometry cylinder radius must be > 0") != std::string::npos;
    }
    assertTrue(threw, "non-positive cylinder radius should fail validation");
}

void testOptionalExport() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 5;
    config.fold_index = 0;
    config.cylinder_radius = 12.0;
    config.export_rotated_capsid = true;
    config.output_prefix = "geometry_test_output";

    const auto result = prepareGeometryAnalysisStage1(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1 export run should succeed");
    assertTrue(!result.export_path.empty(), "export path should be reported");
    assertTrue(std::filesystem::exists(result.export_path), "rotated capsid export should exist");

    std::filesystem::remove(result.export_path);
}

} // namespace

int main() {
    try {
        testStage1IdentityFor2_0();
        testStage1NonIdentityFor3_0();
        testInvalidFoldIndex();
        testInvalidCylinderRadius();
        testOptionalExport();
        std::cout << "All geometry analysis Stage 1 tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Geometry analysis test failure: " << e.what() << '\n';
        return 1;
    }
}
