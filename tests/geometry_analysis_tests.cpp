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
    residue.addAtom(Atom(1, "CA", ' ', "ALA", 'A', 1, ' ', 1.0, 0.0, 1.0, 1.0, 20.0, "C", "", false));
    residue.addAtom(Atom(2, "CB", ' ', "ALA", 'A', 1, ' ', 0.0, 1.0, 1.0, 1.0, 20.0, "C", "", false));
    chain.addResidue(residue);
    capsid.addChain(chain);
    capsid.finalizeCounts();
    return capsid;
}

Capsid makeStage2Capsid() {
    Capsid capsid("stage2_constructed");
    Chain chain(1, 'A');
    Residue residue("GLY", 10, ' ', 'A', 1);
    residue.addAtom(Atom(1, "N", ' ', "GLY", 'A', 10, ' ', 1.0, 1.0, 1.0, 1.0, 20.0, "N", "", false));  // selected
    residue.addAtom(Atom(2, "CA", ' ', "GLY", 'A', 10, ' ', 2.0, 0.0, 2.0, 1.0, 20.0, "C", "", false)); // selected (edge)
    residue.addAtom(Atom(3, "C", ' ', "GLY", 'A', 10, ' ', 0.5, 0.5, 0.0, 1.0, 20.0, "C", "", false));  // reject z==0
    residue.addAtom(Atom(4, "O", ' ', "GLY", 'A', 10, ' ', 0.2, 0.2, -1.0, 1.0, 20.0, "O", "", false)); // reject z<0
    residue.addAtom(Atom(5, "CB", ' ', "GLY", 'A', 10, ' ', 2.1, 0.0, 3.0, 1.0, 20.0, "C", "", false)); // reject radial
    chain.addResidue(residue);
    capsid.addChain(chain);
    capsid.finalizeCounts();

    Capsid::OrientationState state;
    state.in_original_parsed_frame = false;
    state.reoriented_in_place = true;
    state.source_mode = Capsid::OrientationSourceMode::fold;
    state.source_description = "2_0";
    state.requested_target_axis = 'z';
    state.target_direction = {0.0, 0.0, 1.0};
    capsid.setOrientationState(state);

    return capsid;
}

Capsid makeStage2CapsidWithUnknownElement() {
    Capsid capsid("stage2_unknown_element");
    Chain chain(1, 'A');
    Residue residue("GLY", 11, ' ', 'A', 1);
    residue.addAtom(Atom(1, "N", ' ', "GLY", 'A', 11, ' ', 0.5, 0.5, 1.0, 1.0, 20.0, "N", "", false));
    residue.addAtom(Atom(2, "X1", ' ', "GLY", 'A', 11, ' ', 1.0, 0.0, 1.5, 1.0, 20.0, " Xx ", "", false));
    chain.addResidue(residue);
    capsid.addChain(chain);
    capsid.finalizeCounts();

    Capsid::OrientationState state;
    state.in_original_parsed_frame = false;
    state.reoriented_in_place = true;
    state.source_mode = Capsid::OrientationSourceMode::fold;
    state.source_description = "2_0";
    state.requested_target_axis = 'z';
    state.target_direction = {0.0, 0.0, 1.0};
    capsid.setOrientationState(state);

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

void testCylinderClassifier() {
    const auto inside = classifyPatchCylinder({1.0, 1.0, 0.1}, 2.0);
    assertTrue(inside.selected, "z>0 and radial<radius should be selected");

    const auto zeroZ = classifyPatchCylinder({0.0, 0.0, 0.0}, 2.0);
    assertTrue(!zeroZ.selected, "z==0 should be excluded");

    const auto negativeZ = classifyPatchCylinder({0.0, 0.0, -0.1}, 2.0);
    assertTrue(!negativeZ.selected, "z<0 should be excluded");

    const auto onEdge = classifyPatchCylinder({2.0, 0.0, 1.0}, 2.0);
    assertTrue(onEdge.selected, "r==radius with z>0 should be selected");

    const auto outside = classifyPatchCylinder({2.01, 0.0, 1.0}, 2.0);
    assertTrue(!outside.selected, "r>radius should be excluded");
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

void testStage2PatchSelectionAndTraceability() {
    Capsid capsid = makeStage2Capsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "stage2_patch_traceability";

    GeometryPreparationResult stage1;
    stage1.success = true;

    const auto result = runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);

    assertTrue(result.success, "Stage 2 should succeed");
    assertTrue(result.total_atoms_examined == 5, "all atoms should be examined");
    assertTrue(result.selected_atoms_count == 2, "selected count should match rule");
    assertTrue(result.rejected_non_positive_z_count == 2, "z rejection count should match");
    assertTrue(result.rejected_outside_radius_count == 1, "radius rejection count should match");
    assertTrue(result.patch_atoms.size() == result.selected_atom_refs.size(), "patch and refs sizes must match");
    assertTrue(result.patch_atoms[0].original_atom == result.selected_atom_refs[0], "first ref should match");
    assertTrue(result.patch_atoms[1].original_atom == result.selected_atom_refs[1], "second ref should match");
    assertTrue(result.patch_atoms[0].in_positive_z, "selected atom should record positive z");
    assertTrue(result.patch_atoms[0].in_cylinder_radius, "selected atom should record radial inclusion");

    assertTrue(std::filesystem::exists(result.export_path), "Stage 2 export should exist");
    std::filesystem::remove(result.export_path);
}

void testStage2MinAtomsFailure() {
    Capsid capsid = makeStage2Capsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 3;

    GeometryPreparationResult stage1;
    stage1.success = true;

    bool threw = false;
    try {
        (void)runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("below min_atoms_in_patch=3") != std::string::npos;
    }
    assertTrue(threw, "min_atoms_in_patch should be enforced with clear error");
}

void testStage2RequiresStage1() {
    Capsid capsid = makeStage2Capsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 1;

    GeometryPreparationResult stage1;
    stage1.success = false;

    bool threw = false;
    try {
        (void)runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("cannot run before successful Stage 1") != std::string::npos;
    }
    assertTrue(threw, "Stage 2 should require successful Stage 1");
}

void testVdwRadiusLookupAndNormalization() {
    assertTrue(near(vdwRadius("C"), 1.70), "C vdW radius should be 1.70");
    assertTrue(near(vdwRadius("N"), 1.55), "N vdW radius should be 1.55");
    assertTrue(near(vdwRadius("O"), 1.52), "O vdW radius should be 1.52");
    assertTrue(near(vdwRadius("S"), 1.80), "S vdW radius should be 1.80");
    assertTrue(near(vdwRadius("P"), 1.80), "P vdW radius should be 1.80");
    assertTrue(near(vdwRadius(normalizeElementSymbol("c")), 1.70), "lowercase C should normalize");
    assertTrue(near(vdwRadius(normalizeElementSymbol(" se ")), 1.90), "trimmed se should normalize to SE");
    assertTrue(near(vdwRadius(normalizeElementSymbol("XX")), 1.70), "unknown should fallback to 1.70");
}

void testPatchAtomBuilderNormalization() {
    const Atom atom(99, "CA", ' ', "ALA", 'A', 1, ' ', 1.2, -0.4, 3.1, 1.0, 20.0, " c ", "", false);
    const geometry_symmetry::Vector3 rotated{2.5, -1.0, 4.0};
    const CylinderMembership membership = classifyPatchCylinder(rotated, 3.0);

    const PatchAtom patch_atom = makePatchAtom(atom, rotated, membership);
    assertTrue(near(patch_atom.position.x, rotated.x), "builder should preserve rotated x");
    assertTrue(near(patch_atom.position.y, rotated.y), "builder should preserve rotated y");
    assertTrue(near(patch_atom.position.z, rotated.z), "builder should preserve rotated z");
    assertTrue(patch_atom.element == "C", "builder should normalize element");
    assertTrue(near(patch_atom.vdw_radius, 1.70), "builder should assign vdW radius");
    assertTrue(patch_atom.original_atom == &atom, "builder should preserve original atom reference");
    assertTrue(near(patch_atom.radial_xy, membership.radial_xy), "builder should preserve radial xy");
    assertTrue(patch_atom.membership.selected == membership.selected, "builder should preserve membership facts");
}

void testStage3FailsWithoutSuccessfulStage2() {
    GeometryPatchSelectionResult stage2_result;
    stage2_result.success = false;
    bool threw = false;
    try {
        (void)runGeometryAnalysisStage3PatchNormalization(stage2_result, nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("cannot run before successful Stage 2") != std::string::npos;
    }
    assertTrue(threw, "Stage 3 should require successful Stage 2");
}

void testStage2ToStage3Integration() {
    Capsid capsid = makeStage2Capsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "stage3_integration";

    GeometryPreparationResult stage1;
    stage1.success = true;

    const auto stage2 = runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);
    const auto stage3 = runGeometryAnalysisStage3PatchNormalization(stage2, nullptr);

    assertTrue(stage3.success, "Stage 3 should succeed from valid Stage 2 result");
    assertTrue(stage3.analytical_patch.atom_count == stage2.selected_atom_refs.size(),
               "Stage 3 atom count should match selected refs count");
    assertTrue(stage3.analytical_patch.atoms.size() == stage3.analytical_patch.original_atom_refs.size(),
               "Stage 3 atom/reference counts should match");
    for (const PatchAtom& atom : stage3.analytical_patch.atoms) {
        assertTrue(atom.original_atom != nullptr, "each Stage 3 PatchAtom should keep original reference");
        assertTrue(atom.vdw_radius > 0.0, "each Stage 3 PatchAtom should have positive vdW radius");
    }

    assertTrue(std::filesystem::exists(stage2.export_path), "Stage 2 canonical patch export should exist");
    std::filesystem::remove(stage2.export_path);
}

void testStage3FallbackRadiusIntegration() {
    Capsid capsid = makeStage2CapsidWithUnknownElement();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "stage3_fallback";

    GeometryPreparationResult stage1;
    stage1.success = true;

    const auto stage2 = runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);
    const auto stage3 = runGeometryAnalysisStage3PatchNormalization(stage2, nullptr);
    assertTrue(stage3.success, "Stage 3 should succeed with unknown element via fallback radius");
    assertTrue(stage3.analytical_patch.fallback_vdw_radius_count == 1, "unknown element should use fallback count");
    assertTrue(stage3.analytical_patch.atoms[1].element == "XX", "unknown element should still be normalized and kept");
    assertTrue(near(stage3.analytical_patch.atoms[1].vdw_radius, 1.70), "unknown element should use 1.70 fallback");

    std::filesystem::remove(stage2.export_path);
}

void testPatchExportUsesOriginalSubsetNonRegression() {
    Capsid capsid = makeStage2Capsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "stage3_export_nonreg";

    GeometryPreparationResult stage1;
    stage1.success = true;
    const auto stage2 = runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);
    const auto stage3 = runGeometryAnalysisStage3PatchNormalization(stage2, nullptr);

    assertTrue(stage3.success, "Stage 3 should succeed for export non-regression");
    assertTrue(stage3.analytical_patch.export_path == stage2.export_path,
               "Stage 3 should preserve canonical Stage 2 export path");
    assertTrue(stage2.selected_atom_refs.size() == 2, "Stage 2 should provide expected original subset");
    assertTrue(stage2.selected_atom_refs[0]->serial() == 1, "first selected original atom serial should be preserved");
    assertTrue(stage2.selected_atom_refs[1]->serial() == 2, "second selected original atom serial should be preserved");

    std::filesystem::remove(stage2.export_path);
}

void testEndToEndStage1Stage2Integration() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "stage2_integration";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1 + Stage 2 integration should succeed");
    assertTrue(result.preparation.success, "Stage 1 should succeed in integration");
    assertTrue(result.stage2_patch.success, "Stage 2 should succeed in integration");
    assertTrue(result.stage3_patch.success, "Stage 3 should succeed in integration");
    assertTrue(result.stage2_patch.selected_atoms_count == result.stage2_patch.patch_atoms.size(),
               "selected count should match patch atoms size");
    assertTrue(result.stage3_patch.analytical_patch.atom_count == result.stage2_patch.selected_atom_refs.size(),
               "Stage 3 analytical patch count should match selected refs");
    assertTrue(std::filesystem::exists(result.stage2_patch.export_path), "integration should export patch file");

    std::filesystem::remove(result.stage2_patch.export_path);
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
}

} // namespace

int main() {
    try {
        testCylinderClassifier();
        testVdwRadiusLookupAndNormalization();
        testPatchAtomBuilderNormalization();
        testStage1IdentityFor2_0();
        testStage1NonIdentityFor3_0();
        testStage2PatchSelectionAndTraceability();
        testStage2MinAtomsFailure();
        testStage2RequiresStage1();
        testStage3FailsWithoutSuccessfulStage2();
        testStage2ToStage3Integration();
        testStage3FallbackRadiusIntegration();
        testPatchExportUsesOriginalSubsetNonRegression();
        testEndToEndStage1Stage2Integration();
        std::cout << "All geometry analysis Stage 1/2/3 tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Geometry analysis test failure: " << e.what() << '\n';
        return 1;
    }
}
