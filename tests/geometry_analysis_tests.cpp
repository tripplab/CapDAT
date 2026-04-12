#include "geometry_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

bool near(double a, double b, double eps = 1e-8) { return std::fabs(a - b) <= eps; }
std::size_t nodeIndex(std::size_t i, std::size_t j, std::size_t nx) { return (j * nx) + i; }

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
    residue.addAtom(Atom(1, "N", ' ', "GLY", 'A', 11, ' ', 0.5, 0.5, 1.0, 1.0, 20.0, "", "", false));
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

GeometryStage4RawSheetResult makeSyntheticStage4ResultForStage5() {
    GeometryStage4RawSheetResult stage4;
    stage4.success = true;
    stage4.grid = buildStage4RegularGrid(2.0, 1.0);
    stage4.node_count = stage4.grid.nx * stage4.grid.ny;
    stage4.z_outer_raw.assign(stage4.node_count, std::numeric_limits<double>::quiet_NaN());
    stage4.z_inner_raw.assign(stage4.node_count, std::numeric_limits<double>::quiet_NaN());
    stage4.inside_disk_mask.assign(stage4.node_count, 0);
    stage4.valid_mask.assign(stage4.node_count, 0);
    stage4.outer_contact_serial_numbers.assign(stage4.node_count, 0);
    stage4.inner_contact_serial_numbers.assign(stage4.node_count, 0);
    stage4.candidate_patch_atom_counts.assign(stage4.node_count, 0);

    const double radius2 = 4.0;
    for (std::size_t j = 0; j < stage4.grid.ny; ++j) {
        for (std::size_t i = 0; i < stage4.grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, stage4.grid.nx);
            const double x = stage4.grid.x_values[i];
            const double y = stage4.grid.y_values[j];
            if ((x * x) + (y * y) <= radius2 + 1e-12) {
                stage4.inside_disk_mask[idx] = 1;
                ++stage4.inside_disk_count;
            }
        }
    }

    auto markValid = [&](std::size_t i, std::size_t j, double z_outer, double z_inner, int outer_serial,
                         int inner_serial) {
        const std::size_t idx = nodeIndex(i, j, stage4.grid.nx);
        stage4.valid_mask[idx] = 1;
        stage4.z_outer_raw[idx] = z_outer;
        stage4.z_inner_raw[idx] = z_inner;
        stage4.outer_contact_serial_numbers[idx] = outer_serial;
        stage4.inner_contact_serial_numbers[idx] = inner_serial;
        ++stage4.valid_node_count;
    };

    markValid(2, 2, 8.0, 4.0, 101, 201);
    markValid(1, 2, 8.1, 4.1, 102, 202);
    markValid(3, 2, 8.2, 4.2, 103, 203);
    markValid(2, 1, 8.3, 4.3, 104, 204);
    markValid(2, 3, 8.4, 4.4, 105, 205);
    markValid(1, 1, 8.5, 4.5, 101, 201);
    markValid(3, 3, 8.6, 4.6, 103, 203);

    stage4.invalid_node_count = stage4.inside_disk_count - stage4.valid_node_count;
    return stage4;
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

void testLineSphereIntersectionHelper() {
    PatchAtom atom;
    atom.position = {0.0, 0.0, 5.0};
    atom.vdw_radius = 2.0;

    const auto miss = intersectVerticalLineWithSphere(3.0, 0.0, atom);
    assertTrue(!miss.intersects, "d2 > r2 should not intersect");

    const auto tangent = intersectVerticalLineWithSphere(2.0, 0.0, atom);
    assertTrue(tangent.intersects, "tangent should intersect");
    assertTrue(near(tangent.z_low, 5.0), "tangent z_low should match center z");
    assertTrue(near(tangent.z_high, 5.0), "tangent z_high should match center z");

    const auto proper = intersectVerticalLineWithSphere(1.0, 0.0, atom);
    assertTrue(proper.intersects, "d2 < r2 should intersect");
    assertTrue(near(proper.z_low, 5.0 - std::sqrt(3.0)), "proper z_low should be center-delta");
    assertTrue(near(proper.z_high, 5.0 + std::sqrt(3.0)), "proper z_high should be center+delta");

    PatchAtom shifted;
    shifted.position = {0.0, 0.0, 1000.0};
    shifted.vdw_radius = 2.0;
    const auto shifted_hit = intersectVerticalLineWithSphere(0.0, 0.0, shifted);
    assertTrue(shifted_hit.intersects, "large positive z sphere should intersect identically");
    assertTrue(near(shifted_hit.z_low, 998.0), "shifted sphere z_low should be geometric only");
    assertTrue(near(shifted_hit.z_high, 1002.0), "shifted sphere z_high should be geometric only");

    PatchAtom low;
    low.position = {0.0, 0.0, -1000.0};
    low.vdw_radius = 2.0;
    const auto low_hit = intersectVerticalLineWithSphere(2.0, 0.0, low);
    assertTrue(low_hit.intersects, "tangent at large negative z should intersect");
    assertTrue(near(low_hit.z_low, -1000.0), "negative-z tangent z_low should match center");
    assertTrue(near(low_hit.z_high, -1000.0), "negative-z tangent z_high should match center");

    const auto near_miss = intersectVerticalLineWithSphere(2.000001, 0.0, atom, 1e-12);
    assertTrue(!near_miss.intersects, "just outside radius should miss with strict tolerance");
}

void testSingleNodeFirstContact() {
    std::vector<PatchAtom> atoms;

    PatchAtom a;
    a.position = {0.0, 0.0, 5.0};
    a.vdw_radius = 2.0;
    atoms.push_back(a);

    const auto single = detectRawFirstContactAtNode(0.0, 0.0, atoms);
    assertTrue(single.valid, "one intersecting sphere should be valid");
    assertTrue(near(single.z_outer_raw, 5.0), "single sphere outer should be atom z");
    assertTrue(near(single.z_inner_raw, 5.0), "single sphere inner should be atom z");
    assertTrue(single.candidate_patch_atom_count == 1, "single intersecting atom should be counted as candidate");

    PatchAtom b;
    b.position = {0.0, 0.0, 8.0};
    b.vdw_radius = 4.0;
    atoms.push_back(b);

    const auto multi = detectRawFirstContactAtNode(0.0, 0.0, atoms);
    assertTrue(multi.valid, "multiple intersecting spheres should be valid");
    assertTrue(near(multi.z_outer_raw, 8.0), "outer should be maximum atom z");
    assertTrue(near(multi.z_inner_raw, 5.0), "inner should be minimum atom z");
    assertTrue(multi.candidate_patch_atom_count == 2, "both intersecting atoms should be counted as candidates");

    const auto none = detectRawFirstContactAtNode(20.0, 20.0, atoms);
    assertTrue(!none.valid, "no intersection should be invalid");

    std::vector<PatchAtom> tie_atoms;
    PatchAtom t0;
    t0.position = {0.0, 0.0, 5.0};
    t0.vdw_radius = 2.0;
    tie_atoms.push_back(t0);
    PatchAtom t1;
    t1.position = {0.0, 0.0, 5.0};
    t1.vdw_radius = 2.0;
    tie_atoms.push_back(t1);
    const auto tie = detectRawFirstContactAtNode(0.0, 0.0, tie_atoms);
    assertTrue(tie.valid, "tie case should remain valid");
    assertTrue(tie.outer_patch_atom_index == 0, "tie should prefer earlier patch atom for outer");
    assertTrue(tie.inner_patch_atom_index == 0, "tie should prefer earlier patch atom for inner");

    std::vector<PatchAtom> disjoint;
    PatchAtom d0;
    d0.position = {0.0, 0.0, 91.0}; // [90,92]
    d0.vdw_radius = 1.0;
    disjoint.push_back(d0);
    PatchAtom d1;
    d1.position = {0.0, 0.0, 101.0}; // [100,102]
    d1.vdw_radius = 1.0;
    disjoint.push_back(d1);
    const auto disjoint_case = detectRawFirstContactAtNode(0.0, 0.0, disjoint);
    assertTrue(disjoint_case.valid, "non-overlapping intervals should remain valid with max/min envelope rule");
    assertTrue(near(disjoint_case.z_outer_raw, 101.0), "disjoint outer should be max atom z");
    assertTrue(near(disjoint_case.z_inner_raw, 91.0), "disjoint inner should be min atom z");
}

void testGridConstruction() {
    const auto grid = buildStage4RegularGrid(2.0, 1.0);
    assertTrue(grid.nx == 5 && grid.ny == 5, "grid dimensions should include both -R and +R");
    assertTrue(near(grid.x_values[0], -2.0) && near(grid.x_values[4], 2.0), "x bounds should be correct");
    assertTrue(near(grid.y_values[0], -2.0) && near(grid.y_values[4], 2.0), "y bounds should be correct");

    const bool center_inside = ((grid.x_values[2] * grid.x_values[2]) + (grid.y_values[2] * grid.y_values[2])) <= 4.0;
    const bool corner_inside = ((grid.x_values[0] * grid.x_values[0]) + (grid.y_values[0] * grid.y_values[0])) <= 4.0;
    assertTrue(center_inside, "center should be inside disk");
    assertTrue(!corner_inside, "corner should be outside disk");
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

    const PatchAtom patch_atom = makePatchAtom(atom, rotated, membership, 0.0);
    assertTrue(near(patch_atom.position.x, rotated.x), "builder should preserve rotated x");
    assertTrue(near(patch_atom.position.y, rotated.y), "builder should preserve rotated y");
    assertTrue(near(patch_atom.position.z, rotated.z), "builder should preserve rotated z");
    assertTrue(patch_atom.element == "C", "builder should normalize element");
    assertTrue(near(patch_atom.vdw_radius, 1.70), "builder should assign vdW radius");
    assertTrue(patch_atom.original_atom == &atom, "builder should preserve original atom reference");
    assertTrue(near(patch_atom.radial_xy, membership.radial_xy), "builder should preserve radial xy");
    assertTrue(patch_atom.membership.selected == membership.selected, "builder should preserve membership facts");
}

void testPatchAtomBuilderInfersElementFromNameWhenMissing() {
    const Atom atom(101, "NZ", ' ', "LYS", 'A', 3, ' ', 1.0, 2.0, 3.0, 1.0, 20.0, "", "", false);
    const geometry_symmetry::Vector3 rotated{1.0, 2.0, 3.0};
    const CylinderMembership membership = classifyPatchCylinder(rotated, 5.0);

    const PatchAtom patch_atom = makePatchAtom(atom, rotated, membership, 0.2);
    assertTrue(patch_atom.element == "N", "builder should infer N from atom name when element is missing");
    assertTrue(near(patch_atom.vdw_radius, 1.75), "inferred N should include configured vdW delta");
}

void testStage3FailsWithoutSuccessfulStage2() {
    GeometryPatchSelectionResult stage2_result;
    stage2_result.success = false;
    FoldPatchAnalysisConfig config;
    bool threw = false;
    try {
        (void)runGeometryAnalysisStage3PatchNormalization(stage2_result, config, nullptr);
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
    const auto stage3 = runGeometryAnalysisStage3PatchNormalization(stage2, config, nullptr);

    assertTrue(stage3.success, "Stage 3 should succeed from valid Stage 2 result");
    assertTrue(stage3.analytical_patch.atom_count == stage2.selected_atom_refs.size(),
               "Stage 3 atom count should match selected refs count");
    assertTrue(stage3.analytical_patch.atoms.size() == stage3.analytical_patch.original_atom_refs.size(),
               "Stage 3 atom/reference counts should match");
    assertTrue(stage3.analytical_patch.explicit_vdw_radius_count > 0,
               "Stage 3 should classify explicit vdW assignments");
    assertTrue(stage3.analytical_patch.fallback_vdw_radius_count == 0,
               "known explicit elements should not require fallback assignments");
    for (const PatchAtom& atom : stage3.analytical_patch.atoms) {
        assertTrue(atom.original_atom != nullptr, "each Stage 3 PatchAtom should keep original reference");
        assertTrue(atom.vdw_radius > 0.0, "each Stage 3 PatchAtom should have positive vdW radius");
    }

    assertTrue(std::filesystem::exists(stage2.export_path), "Stage 2 canonical patch export should exist");
    std::filesystem::remove(stage2.export_path);
}

void testStage3InfersAndFallsBackWhenElementMissing() {
    Capsid capsid = makeStage2CapsidWithUnknownElement();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "stage3_unknown";

    GeometryPreparationResult stage1;
    stage1.success = true;

    const auto stage2 = runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);
    const auto stage3 = runGeometryAnalysisStage3PatchNormalization(stage2, config, nullptr);

    assertTrue(stage3.success, "Stage 3 should succeed with missing/unknown element symbols");
    assertTrue(stage3.analytical_patch.explicit_vdw_radius_count == 0,
               "blank element fields should not be counted as explicit assignments");
    assertTrue(stage3.analytical_patch.inferred_vdw_radius_count == 1,
               "missing element should be inferred from supported atom name");
    assertTrue(stage3.analytical_patch.fallback_vdw_radius_count == 1,
               "unsupported unknown element should fall back");

    std::filesystem::remove(stage2.export_path);
}

void testStage4RoleClassificationAndCsvAndPdb() {
    Capsid capsid = makeStage2Capsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 2.0;
    config.min_atoms_in_patch = 2;
    config.debug = true;
    config.output_prefix = "stage4_contacts";

    GeometryPreparationResult stage1;
    stage1.success = true;

    const auto stage2 = runGeometryAnalysisStage2PatchSelection(capsid, config, makeParserConfig(), stage1, nullptr);
    const auto stage3 = runGeometryAnalysisStage3PatchNormalization(stage2, config, nullptr);
    const auto stage4 = runGeometryAnalysisStage4RawSheetDetection(capsid, config, makeParserConfig(), stage3, nullptr);

    assertTrue(stage4.success, "Stage 4 should succeed");
    assertTrue(stage4.valid_node_count > 0, "Stage 4 should have valid nodes");
    assertTrue(stage4.unique_outer_contact_atom_count > 0, "Stage 4 should classify outer contacts");
    assertTrue(stage4.unique_inner_contact_atom_count > 0, "Stage 4 should classify inner contacts");

    bool seen_both = false;
    for (PatchAtomContactRole role : stage4.atom_roles) {
        if (role == PatchAtomContactRole::both) {
            seen_both = true;
            break;
        }
    }
    assertTrue(seen_both, "at least one contact atom should be classified as both");

    assertTrue(std::filesystem::exists(stage4.outer_csv_path), "outer csv should exist");
    assertTrue(std::filesystem::exists(stage4.inner_csv_path), "inner csv should exist");
    assertTrue(std::filesystem::exists(stage4.valid_mask_csv_path), "valid mask csv should exist");
    assertTrue(std::filesystem::exists(stage4.outer_only_mask_csv_path), "outer-only mask csv should exist");
    assertTrue(std::filesystem::exists(stage4.inner_only_mask_csv_path), "inner-only mask csv should exist");
    assertTrue(std::filesystem::exists(stage4.negative_thickness_mask_csv_path),
               "negative-thickness mask csv should exist");
    assertTrue(std::filesystem::exists(stage4.summary_csv_path), "summary csv should exist");
    assertTrue(std::filesystem::exists(stage4.contact_atoms_pdb_path), "contact atom pdb should exist");
    assertTrue(stage4.negative_thickness_node_count == 0, "negative-thickness valid nodes should be zero");
    assertTrue(stage4.both_hit_node_count == stage4.valid_node_count,
               "current policy should classify all valid nodes as both-hit");

    std::ifstream outer(stage4.outer_csv_path);
    std::string header;
    std::getline(outer, header);
    assertTrue(header == "i,j,x,y,inside_disk,valid,z_outer_raw", "outer csv header should match");
    std::string line;
    bool saw_nan = false;
    while (std::getline(outer, line)) {
        if (line.find(",0,nan") != std::string::npos) {
            saw_nan = true;
            break;
        }
    }
    assertTrue(saw_nan, "outer csv should contain invalid nodes with nan");

    std::ifstream pdb(stage4.contact_atoms_pdb_path);
    int atom_records = 0;
    std::string pdb_line;
    while (std::getline(pdb, pdb_line)) {
        if (pdb_line.rfind("ATOM", 0) == 0 || pdb_line.rfind("HETATM", 0) == 0) {
            ++atom_records;
        }
    }
    std::size_t expected_contact_atoms = 0;
    for (PatchAtomContactRole role : stage4.atom_roles) {
        if (role != PatchAtomContactRole::none) {
            ++expected_contact_atoms;
        }
    }
    assertTrue(atom_records == static_cast<int>(expected_contact_atoms),
               "contact PDB should export only original atoms classified as outer/inner contacts");

    std::ifstream summary(stage4.summary_csv_path);
    std::getline(summary, header);
    assertTrue(header.find("stage4_start_utc") != std::string::npos, "summary csv should include timestamps");
    assertTrue(header.find("negative_thickness_nodes") != std::string::npos,
               "summary csv should include thickness counters");

    std::ifstream valid_mask(stage4.valid_mask_csv_path);
    std::getline(valid_mask, header);
    assertTrue(header.find("candidate_patch_atom_count") != std::string::npos,
               "valid mask csv should include candidate atom count");
    assertTrue(header.find("inner_contact_serial_number") != std::string::npos,
               "valid mask csv should include inner serial column");
    assertTrue(header.find("outer_contact_serial_number") != std::string::npos,
               "valid mask csv should include outer serial column");

    std::filesystem::remove(stage2.export_path);
    std::filesystem::remove(stage4.outer_csv_path);
    std::filesystem::remove(stage4.inner_csv_path);
    std::filesystem::remove(stage4.valid_mask_csv_path);
    std::filesystem::remove(stage4.outer_only_mask_csv_path);
    std::filesystem::remove(stage4.inner_only_mask_csv_path);
    std::filesystem::remove(stage4.negative_thickness_mask_csv_path);
    std::filesystem::remove(stage4.summary_csv_path);
    std::filesystem::remove(stage4.contact_atoms_pdb_path);
}

void testStage4DeterministicOutputs() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.min_atoms_in_patch = 2;
    config.debug = true;
    config.output_prefix = "stage4_deterministic";

    const auto first = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(first.success, "first Stage 1-4 run should succeed");
    std::ifstream f1(first.stage4_raw.outer_csv_path);
    const std::string csv1((std::istreambuf_iterator<char>(f1)), std::istreambuf_iterator<char>());

    const auto second = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(second.success, "second Stage 1-4 run should succeed");
    std::ifstream f2(second.stage4_raw.outer_csv_path);
    const std::string csv2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());

    assertTrue(csv1 == csv2, "outer raw csv should be byte-stable across identical runs");
    assertTrue(first.stage4_raw.valid_node_count == second.stage4_raw.valid_node_count,
               "valid-node count should be deterministic");

    std::filesystem::remove(first.stage2_patch.export_path);
    std::filesystem::remove(first.stage4_raw.outer_csv_path);
    std::filesystem::remove(first.stage4_raw.inner_csv_path);
    std::filesystem::remove(first.stage4_raw.valid_mask_csv_path);
    std::filesystem::remove(first.stage4_raw.outer_only_mask_csv_path);
    std::filesystem::remove(first.stage4_raw.inner_only_mask_csv_path);
    std::filesystem::remove(first.stage4_raw.negative_thickness_mask_csv_path);
    std::filesystem::remove(first.stage4_raw.summary_csv_path);
    std::filesystem::remove(first.stage4_raw.contact_atoms_pdb_path);
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

void testStage1ToStage4Integration() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.min_atoms_in_patch = 2;
    config.debug = true;
    config.output_prefix = "stage4_integration";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1-4 integration should succeed");
    assertTrue(result.preparation.success, "Stage 1 should succeed");
    assertTrue(result.stage2_patch.success, "Stage 2 should succeed");
    assertTrue(result.stage3_patch.success, "Stage 3 should succeed");
    assertTrue(result.stage4_raw.success, "Stage 4 should succeed");
    assertTrue(result.stage4_raw.valid_node_count > 0, "Stage 4 valid mask should include valid nodes");
    assertTrue(std::filesystem::exists(result.stage4_raw.contact_atoms_pdb_path), "Stage 4 should export contact atoms");

    std::filesystem::remove(result.stage2_patch.export_path);
    std::filesystem::remove(result.stage4_raw.outer_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_csv_path);
    std::filesystem::remove(result.stage4_raw.valid_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.outer_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.negative_thickness_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.summary_csv_path);
    std::filesystem::remove(result.stage4_raw.contact_atoms_pdb_path);
}

void testStage5RequiresSuccessfulStage4() {
    GeometryStage4RawSheetResult stage4;
    stage4.success = false;
    FoldPatchAnalysisConfig config;
    bool threw = false;
    try {
        (void)runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("before successful Stage 4") != std::string::npos;
    }
    assertTrue(threw, "Stage 5 should require successful Stage 4");
}

void testStage5SeedPromotionAndUniqueSerials() {
    const GeometryStage4RawSheetResult stage4 = makeSyntheticStage4ResultForStage5();
    FoldPatchAnalysisConfig config;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    const auto stage5 = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);

    assertTrue(stage5.success, "Stage 5 should succeed");
    const std::size_t center = nodeIndex(2, 2, stage4.grid.nx);
    assertTrue(stage5.paired_seed_mask[center] == 1, "valid Stage 4 node should be a paired seed");
    assertTrue(near(stage5.z_outer_seed[center], stage4.z_outer_raw[center]),
               "outer seed value should match Stage 4 raw value");
    assertTrue(near(stage5.z_inner_seed[center], stage4.z_inner_raw[center]),
               "inner seed value should match Stage 4 raw value");

    const std::size_t invalid_idx = nodeIndex(0, 2, stage4.grid.nx);
    assertTrue(stage5.paired_seed_mask[invalid_idx] == 0, "invalid Stage 4 node should not be a paired seed");
    assertTrue(std::isnan(stage5.z_outer_seed[invalid_idx]), "invalid node outer seed should remain NaN");
    assertTrue(std::isnan(stage5.z_inner_seed[invalid_idx]), "invalid node inner seed should remain NaN");

    assertTrue(stage5.unique_outer_seed_atom_serials.size() == 5,
               "outer seed serials should be unique and deduplicated");
    assertTrue(stage5.unique_inner_seed_atom_serials.size() == 5,
               "inner seed serials should be unique and deduplicated");
    assertTrue(stage5.unique_outer_seed_atom_serials.front() == 101 &&
                   stage5.unique_outer_seed_atom_serials.back() == 105,
               "outer serial list should be sorted ascending");
    assertTrue(stage5.unique_inner_seed_atom_serials.front() == 201 &&
                   stage5.unique_inner_seed_atom_serials.back() == 205,
               "inner serial list should be sorted ascending");
}

void testStage5SeedSerialTraceabilityUsesAllPairedSeeds() {
    GeometryStage4RawSheetResult stage4 = makeSyntheticStage4ResultForStage5();
    auto markValid = [&](std::size_t i, std::size_t j, double z_outer, double z_inner, int outer_serial,
                         int inner_serial) {
        const std::size_t idx = nodeIndex(i, j, stage4.grid.nx);
        if (stage4.valid_mask[idx] == 0) {
            ++stage4.valid_node_count;
        }
        stage4.valid_mask[idx] = 1;
        stage4.z_outer_raw[idx] = z_outer;
        stage4.z_inner_raw[idx] = z_inner;
        stage4.outer_contact_serial_numbers[idx] = outer_serial;
        stage4.inner_contact_serial_numbers[idx] = inner_serial;
    };

    // Add additional valid nodes near/outside the reliable-core region and at boundary-prone positions.
    markValid(2, 0, 8.7, 4.7, 120, 220);
    markValid(0, 2, 8.8, 4.8, 130, 230);
    markValid(4, 2, 8.9, 4.9, 120, 240); // duplicate outer serial, distinct inner serial.
    stage4.invalid_node_count = stage4.inside_disk_count - stage4.valid_node_count;

    FoldPatchAnalysisConfig config;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.stage5_boundary_margin = 0.75;
    config.stage5_reliable_radius = 0.5;
    config.stage5_support_radius = 1.5;
    config.stage5_min_support_nodes = 4;
    const auto stage5 = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);
    assertTrue(stage5.success, "Stage 5 should succeed");

    std::vector<int> expected_outer;
    std::vector<int> expected_inner;
    for (std::size_t idx = 0; idx < stage5.node_count; ++idx) {
        if (stage5.paired_seed_mask[idx] == 0) {
            continue;
        }
        if (stage4.outer_contact_serial_numbers[idx] > 0) {
            expected_outer.push_back(stage4.outer_contact_serial_numbers[idx]);
        }
        if (stage4.inner_contact_serial_numbers[idx] > 0) {
            expected_inner.push_back(stage4.inner_contact_serial_numbers[idx]);
        }
    }
    std::sort(expected_outer.begin(), expected_outer.end());
    expected_outer.erase(std::unique(expected_outer.begin(), expected_outer.end()), expected_outer.end());
    std::sort(expected_inner.begin(), expected_inner.end());
    expected_inner.erase(std::unique(expected_inner.begin(), expected_inner.end()), expected_inner.end());

    assertTrue(stage5.unique_outer_seed_atom_serials == expected_outer,
               "Stage 5 outer seed serials should match all paired-seed nodes");
    assertTrue(stage5.unique_inner_seed_atom_serials == expected_inner,
               "Stage 5 inner seed serials should match all paired-seed nodes");
}

void testStage5BoundaryExclusionBehavior() {
    const GeometryStage4RawSheetResult stage4 = makeSyntheticStage4ResultForStage5();
    FoldPatchAnalysisConfig config;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.stage5_boundary_margin = 0.5;
    const auto stage5 = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);

    const std::size_t edge = nodeIndex(2, 0, stage4.grid.nx); // x=0, y=-2, r=2
    const std::size_t center = nodeIndex(2, 2, stage4.grid.nx); // r=0
    assertTrue(stage5.boundary_exclusion_mask[edge] == 1, "node with r > R-boundary_margin should be excluded");
    assertTrue(stage5.boundary_exclusion_mask[center] == 0, "center node should not be boundary-excluded");
}

void testStage5InterpolationAdmissibilityBehavior() {
    GeometryStage4RawSheetResult stage4 = makeSyntheticStage4ResultForStage5();
    FoldPatchAnalysisConfig config;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.stage5_boundary_margin = 0.1;
    config.stage5_support_radius = 2.1;
    config.stage5_min_support_nodes = 4;
    const auto stage5 = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);

    const std::size_t fillable = nodeIndex(3, 1, stage4.grid.nx);
    assertTrue(stage5.paired_seed_mask[fillable] == 0, "target node should be non-seed");
    assertTrue(stage5.paired_interp_allowed_mask[fillable] == 1,
               "supported interior non-seed node should be interpolation-allowed");
    assertTrue(stage5.hard_invalid_mask[fillable] == 0, "supported node should not be hard-invalid");

    const std::size_t unsupported = nodeIndex(0, 2, stage4.grid.nx);
    assertTrue(stage5.paired_interp_allowed_mask[unsupported] == 0,
               "unsupported non-seed node should not be interpolation-allowed");
    assertTrue(stage5.hard_invalid_mask[unsupported] == 1, "unsupported node should be hard-invalid");
}

void testStage5ReliableCoreBehavior() {
    const GeometryStage4RawSheetResult stage4 = makeSyntheticStage4ResultForStage5();
    FoldPatchAnalysisConfig config;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.stage5_boundary_margin = 0.1;
    config.stage5_support_radius = 1.5;
    config.stage5_min_support_nodes = 4;
    config.stage5_reliable_radius = 1.0;
    const auto stage5 = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);

    const std::size_t center = nodeIndex(2, 2, stage4.grid.nx);
    const std::size_t outer_ring = nodeIndex(2, 0, stage4.grid.nx);
    const std::size_t hard_invalid = nodeIndex(0, 2, stage4.grid.nx);
    assertTrue(stage5.reliable_core_mask[center] == 1, "center valid seed should be in reliable core");
    assertTrue(stage5.reliable_core_mask[outer_ring] == 0, "node outside reliable radius should be excluded");
    assertTrue(stage5.hard_invalid_mask[hard_invalid] == 1, "sanity: chosen node should be hard-invalid");
    assertTrue(stage5.reliable_core_mask[hard_invalid] == 0, "hard-invalid node should not be reliable core");
}

void testStage5DebugCsvExportAndDeterminism() {
    const GeometryStage4RawSheetResult stage4 = makeSyntheticStage4ResultForStage5();
    FoldPatchAnalysisConfig config;
    config.debug = true;
    config.output_prefix = "stage5_debug";
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.stage5_support_radius = 1.5;
    config.stage5_min_support_nodes = 4;

    const auto first = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);
    assertTrue(std::filesystem::exists(first.outer_seed_csv_path), "outer seed csv should exist");
    assertTrue(std::filesystem::exists(first.inner_seed_csv_path), "inner seed csv should exist");
    assertTrue(std::filesystem::exists(first.paired_seed_mask_csv_path), "paired seed mask csv should exist");
    assertTrue(std::filesystem::exists(first.boundary_exclusion_mask_csv_path),
               "boundary exclusion mask csv should exist");
    assertTrue(std::filesystem::exists(first.interp_allowed_mask_csv_path), "interp allowed mask csv should exist");
    assertTrue(std::filesystem::exists(first.hard_invalid_mask_csv_path), "hard invalid mask csv should exist");
    assertTrue(std::filesystem::exists(first.reliable_core_mask_csv_path), "reliable core mask csv should exist");
    assertTrue(std::filesystem::exists(first.summary_csv_path), "summary csv should exist");

    std::ifstream outer_seed(first.outer_seed_csv_path);
    std::string header;
    std::getline(outer_seed, header);
    assertTrue(header == "i,j,x,y,inside_disk,raw_valid,paired_seed,z_outer_seed",
               "outer seed csv header should match exactly");

    std::ifstream f1(first.interp_allowed_mask_csv_path);
    const std::string csv1((std::istreambuf_iterator<char>(f1)), std::istreambuf_iterator<char>());
    const auto second = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);
    std::ifstream f2(second.interp_allowed_mask_csv_path);
    const std::string csv2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());
    assertTrue(csv1 == csv2, "Stage 5 interpolation mask csv should be byte-stable across identical runs");

    std::filesystem::remove(first.outer_seed_csv_path);
    std::filesystem::remove(first.inner_seed_csv_path);
    std::filesystem::remove(first.paired_seed_mask_csv_path);
    std::filesystem::remove(first.boundary_exclusion_mask_csv_path);
    std::filesystem::remove(first.interp_allowed_mask_csv_path);
    std::filesystem::remove(first.hard_invalid_mask_csv_path);
    std::filesystem::remove(first.reliable_core_mask_csv_path);
    std::filesystem::remove(first.summary_csv_path);
}

void testStage1ToStage5Integration() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.min_atoms_in_patch = 2;
    config.debug = true;
    config.output_prefix = "stage5_integration";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1-5 integration should succeed");
    assertTrue(result.stage5_prep.success, "Stage 5 should succeed in full pipeline");
    assertTrue(result.stage5_prep.paired_seed_node_count > 0, "Stage 5 should contain paired seed nodes");
    assertTrue(result.stage5_prep.reliable_core_node_count > 0, "Stage 5 should contain reliable-core nodes");
    assertTrue(std::filesystem::exists(result.stage5_prep.summary_csv_path), "Stage 5 summary csv should exist");

    std::filesystem::remove(result.stage2_patch.export_path);
    std::filesystem::remove(result.stage4_raw.outer_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_csv_path);
    std::filesystem::remove(result.stage4_raw.valid_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.outer_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.negative_thickness_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.summary_csv_path);
    std::filesystem::remove(result.stage4_raw.contact_atoms_pdb_path);
    std::filesystem::remove(result.stage5_prep.outer_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.inner_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.paired_seed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.boundary_exclusion_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.interp_allowed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.hard_invalid_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.reliable_core_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.summary_csv_path);
}

} // namespace

int main() {
    try {
        testCylinderClassifier();
        testVdwRadiusLookupAndNormalization();
        testPatchAtomBuilderNormalization();
        testPatchAtomBuilderInfersElementFromNameWhenMissing();
        testLineSphereIntersectionHelper();
        testSingleNodeFirstContact();
        testGridConstruction();
        testStage1IdentityFor2_0();
        testStage1NonIdentityFor3_0();
        testStage2PatchSelectionAndTraceability();
        testStage2MinAtomsFailure();
        testStage2RequiresStage1();
        testStage3FailsWithoutSuccessfulStage2();
        testStage2ToStage3Integration();
        testStage3InfersAndFallsBackWhenElementMissing();
        testStage4RoleClassificationAndCsvAndPdb();
        testStage4DeterministicOutputs();
        testStage1ToStage4Integration();
        testStage5RequiresSuccessfulStage4();
        testStage5SeedPromotionAndUniqueSerials();
        testStage5SeedSerialTraceabilityUsesAllPairedSeeds();
        testStage5BoundaryExclusionBehavior();
        testStage5InterpolationAdmissibilityBehavior();
        testStage5ReliableCoreBehavior();
        testStage5DebugCsvExportAndDeterminism();
        testStage1ToStage5Integration();
        std::cout << "All geometry analysis Stage 1/2/3/4/5 tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Geometry analysis test failure: " << e.what() << '\n';
        return 1;
    }
}
