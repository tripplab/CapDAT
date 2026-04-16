#include "geometry_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
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

void removeIfExists(const std::string& path) {
    if (!path.empty()) {
        std::filesystem::remove(path);
    }
}

std::vector<std::string> splitComma(const std::string& line) {
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) {
        parts.push_back(item);
    }
    return parts;
}

std::vector<std::vector<std::string>> readCsvRows(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(in, line)) {
        rows.push_back(splitComma(line));
    }
    return rows;
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

Capsid makeDenseCurvatureCapsid() {
    Capsid capsid("dense_curvature_constructed");
    Chain chain(1, 'A');
    int serial = 1;
    int residue_seq = 1;
    for (int jy = -2; jy <= 2; ++jy) {
        for (int ix = -2; ix <= 2; ++ix) {
            Residue residue("ALA", residue_seq++, ' ', 'A', 1);
            const double x = static_cast<double>(ix) * 0.8;
            const double y = static_cast<double>(jy) * 0.8;
            const double z = 2.0 + (0.1 * x * x) + (0.05 * y * y);
            residue.addAtom(Atom(serial++, "CA", ' ', "ALA", 'A', residue_seq, ' ', x, y, z, 1.0, 20.0, "C", "", false));
            chain.addResidue(residue);
        }
    }
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
    stage4.outer_contact_patch_atom_indices.assign(stage4.node_count, -1);
    stage4.inner_contact_patch_atom_indices.assign(stage4.node_count, -1);
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

    auto markValid = [&](std::size_t i,
                         std::size_t j,
                         double z_outer,
                         double z_inner,
                         int outer_serial,
                         int inner_serial,
                         int outer_patch_idx,
                         int inner_patch_idx) {
        const std::size_t idx = nodeIndex(i, j, stage4.grid.nx);
        stage4.valid_mask[idx] = 1;
        stage4.z_outer_raw[idx] = z_outer;
        stage4.z_inner_raw[idx] = z_inner;
        stage4.outer_contact_serial_numbers[idx] = outer_serial;
        stage4.inner_contact_serial_numbers[idx] = inner_serial;
        stage4.outer_contact_patch_atom_indices[idx] = outer_patch_idx;
        stage4.inner_contact_patch_atom_indices[idx] = inner_patch_idx;
        ++stage4.valid_node_count;
    };

    markValid(2, 2, 8.0, 4.0, 101, 201, 10, 20);
    markValid(1, 2, 8.1, 4.1, 102, 202, 11, 21);
    markValid(3, 2, 8.2, 4.2, 103, 203, 12, 22);
    markValid(2, 1, 8.3, 4.3, 104, 204, 13, 23);
    markValid(2, 3, 8.4, 4.4, 105, 205, 14, 24);
    markValid(1, 1, 8.5, 4.5, 101, 201, 10, 20);
    markValid(3, 3, 8.6, 4.6, 103, 203, 12, 22);

    stage4.invalid_node_count = stage4.inside_disk_count - stage4.valid_node_count;
    return stage4;
}

GeometryStage4RawSheetResult makeStage4ResultWithIndexSerialMismatch() {
    GeometryStage4RawSheetResult stage4 = makeSyntheticStage4ResultForStage5();

    auto markValid = [&](std::size_t i,
                         std::size_t j,
                         double z_outer,
                         double z_inner,
                         int outer_serial,
                         int inner_serial,
                         int outer_patch_idx,
                         int inner_patch_idx) {
        const std::size_t idx = nodeIndex(i, j, stage4.grid.nx);
        if (stage4.valid_mask[idx] == 0) {
            ++stage4.valid_node_count;
        }
        stage4.valid_mask[idx] = 1;
        stage4.z_outer_raw[idx] = z_outer;
        stage4.z_inner_raw[idx] = z_inner;
        stage4.outer_contact_serial_numbers[idx] = outer_serial;
        stage4.inner_contact_serial_numbers[idx] = inner_serial;
        stage4.outer_contact_patch_atom_indices[idx] = outer_patch_idx;
        stage4.inner_contact_patch_atom_indices[idx] = inner_patch_idx;
    };

    // Introduce serial collisions across distinct patch indices.
    markValid(0, 2, 9.1, 4.1, 7, -3, 30, 40);
    markValid(4, 2, 9.2, 4.2, 7, -3, 31, 41);
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

    const auto center_hit = detectRawFirstContactAtNode(0.0, 0.0, atoms);
    assertTrue(center_hit.valid, "one intersecting sphere should be valid");
    assertTrue(center_hit.candidate_patch_atom_count == 1, "single intersecting atom should be counted as candidate");
    assertTrue(near(center_hit.z_inner_raw, 3.0), "single-sphere center ray should use z_low=3");
    assertTrue(near(center_hit.z_outer_raw, 7.0), "single-sphere center ray should use z_high=7");
    assertTrue(center_hit.outer_patch_atom_index == 0 && center_hit.inner_patch_atom_index == 0,
               "single-sphere winner provenance should point to the only atom");

    const auto tangent_hit = detectRawFirstContactAtNode(2.0, 0.0, atoms);
    assertTrue(tangent_hit.valid, "tangent intersection should be valid");
    assertTrue(near(tangent_hit.z_inner_raw, 5.0), "tangent z_low should equal center z");
    assertTrue(near(tangent_hit.z_outer_raw, 5.0), "tangent z_high should equal center z");

    const auto proper_off_center = detectRawFirstContactAtNode(1.0, 0.0, atoms);
    assertTrue(proper_off_center.valid, "proper off-center hit should be valid");
    assertTrue(near(proper_off_center.z_inner_raw, 5.0 - std::sqrt(3.0)), "off-center inner should use z_low");
    assertTrue(near(proper_off_center.z_outer_raw, 5.0 + std::sqrt(3.0)), "off-center outer should use z_high");

    const auto none = detectRawFirstContactAtNode(20.0, 20.0, atoms);
    assertTrue(!none.valid, "no intersection should be invalid");
}

void testStage4FirstContactEnvelopeCompetingWinners() {
    std::vector<PatchAtom> atoms;

    PatchAtom high_center_small_radius;
    high_center_small_radius.position = {0.0, 0.0, 7.0};
    high_center_small_radius.vdw_radius = 0.6;
    atoms.push_back(high_center_small_radius); // z_high = 7.6, z_low = 6.4 at x=y=0

    PatchAtom lower_center_large_radius;
    lower_center_large_radius.position = {0.0, 0.0, 6.0};
    lower_center_large_radius.vdw_radius = 2.0;
    atoms.push_back(lower_center_large_radius); // z_high = 8.0, z_low = 4.0 at x=y=0

    const auto node = detectRawFirstContactAtNode(0.0, 0.0, atoms);
    assertTrue(node.valid, "two intersecting spheres should produce a valid envelope node");
    assertTrue(node.candidate_patch_atom_count == 2, "both intersecting atoms should be counted as candidates");
    assertTrue(near(node.z_outer_raw, 8.0), "outer winner must be the largest z_high envelope contact");
    assertTrue(near(node.z_inner_raw, 4.0), "inner winner must be the smallest z_low envelope contact");
    assertTrue(node.outer_patch_atom_index == 1,
               "outer winner should be selected by z_high, not by highest atom center z");
    assertTrue(node.inner_patch_atom_index == 1, "inner winner should be selected by z_low");
}

void testStage4InnerWinnerCanDisagreeWithLowestCenterZ() {
    std::vector<PatchAtom> atoms;

    PatchAtom lowest_center_but_shallow;
    lowest_center_but_shallow.position = {0.0, 0.0, 4.0};
    lowest_center_but_shallow.vdw_radius = 1.0;
    atoms.push_back(lowest_center_but_shallow); // z_low = 3.0

    PatchAtom higher_center_but_deeper;
    higher_center_but_deeper.position = {0.0, 0.0, 5.0};
    higher_center_but_deeper.vdw_radius = 2.5;
    atoms.push_back(higher_center_but_deeper); // z_low = 2.5

    const auto node = detectRawFirstContactAtNode(0.0, 0.0, atoms);
    assertTrue(node.valid, "node should be valid when both atoms intersect");
    assertTrue(node.inner_patch_atom_index == 1,
               "inner winner should be selected by z_low envelope, not by lowest center z");
    assertTrue(near(node.z_inner_raw, 2.5), "inner raw contact should match winning z_low");
}

void testStage4FirstContactTieStability() {
    std::vector<PatchAtom> atoms;
    PatchAtom a;
    a.position = {0.0, 0.0, 5.0};
    a.vdw_radius = 2.0;
    atoms.push_back(a);

    PatchAtom b = a;
    atoms.push_back(b);

    const auto tie = detectRawFirstContactAtNode(0.0, 0.0, atoms);
    assertTrue(tie.valid, "tie case should remain valid");
    assertTrue(tie.outer_patch_atom_index == 0, "outer tie should preserve first-wins behavior");
    assertTrue(tie.inner_patch_atom_index == 0, "inner tie should preserve first-wins behavior");
}

void testStage4PipelineUsesEnvelopeContactsAndProvenance() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 3.0;
    config.grid_spacing = 0.5;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "stage4_envelope_regression";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1-7 pipeline should succeed for Stage 4 envelope regression fixture");
    assertTrue(result.stage4_raw.valid_node_count > 0, "Stage 4 should produce valid nodes");

    const auto& stage3_atoms = result.stage3_patch.analytical_patch.atoms;
    const auto& stage4 = result.stage4_raw;
    bool checked_node = false;
    for (std::size_t j = 0; j < stage4.grid.ny && !checked_node; ++j) {
        for (std::size_t i = 0; i < stage4.grid.nx && !checked_node; ++i) {
            const std::size_t idx = nodeIndex(i, j, stage4.grid.nx);
            if (stage4.valid_mask[idx] == 0) {
                continue;
            }
            const double x = stage4.grid.x_values[i];
            const double y = stage4.grid.y_values[j];
            double best_high = -std::numeric_limits<double>::infinity();
            double best_low = std::numeric_limits<double>::infinity();
            std::size_t best_high_idx = 0;
            std::size_t best_low_idx = 0;
            std::size_t hits = 0;
            for (std::size_t aidx = 0; aidx < stage3_atoms.size(); ++aidx) {
                const auto intersection = intersectVerticalLineWithSphere(x, y, stage3_atoms[aidx]);
                if (!intersection.intersects) {
                    continue;
                }
                ++hits;
                if (intersection.z_high > best_high + 1e-12) {
                    best_high = intersection.z_high;
                    best_high_idx = aidx;
                }
                if (intersection.z_low < best_low - 1e-12) {
                    best_low = intersection.z_low;
                    best_low_idx = aidx;
                }
            }
            assertTrue(hits == stage4.candidate_patch_atom_counts[idx],
                       "candidate_patch_atom_count should match intersected analytical patch atoms");
            assertTrue(near(stage4.z_outer_raw[idx], best_high),
                       "Stage 4 outer raw value should equal winning ray-sphere z_high envelope contact");
            assertTrue(near(stage4.z_inner_raw[idx], best_low),
                       "Stage 4 inner raw value should equal winning ray-sphere z_low envelope contact");
            assertTrue(stage4.outer_contact_patch_atom_indices[idx] == static_cast<int>(best_high_idx),
                       "outer Stage 4 provenance should point to the winning z_high atom index");
            assertTrue(stage4.inner_contact_patch_atom_indices[idx] == static_cast<int>(best_low_idx),
                       "inner Stage 4 provenance should point to the winning z_low atom index");
            assertTrue(stage4.inner_contact_serial_numbers[idx] ==
                           stage3_atoms[best_low_idx].original_atom->serial() &&
                           stage4.outer_contact_serial_numbers[idx] ==
                           stage3_atoms[best_high_idx].original_atom->serial(),
                       "serial provenance should follow winning z_low/z_high envelope atoms");
            checked_node = true;
        }
    }
    assertTrue(checked_node, "regression fixture must include at least one valid Stage 4 node");

    removeIfExists(result.stage2_patch.export_path);
    removeIfExists(result.stage4_raw.outer_csv_path);
    removeIfExists(result.stage4_raw.inner_csv_path);
    removeIfExists(result.stage4_raw.valid_mask_csv_path);
    removeIfExists(result.stage4_raw.outer_only_mask_csv_path);
    removeIfExists(result.stage4_raw.inner_only_mask_csv_path);
    removeIfExists(result.stage4_raw.negative_thickness_mask_csv_path);
    removeIfExists(result.stage4_raw.summary_csv_path);
    removeIfExists(result.stage4_raw.contact_atoms_pdb_path);
    removeIfExists(result.stage5_prep.outer_seed_csv_path);
    removeIfExists(result.stage5_prep.inner_seed_csv_path);
    removeIfExists(result.stage5_prep.paired_seed_mask_csv_path);
    removeIfExists(result.stage5_prep.boundary_exclusion_mask_csv_path);
    removeIfExists(result.stage5_prep.interp_allowed_mask_csv_path);
    removeIfExists(result.stage5_prep.hard_invalid_mask_csv_path);
    removeIfExists(result.stage5_prep.reliable_core_mask_csv_path);
    removeIfExists(result.stage5_prep.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.inner_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.reconstructed_mask_csv_path);
    removeIfExists(result.stage6_surfaces.final_valid_analysis_mask_csv_path);
    removeIfExists(result.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage6_surfaces.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_obj_path);
    removeIfExists(result.stage6_surfaces.inner_obj_path);
    removeIfExists(result.stage7_smooth.outer_smooth_csv_path);
    removeIfExists(result.stage7_smooth.inner_smooth_csv_path);
    removeIfExists(result.stage7_smooth.smooth_valid_mask_csv_path);
    removeIfExists(result.stage7_smooth.metric_domain_mask_csv_path);
    removeIfExists(result.stage7_smooth.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage7_smooth.summary_csv_path);
    removeIfExists(result.stage7_smooth.outer_mesh_path);
    removeIfExists(result.stage7_smooth.inner_mesh_path);
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
    bool found_identity_frame_note = false;
    for (const std::string& message : result.messages) {
        if (message.find("final working frame: derived_reoriented_frame (already aligned to +Z; workflow recorded without coordinate changes)") !=
            std::string::npos) {
            found_identity_frame_note = true;
            break;
        }
    }
    assertTrue(found_identity_frame_note, "identity Stage 1 logs should clarify frame provenance without coordinate edits");

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

    bool seen_any_contact_role = false;
    for (PatchAtomContactRole role : stage4.atom_roles) {
        if (role != PatchAtomContactRole::none) {
            seen_any_contact_role = true;
            break;
        }
    }
    assertTrue(seen_any_contact_role, "at least one patch atom should be classified as an outer/inner contact");

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
    removeIfExists(stage4.stage3_normalized_atoms_csv_path);
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
    removeIfExists(first.stage4_raw.stage3_normalized_atoms_csv_path);
    removeIfExists(first.stage5_prep.outer_seed_csv_path);
    removeIfExists(first.stage5_prep.inner_seed_csv_path);
    removeIfExists(first.stage5_prep.paired_seed_mask_csv_path);
    removeIfExists(first.stage5_prep.boundary_exclusion_mask_csv_path);
    removeIfExists(first.stage5_prep.interp_allowed_mask_csv_path);
    removeIfExists(first.stage5_prep.hard_invalid_mask_csv_path);
    removeIfExists(first.stage5_prep.reliable_core_mask_csv_path);
    removeIfExists(first.stage5_prep.summary_csv_path);
    removeIfExists(first.stage6_surfaces.outer_reconstructed_csv_path);
    removeIfExists(first.stage6_surfaces.inner_reconstructed_csv_path);
    removeIfExists(first.stage6_surfaces.reconstructed_mask_csv_path);
    removeIfExists(first.stage6_surfaces.final_valid_analysis_mask_csv_path);
    removeIfExists(first.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    removeIfExists(first.stage6_surfaces.summary_csv_path);
    removeIfExists(first.stage6_surfaces.outer_obj_path);
    removeIfExists(first.stage6_surfaces.inner_obj_path);
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
    removeIfExists(result.stage4_raw.stage3_normalized_atoms_csv_path);
    removeIfExists(result.stage5_prep.outer_seed_csv_path);
    removeIfExists(result.stage5_prep.inner_seed_csv_path);
    removeIfExists(result.stage5_prep.paired_seed_mask_csv_path);
    removeIfExists(result.stage5_prep.boundary_exclusion_mask_csv_path);
    removeIfExists(result.stage5_prep.interp_allowed_mask_csv_path);
    removeIfExists(result.stage5_prep.hard_invalid_mask_csv_path);
    removeIfExists(result.stage5_prep.reliable_core_mask_csv_path);
    removeIfExists(result.stage5_prep.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.inner_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.reconstructed_mask_csv_path);
    removeIfExists(result.stage6_surfaces.final_valid_analysis_mask_csv_path);
    removeIfExists(result.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage6_surfaces.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_obj_path);
    removeIfExists(result.stage6_surfaces.inner_obj_path);
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
        expected_outer.push_back(stage4.outer_contact_serial_numbers[idx]);
        expected_inner.push_back(stage4.inner_contact_serial_numbers[idx]);
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

void testStage4SummaryCsvIncludesExplicitPatchAndSerialProvenanceColumns() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.debug = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.min_atoms_in_patch = 2;
    config.output_prefix = "test_stage4_provenance_columns";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1-5 pipeline should succeed");
    std::ifstream summary(result.stage4_raw.summary_csv_path);
    assertTrue(summary.good(), "Stage 4 summary csv should be readable");

    std::string header;
    std::getline(summary, header);
    assertTrue(header.find("unique_outer_contact_patch_atoms") != std::string::npos,
               "Stage 4 summary should include index-based outer provenance");
    assertTrue(header.find("unique_outer_contact_serials") != std::string::npos,
               "Stage 4 summary should include serial-based outer provenance");
    assertTrue(header.find("unique_contact_serial_union") != std::string::npos,
               "Stage 4 summary should include serial union provenance");

    std::filesystem::remove(result.stage2_patch.export_path);
    std::filesystem::remove(result.stage4_raw.outer_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_csv_path);
    std::filesystem::remove(result.stage4_raw.valid_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.outer_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.negative_thickness_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.summary_csv_path);
    std::filesystem::remove(result.stage4_raw.contact_atoms_pdb_path);
    removeIfExists(result.stage4_raw.stage3_normalized_atoms_csv_path);
    std::filesystem::remove(result.stage5_prep.outer_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.inner_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.paired_seed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.boundary_exclusion_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.interp_allowed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.hard_invalid_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.reliable_core_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.inner_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.reconstructed_mask_csv_path);
    removeIfExists(result.stage6_surfaces.final_valid_analysis_mask_csv_path);
    removeIfExists(result.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage6_surfaces.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_obj_path);
    removeIfExists(result.stage6_surfaces.inner_obj_path);
}

void testStage5SeedProvenanceSeparatesPatchIndicesFromSerials() {
    const GeometryStage4RawSheetResult stage4 = makeStage4ResultWithIndexSerialMismatch();
    FoldPatchAnalysisConfig config;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;

    const auto stage5 = runGeometryAnalysisStage5SurfacePreparation(stage4, config, nullptr);
    assertTrue(stage5.success, "Stage 5 should succeed");
    assertTrue(stage5.unique_outer_seed_patch_atom_count > stage5.unique_outer_seed_atom_serials.size(),
               "outer index-based provenance should exceed serial-based provenance in mismatch fixture");
    assertTrue(stage5.unique_inner_seed_patch_atom_count > stage5.unique_inner_seed_atom_serials.size(),
               "inner index-based provenance should exceed serial-based provenance in mismatch fixture");
    assertTrue(std::is_sorted(stage5.unique_outer_seed_patch_atom_indices.begin(),
                              stage5.unique_outer_seed_patch_atom_indices.end()),
               "outer seed patch indices should be sorted");
    assertTrue(std::is_sorted(stage5.unique_outer_seed_atom_serials.begin(), stage5.unique_outer_seed_atom_serials.end()),
               "outer seed serials should be sorted");
    assertTrue(std::find(stage5.unique_outer_seed_atom_serials.begin(), stage5.unique_outer_seed_atom_serials.end(), 7) !=
                   stage5.unique_outer_seed_atom_serials.end(),
               "serial-based provenance should preserve non-positive/positive values from Stage 4");
    assertTrue(std::find(stage5.unique_inner_seed_atom_serials.begin(), stage5.unique_inner_seed_atom_serials.end(), -3) !=
                   stage5.unique_inner_seed_atom_serials.end(),
               "serial-based provenance should include negative serials");
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
    removeIfExists(result.stage4_raw.stage3_normalized_atoms_csv_path);
    std::filesystem::remove(result.stage5_prep.outer_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.inner_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.paired_seed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.boundary_exclusion_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.interp_allowed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.hard_invalid_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.reliable_core_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.inner_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.reconstructed_mask_csv_path);
    removeIfExists(result.stage6_surfaces.final_valid_analysis_mask_csv_path);
    removeIfExists(result.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage6_surfaces.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_obj_path);
    removeIfExists(result.stage6_surfaces.inner_obj_path);
}

GeometryStage5SurfacePrepResult makeSyntheticStage5ResultForStage6() {
    GeometryStage5SurfacePrepResult stage5;
    stage5.success = true;
    stage5.grid = buildStage4RegularGrid(2.0, 1.0);
    stage5.node_count = stage5.grid.nx * stage5.grid.ny;
    stage5.inside_disk_mask.assign(stage5.node_count, 0);
    stage5.paired_seed_mask.assign(stage5.node_count, 0);
    stage5.paired_interp_allowed_mask.assign(stage5.node_count, 0);
    stage5.hard_invalid_mask.assign(stage5.node_count, 0);
    stage5.reliable_core_mask.assign(stage5.node_count, 0);
    stage5.z_outer_seed.assign(stage5.node_count, std::numeric_limits<double>::quiet_NaN());
    stage5.z_inner_seed.assign(stage5.node_count, std::numeric_limits<double>::quiet_NaN());

    for (std::size_t j = 0; j < stage5.grid.ny; ++j) {
        for (std::size_t i = 0; i < stage5.grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, stage5.grid.nx);
            const double x = stage5.grid.x_values[i];
            const double y = stage5.grid.y_values[j];
            if ((x * x) + (y * y) <= 4.0 + 1e-12) {
                stage5.inside_disk_mask[idx] = 1;
                ++stage5.inside_disk_count;
            }
        }
    }

    auto setSeed = [&](std::size_t i, std::size_t j, double zo, double zi) {
        const std::size_t idx = nodeIndex(i, j, stage5.grid.nx);
        stage5.paired_seed_mask[idx] = 1;
        stage5.z_outer_seed[idx] = zo;
        stage5.z_inner_seed[idx] = zi;
        ++stage5.paired_seed_node_count;
    };
    setSeed(2, 2, 10.0, 6.0);
    setSeed(1, 2, 9.0, 5.0);
    setSeed(3, 2, 9.0, 5.0);
    setSeed(2, 1, 9.0, 5.0);
    setSeed(2, 3, 9.0, 5.0);

    const std::size_t center_gap = nodeIndex(1, 1, stage5.grid.nx);
    stage5.paired_interp_allowed_mask[center_gap] = 1;
    stage5.paired_interp_allowed_node_count = 1;
    return stage5;
}

struct SyntheticStage7Inputs {
    GeometryStage5SurfacePrepResult stage5;
    GeometryStage6SurfaceReconstructionResult stage6;
};

SyntheticStage7Inputs makeSyntheticStage7Inputs() {
    SyntheticStage7Inputs inputs;
    inputs.stage5 = makeSyntheticStage5ResultForStage6();
    inputs.stage5.reliable_core_mask.assign(inputs.stage5.node_count, 0);
    inputs.stage5.reliable_core_node_count = 0;
    for (std::size_t idx = 0; idx < inputs.stage5.node_count; ++idx) {
        if (inputs.stage5.inside_disk_mask[idx] != 0) {
            inputs.stage5.reliable_core_mask[idx] = 1;
            ++inputs.stage5.reliable_core_node_count;
        }
    }

    FoldPatchAnalysisConfig stage6_config;
    stage6_config.stage6_export_obj_meshes = false;
    inputs.stage6 = runGeometryAnalysisStage6SurfaceReconstruction(inputs.stage5, stage6_config, nullptr);
    return inputs;
}

void testStage6RequiresSuccessfulStage5() {
    GeometryStage5SurfacePrepResult stage5;
    stage5.success = false;
    FoldPatchAnalysisConfig config;
    bool threw = false;
    try {
        (void)runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("before successful Stage 5") != std::string::npos;
    }
    assertTrue(threw, "Stage 6 should require successful Stage 5");
}

void testStage6MinSeparationDefaultIsZero() {
    FoldPatchAnalysisConfig config;
    assertTrue(near(config.stage6_min_separation, 0.0), "Stage 6 min separation should default to 0.0");
}

void testStage6SeedPreservationAndInterpolation() {
    const GeometryStage5SurfacePrepResult stage5 = makeSyntheticStage5ResultForStage6();
    FoldPatchAnalysisConfig config;
    config.stage6_smoothing_weight = 1.0;
    config.stage6_max_iterations = 200;
    config.stage6_convergence_tolerance = 1e-8;
    config.stage6_export_obj_meshes = false;

    const auto stage6 = runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);
    assertTrue(stage6.success, "Stage 6 should succeed on synthetic Stage 5");

    for (std::size_t idx = 0; idx < stage5.node_count; ++idx) {
        if (stage5.paired_seed_mask[idx] == 0) {
            continue;
        }
        assertTrue(stage6.z_outer_reconstructed[idx] == stage5.z_outer_seed[idx], "outer seed nodes must remain fixed");
        assertTrue(stage6.z_inner_reconstructed[idx] == stage5.z_inner_seed[idx], "inner seed nodes must remain fixed");
    }

    const std::size_t center_gap = nodeIndex(1, 1, stage5.grid.nx);
    assertTrue(std::isfinite(stage6.z_outer_reconstructed[center_gap]), "center interpolation node should be solved");
    assertTrue(std::isfinite(stage6.z_inner_reconstructed[center_gap]), "center interpolation node should be solved");
    assertTrue(near(stage6.z_outer_reconstructed[center_gap], 9.0, 1e-6),
               "center outer reconstruction should approach neighbor average");
    assertTrue(near(stage6.z_inner_reconstructed[center_gap], 5.0, 1e-6),
               "center inner reconstruction should approach neighbor average");
}

void testStage6HardInvalidAndMasks() {
    GeometryStage5SurfacePrepResult stage5 = makeSyntheticStage5ResultForStage6();
    const std::size_t center_gap = nodeIndex(1, 1, stage5.grid.nx);
    stage5.hard_invalid_mask[center_gap] = 1;
    FoldPatchAnalysisConfig config;
    config.stage6_export_obj_meshes = false;
    const auto stage6 = runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);
    assertTrue(std::isnan(stage6.z_outer_reconstructed[center_gap]), "hard-invalid outer node should remain NaN");
    assertTrue(std::isnan(stage6.z_inner_reconstructed[center_gap]), "hard-invalid inner node should remain NaN");
    assertTrue(stage6.reconstructed_mask[center_gap] == 0, "hard-invalid node should not be reconstructed");
    assertTrue(stage6.final_valid_analysis_mask[center_gap] == 0, "hard-invalid node should not be final-valid");
}

void testStage6NonCrossingEnforcementAndConvergenceMetadata() {
    GeometryStage5SurfacePrepResult stage5 = makeSyntheticStage5ResultForStage6();
    const std::size_t idx = nodeIndex(2, 2, stage5.grid.nx);
    stage5.z_outer_seed[idx] = 4.0;
    stage5.z_inner_seed[idx] = 8.0;

    FoldPatchAnalysisConfig config;
    config.stage6_min_separation = 1.5;
    config.stage6_enforce_non_crossing = true;
    config.stage6_export_obj_meshes = false;
    const auto stage6 = runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);

    assertTrue(stage6.non_crossing_adjustment_mask[idx] == 1, "violating node should be adjusted");
    assertTrue(stage6.z_outer_reconstructed[idx] >= stage6.z_inner_reconstructed[idx] + config.stage6_min_separation - 1e-12,
               "enforced separation should hold");
    assertTrue(stage6.outer_iterations_used > 0, "outer iterations should be recorded");
    assertTrue(stage6.inner_iterations_used > 0, "inner iterations should be recorded");
    assertTrue(std::isfinite(stage6.outer_final_max_update), "outer final update should be finite");
    assertTrue(std::isfinite(stage6.inner_final_max_update), "inner final update should be finite");
    assertTrue(stage6.outer_iterations_used <= config.stage6_max_iterations, "outer iterations should honor cap");
    assertTrue(stage6.inner_iterations_used <= config.stage6_max_iterations, "inner iterations should honor cap");
}

void testStage6DebugCsvArtifactsAndObjExportAndDeterminism() {
    GeometryStage5SurfacePrepResult stage5 = makeSyntheticStage5ResultForStage6();
    const std::size_t promoted_seed = nodeIndex(1, 1, stage5.grid.nx);
    stage5.paired_interp_allowed_mask[promoted_seed] = 0;
    stage5.paired_seed_mask[promoted_seed] = 1;
    stage5.z_outer_seed[promoted_seed] = 9.0;
    stage5.z_inner_seed[promoted_seed] = 5.0;
    stage5.paired_interp_allowed_node_count = 0;
    ++stage5.paired_seed_node_count;
    // Enable enough vertices/faces for a valid mesh.
    const std::size_t interp_extra = nodeIndex(3, 1, stage5.grid.nx);
    stage5.paired_interp_allowed_mask[interp_extra] = 1;
    ++stage5.paired_interp_allowed_node_count;

    FoldPatchAnalysisConfig config;
    config.debug = true;
    config.output_prefix = "stage6_debug";
    config.stage6_export_obj_meshes = true;
    config.stage6_mesh_export_format = FoldPatchAnalysisConfig::MeshExportFormat::obj;
    config.stage6_split_in_out_meshes = true;
    const auto first = runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);
    assertTrue(std::filesystem::exists(first.outer_reconstructed_csv_path), "outer reconstructed csv should exist");
    assertTrue(std::filesystem::exists(first.inner_reconstructed_csv_path), "inner reconstructed csv should exist");
    assertTrue(std::filesystem::exists(first.reconstructed_mask_csv_path), "reconstructed mask csv should exist");
    assertTrue(std::filesystem::exists(first.final_valid_analysis_mask_csv_path), "final valid mask csv should exist");
    assertTrue(std::filesystem::exists(first.non_crossing_adjustment_mask_csv_path),
               "non-crossing adjustment mask csv should exist");
    assertTrue(std::filesystem::exists(first.summary_csv_path), "Stage 6 summary csv should exist");
    assertTrue(std::filesystem::exists(first.outer_obj_path), "outer obj should exist");
    assertTrue(std::filesystem::exists(first.inner_obj_path), "inner obj should exist");
    assertTrue(first.outer_obj_vertex_count > 0 && first.inner_obj_vertex_count > 0, "OBJ vertex counts should be nonzero");
    assertTrue(first.outer_obj_face_count >= 0 && first.inner_obj_face_count >= 0, "OBJ face counts should be tracked");

    std::ifstream outer_csv(first.outer_reconstructed_csv_path);
    std::string header;
    std::getline(outer_csv, header);
    assertTrue(header ==
                   "i,j,x,y,inside_disk,paired_seed,paired_interp_allowed,hard_invalid,reconstructed,z_outer_reconstructed",
               "Stage 6 outer reconstructed CSV header should match");

    std::ifstream obj_file(first.outer_obj_path);
    std::string line;
    bool found_v = false;
    bool found_f = false;
    while (std::getline(obj_file, line)) {
        if (line.rfind("v ", 0) == 0) {
            found_v = true;
        }
        if (line.rfind("f ", 0) == 0) {
            found_f = true;
        }
    }
    assertTrue(found_v, "OBJ should contain at least one vertex line");
    assertTrue(found_f, "OBJ should contain at least one face line when valid cells exist");

    const auto second = runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);
    std::ifstream csv1(first.outer_reconstructed_csv_path);
    std::ifstream csv2(second.outer_reconstructed_csv_path);
    const std::string csv_first((std::istreambuf_iterator<char>(csv1)), std::istreambuf_iterator<char>());
    const std::string csv_second((std::istreambuf_iterator<char>(csv2)), std::istreambuf_iterator<char>());
    assertTrue(csv_first == csv_second, "Stage 6 CSV output should be deterministic");
    std::ifstream obj1(first.outer_obj_path);
    std::ifstream obj2(second.outer_obj_path);
    const std::string obj_first((std::istreambuf_iterator<char>(obj1)), std::istreambuf_iterator<char>());
    const std::string obj_second((std::istreambuf_iterator<char>(obj2)), std::istreambuf_iterator<char>());
    assertTrue(obj_first == obj_second, "Stage 6 OBJ output should be deterministic");

    GeometryStage5SurfacePrepResult with_hole = stage5;
    const std::size_t hole_idx = nodeIndex(2, 1, with_hole.grid.nx);
    with_hole.hard_invalid_mask[hole_idx] = 1;
    const auto hole_result = runGeometryAnalysisStage6SurfaceReconstruction(with_hole, config, nullptr);
    assertTrue(hole_result.outer_obj_face_count < first.outer_obj_face_count,
               "introducing an invalid hole should reduce emitted faces by suppressing partial invalid cells");
    std::ifstream hole_obj(hole_result.outer_obj_path);
    std::string hole_line;
    while (std::getline(hole_obj, hole_line)) {
        if (hole_line.rfind("f ", 0) == 0) {
            assertTrue(hole_line.find("-1") == std::string::npos, "OBJ faces should never reference invalid vertices");
        }
    }

    std::filesystem::remove(first.outer_reconstructed_csv_path);
    std::filesystem::remove(first.inner_reconstructed_csv_path);
    std::filesystem::remove(first.reconstructed_mask_csv_path);
    std::filesystem::remove(first.final_valid_analysis_mask_csv_path);
    std::filesystem::remove(first.non_crossing_adjustment_mask_csv_path);
    std::filesystem::remove(first.summary_csv_path);
    std::filesystem::remove(first.outer_obj_path);
    std::filesystem::remove(first.inner_obj_path);
    std::filesystem::remove(hole_result.outer_reconstructed_csv_path);
    std::filesystem::remove(hole_result.inner_reconstructed_csv_path);
    std::filesystem::remove(hole_result.reconstructed_mask_csv_path);
    std::filesystem::remove(hole_result.final_valid_analysis_mask_csv_path);
    std::filesystem::remove(hole_result.non_crossing_adjustment_mask_csv_path);
    std::filesystem::remove(hole_result.summary_csv_path);
    std::filesystem::remove(hole_result.outer_obj_path);
    std::filesystem::remove(hole_result.inner_obj_path);
}

void testStage6StlExport() {
    GeometryStage5SurfacePrepResult stage5 = makeSyntheticStage5ResultForStage6();
    const std::size_t promoted_seed = nodeIndex(1, 1, stage5.grid.nx);
    stage5.paired_interp_allowed_mask[promoted_seed] = 0;
    stage5.paired_seed_mask[promoted_seed] = 1;
    stage5.z_outer_seed[promoted_seed] = 9.0;
    stage5.z_inner_seed[promoted_seed] = 5.0;
    stage5.paired_interp_allowed_node_count = 0;
    ++stage5.paired_seed_node_count;

    FoldPatchAnalysisConfig config;
    config.output_prefix = "stage6_stl";
    config.stage6_export_obj_meshes = true;
    config.stage6_mesh_export_format = FoldPatchAnalysisConfig::MeshExportFormat::stl;
    config.stage6_split_in_out_meshes = true;

    const auto stage6 = runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);
    assertTrue(std::filesystem::exists(stage6.outer_obj_path), "outer stl should exist");
    assertTrue(std::filesystem::exists(stage6.inner_obj_path), "inner stl should exist");
    assertTrue(stage6.outer_obj_path.find(".stl") != std::string::npos, "outer mesh path should use stl extension");

    std::ifstream stl_file(stage6.outer_obj_path);
    std::string first_line;
    std::getline(stl_file, first_line);
    assertTrue(first_line.rfind("solid ", 0) == 0, "STL should begin with a solid header");

    std::filesystem::remove(stage6.outer_obj_path);
    std::filesystem::remove(stage6.inner_obj_path);
}

void testStage6CombinedExportDefault() {
    GeometryStage5SurfacePrepResult stage5 = makeSyntheticStage5ResultForStage6();
    FoldPatchAnalysisConfig config;
    config.output_prefix = "stage6_combined";
    config.stage6_export_obj_meshes = true;
    config.stage6_mesh_export_format = FoldPatchAnalysisConfig::MeshExportFormat::obj;

    const auto stage6 = runGeometryAnalysisStage6SurfaceReconstruction(stage5, config, nullptr);
    assertTrue(stage6.outer_obj_path == stage6.inner_obj_path, "combined export should use one shared mesh path");
    assertTrue(stage6.outer_obj_path.find("_stage6") == std::string::npos,
               "combined export path should not include _stage6 token");
    assertTrue(std::filesystem::exists(stage6.outer_obj_path), "combined mesh should exist");

    std::ifstream mesh_file(stage6.outer_obj_path);
    std::string mesh_text((std::istreambuf_iterator<char>(mesh_file)), std::istreambuf_iterator<char>());
    assertTrue(mesh_text.find("o outer_surface") != std::string::npos, "combined OBJ should contain outer object");
    assertTrue(mesh_text.find("o inner_surface") != std::string::npos, "combined OBJ should contain inner object");

    std::filesystem::remove(stage6.outer_obj_path);
}

void testStage7RequiresSuccessfulStage6() {
    GeometryStage6SurfaceReconstructionResult stage6;
    stage6.success = false;
    GeometryStage5SurfacePrepResult stage5 = makeSyntheticStage5ResultForStage6();
    stage5.reliable_core_mask.assign(stage5.node_count, 1);
    stage5.reliable_core_node_count = stage5.node_count;
    FoldPatchAnalysisConfig config;
    bool threw = false;
    try {
        (void)runGeometryAnalysisStage7SurfaceSmoothing(stage6, stage5, config, nullptr);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("before successful Stage 6") != std::string::npos;
    }
    assertTrue(threw, "Stage 7 should require successful Stage 6");
}

void testStage7SeedPinSwitchAndRoughnessReduction() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.stage7_export_meshes = false;
    config.stage7_smoothing_weight = 2.0;
    config.stage7_max_iterations = 100;

    const std::size_t seed_idx = nodeIndex(2, 2, inputs.stage5.grid.nx);
    const std::size_t noisy_idx = nodeIndex(1, 1, inputs.stage5.grid.nx);
    inputs.stage6.z_outer_reconstructed[noisy_idx] += 8.0;
    inputs.stage6.z_inner_reconstructed[noisy_idx] += 8.0;
    inputs.stage6.reconstructed_mask[noisy_idx] = 1;

    const double before_seed = inputs.stage6.z_outer_reconstructed[seed_idx];
    const double before_rough = std::fabs(inputs.stage6.z_outer_reconstructed[noisy_idx] - before_seed);
    const auto stage7_default =
        runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(stage7_default.z_outer_smooth[seed_idx] != before_seed,
               "Stage 7 should update seed-supported outer values when pin-seed is OFF");
    const double after_rough =
        std::fabs(stage7_default.z_outer_smooth[noisy_idx] - stage7_default.z_outer_smooth[seed_idx]);
    assertTrue(after_rough < before_rough, "Stage 7 smoothing should reduce local roughness on non-seed nodes");

    config.stage7_preserve_seed_values = true;
    const auto stage7_pinned = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(stage7_pinned.z_outer_smooth[seed_idx] == before_seed,
               "Stage 7 should preserve seed-supported outer values when pin-seed is ON");
}

void testStage7ReliableCoreDefinesMetricDomain() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    inputs.stage5.reliable_core_mask.assign(inputs.stage5.node_count, 0);
    inputs.stage5.reliable_core_node_count = 0;
    const std::size_t core_idx = nodeIndex(2, 2, inputs.stage5.grid.nx);
    const std::size_t core_idx2 = nodeIndex(2, 1, inputs.stage5.grid.nx);
    inputs.stage5.reliable_core_mask[core_idx] = 1;
    inputs.stage5.reliable_core_mask[core_idx2] = 1;
    inputs.stage5.reliable_core_node_count = 2;

    FoldPatchAnalysisConfig config;
    config.stage7_export_meshes = false;
    const auto stage7 = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    for (std::size_t idx = 0; idx < stage7.node_count; ++idx) {
        const uint8_t expected =
            (stage7.smooth_valid_mask[idx] != 0 && stage7.reliable_core_mask[idx] != 0) ? static_cast<uint8_t>(1) : 0;
        assertTrue(stage7.metric_domain_mask[idx] == expected, "metric domain must equal smooth_valid AND reliable_core");
    }
}

void testStage7NonCrossingAndMetadata() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    const std::size_t idx = nodeIndex(2, 2, inputs.stage5.grid.nx);
    inputs.stage6.z_outer_reconstructed[idx] = 2.0;
    inputs.stage6.z_inner_reconstructed[idx] = 4.5;

    FoldPatchAnalysisConfig config;
    config.stage7_export_meshes = false;
    config.stage7_preserve_seed_values = true;
    config.stage7_enforce_non_crossing = true;
    config.stage7_min_separation = 1.25;
    const auto stage7 = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);

    assertTrue(stage7.smooth_non_crossing_adjustment_mask[idx] == 1, "Stage 7 should mark adjusted non-crossing nodes");
    assertTrue(stage7.z_outer_smooth[idx] >= stage7.z_inner_smooth[idx] + config.stage7_min_separation - 1e-12,
               "Stage 7 should enforce minimum smooth separation");
    assertTrue(stage7.normal_orientation_outer == "toward_positive_z", "outer normal convention metadata should be explicit");
    assertTrue(stage7.normal_orientation_inner == "toward_negative_z", "inner normal convention metadata should be explicit");
    assertTrue(stage7.local_thickness_definition == "vertical_z_difference",
               "local thickness definition metadata should be vertical Z separation");
    assertTrue(stage7.metric_surface_definition ==
                   "stage6_reconstructed_then_stage7_smoothed_restricted_to_reliable_core",
               "metric-surface definition metadata should match Stage 7 contract");
}

void testStage7MeshExportAndDeterminism() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.debug = true;
    config.output_prefix = "stage7_debug";
    config.stage7_export_meshes = true;
    config.stage6_mesh_export_format = FoldPatchAnalysisConfig::MeshExportFormat::obj;
    config.stage6_split_in_out_meshes = true;

    const auto first = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(std::filesystem::exists(first.outer_mesh_path), "Stage 7 outer mesh export should exist");
    assertTrue(std::filesystem::exists(first.inner_mesh_path), "Stage 7 inner mesh export should exist");
    assertTrue(std::filesystem::exists(first.summary_csv_path), "Stage 7 summary CSV should exist");
    assertTrue(first.outer_mesh_vertex_count > 0, "Stage 7 outer mesh should include vertices");

    const auto second = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    std::ifstream csv1(first.outer_smooth_csv_path);
    std::ifstream csv2(second.outer_smooth_csv_path);
    const std::string csv_first((std::istreambuf_iterator<char>(csv1)), std::istreambuf_iterator<char>());
    const std::string csv_second((std::istreambuf_iterator<char>(csv2)), std::istreambuf_iterator<char>());
    assertTrue(csv_first == csv_second, "Stage 7 CSV output should be deterministic");

    std::ifstream mesh1(first.outer_mesh_path);
    std::ifstream mesh2(second.outer_mesh_path);
    const std::string mesh_first((std::istreambuf_iterator<char>(mesh1)), std::istreambuf_iterator<char>());
    const std::string mesh_second((std::istreambuf_iterator<char>(mesh2)), std::istreambuf_iterator<char>());
    assertTrue(mesh_first == mesh_second, "Stage 7 mesh output should be deterministic");

    removeIfExists(first.outer_smooth_csv_path);
    removeIfExists(first.inner_smooth_csv_path);
    removeIfExists(first.smooth_valid_mask_csv_path);
    removeIfExists(first.metric_domain_mask_csv_path);
    removeIfExists(first.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(first.summary_csv_path);
    removeIfExists(first.outer_mesh_path);
    removeIfExists(first.inner_mesh_path);
}

double computeRoughnessLaplacianEnergy(const Stage4GridDescriptor& grid,
                                       const std::vector<uint8_t>& mask,
                                       const std::vector<double>& values) {
    double energy = 0.0;
    for (std::size_t j = 0; j < grid.ny; ++j) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, grid.nx);
            if (mask[idx] == 0 || !std::isfinite(values[idx])) {
                continue;
            }
            double sum_neighbors = 0.0;
            std::size_t count = 0;
            if (i > 0) {
                const std::size_t nidx = nodeIndex(i - 1, j, grid.nx);
                if (mask[nidx] != 0 && std::isfinite(values[nidx])) {
                    sum_neighbors += values[nidx];
                    ++count;
                }
            }
            if (i + 1 < grid.nx) {
                const std::size_t nidx = nodeIndex(i + 1, j, grid.nx);
                if (mask[nidx] != 0 && std::isfinite(values[nidx])) {
                    sum_neighbors += values[nidx];
                    ++count;
                }
            }
            if (j > 0) {
                const std::size_t nidx = nodeIndex(i, j - 1, grid.nx);
                if (mask[nidx] != 0 && std::isfinite(values[nidx])) {
                    sum_neighbors += values[nidx];
                    ++count;
                }
            }
            if (j + 1 < grid.ny) {
                const std::size_t nidx = nodeIndex(i, j + 1, grid.nx);
                if (mask[nidx] != 0 && std::isfinite(values[nidx])) {
                    sum_neighbors += values[nidx];
                    ++count;
                }
            }
            const double lap = (static_cast<double>(count) * values[idx]) - sum_neighbors;
            energy += lap * lap;
        }
    }
    return energy;
}

void testStage7MethodDispatchAndMetadata() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.stage7_export_meshes = false;
    const auto legacy = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(legacy.stage7_method_label == "smooth", "default Stage 7 method must remain smooth");

    config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;
    const auto fit = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(fit.stage7_method_label == "thin_plate_grid_fit", "explicit Stage 7 fit mode label should be set");
}

void testStage7FitDomainWeightingMetadata() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.stage7_export_meshes = false;
    config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;
    const auto stage7 = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(stage7.stage7_fit_active_node_count > 0, "fit domain should include reconstructed nodes");
    assertTrue(stage7.stage7_fit_seed_like_node_count > 0, "fit metadata should count paired-seed nodes");
    assertTrue(stage7.stage7_fit_interp_like_node_count > 0, "fit metadata should count reconstructed non-seed nodes");
    assertTrue(stage7.stage7_data_weight_policy.find("paired_seed") != std::string::npos,
               "fit metadata should describe fidelity weighting policy");
}

void testStage7BoundaryModes() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.stage7_export_meshes = false;
    config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;
    config.stage7_lambda = 2.0;
    config.stage7_data_weight_seed = 0.8;
    config.stage7_data_weight_interp = 0.2;
    const std::size_t boundary_idx = nodeIndex(1, 2, inputs.stage5.grid.nx);
    const double boundary_stage6 = inputs.stage6.z_outer_reconstructed[boundary_idx];

    config.stage7_boundary_condition_mode = FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::fixed_to_stage6;
    const auto fixed = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(near(fixed.z_outer_smooth[boundary_idx], boundary_stage6, 1e-9),
               "fixed_to_stage6 should preserve boundary reconstructed values");

    config.stage7_boundary_condition_mode = FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::free;
    const auto free = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    config.stage7_boundary_condition_mode = FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::soft_to_stage6;
    const auto soft = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(std::fabs(soft.z_outer_smooth[boundary_idx] - free.z_outer_smooth[boundary_idx]) > 1e-12,
               "soft boundary mode should produce deterministic differences from free mode");
}

void testStage7FitVsLegacyAndRoughness() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    const std::size_t noisy_idx = nodeIndex(1, 1, inputs.stage5.grid.nx);
    inputs.stage6.z_outer_reconstructed[noisy_idx] += 4.0;
    inputs.stage6.z_inner_reconstructed[noisy_idx] += 4.0;
    inputs.stage6.reconstructed_mask[noisy_idx] = 1;

    FoldPatchAnalysisConfig legacy_config;
    legacy_config.stage7_export_meshes = false;
    const auto legacy = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, legacy_config, nullptr);

    FoldPatchAnalysisConfig fit_config = legacy_config;
    fit_config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;
    const auto fit = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, fit_config, nullptr);

    assertTrue(std::fabs(legacy.z_outer_smooth[noisy_idx] - fit.z_outer_smooth[noisy_idx]) > 1e-10,
               "thin-plate fit should differ from legacy smooth output");
    assertTrue(fit.metric_domain_node_count > 0, "fit method should keep metric domain non-empty");

    const double stage6_rough =
        computeRoughnessLaplacianEnergy(inputs.stage6.grid, inputs.stage6.reconstructed_mask, inputs.stage6.z_outer_reconstructed);
    const double fit_rough = computeRoughnessLaplacianEnergy(fit.grid, fit.smooth_valid_mask, fit.z_outer_smooth);
    assertTrue(fit_rough < stage6_rough, "thin-plate fit should reduce roughness relative to Stage 6 target field");
}

void testStage7ThinPlateDeterminism() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.debug = true;
    config.output_prefix = "stage7_fit_debug";
    config.stage7_export_meshes = false;
    config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;

    const auto first = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    const auto second = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    std::ifstream csv1(first.outer_smooth_csv_path);
    std::ifstream csv2(second.outer_smooth_csv_path);
    const std::string csv_first((std::istreambuf_iterator<char>(csv1)), std::istreambuf_iterator<char>());
    const std::string csv_second((std::istreambuf_iterator<char>(csv2)), std::istreambuf_iterator<char>());
    assertTrue(csv_first == csv_second, "thin-plate Stage 7 csv should be byte-deterministic");
    assertTrue(first.outer_fit_final_residual == second.outer_fit_final_residual,
               "thin-plate Stage 7 residual metadata should be deterministic");

    removeIfExists(first.outer_smooth_csv_path);
    removeIfExists(first.inner_smooth_csv_path);
    removeIfExists(first.smooth_valid_mask_csv_path);
    removeIfExists(first.metric_domain_mask_csv_path);
    removeIfExists(first.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(first.summary_csv_path);
}

void testStage7S7S6DeltaCsvsThinPlate() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.debug = true;
    config.output_prefix = "stage7_s7s6_fit";
    config.stage7_export_meshes = false;
    config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;
    config.stage7_export_s7s6_deltas = true;

    const std::size_t valid_idx = nodeIndex(2, 2, inputs.stage5.grid.nx);
    const std::size_t nan_idx = nodeIndex(0, 0, inputs.stage5.grid.nx);
    const std::size_t interp_idx = nodeIndex(1, 1, inputs.stage5.grid.nx);
    inputs.stage5.hard_invalid_mask[valid_idx] = 1;
    inputs.stage6.z_outer_reconstructed[nan_idx] = std::numeric_limits<double>::quiet_NaN();
    inputs.stage6.z_inner_reconstructed[nan_idx] = std::numeric_limits<double>::quiet_NaN();

    const auto first = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(std::filesystem::exists(first.outer_delta_csv_path), "outer delta csv should exist");
    assertTrue(std::filesystem::exists(first.inner_delta_csv_path), "inner delta csv should exist");
    assertTrue(std::filesystem::exists(first.thickness_delta_csv_path), "thickness delta csv should exist");

    std::ifstream outer_header_file(first.outer_delta_csv_path);
    std::ifstream inner_header_file(first.inner_delta_csv_path);
    std::ifstream thick_header_file(first.thickness_delta_csv_path);
    std::string outer_header;
    std::string inner_header;
    std::string thick_header;
    std::getline(outer_header_file, outer_header);
    std::getline(inner_header_file, inner_header);
    std::getline(thick_header_file, thick_header);
    assertTrue(outer_header ==
                   "i,j,x,y,inside_disk,reconstructed,smooth_valid,metric_domain,paired_seed,interp_allowed,hard_invalid,"
                   "reliable_core,z_outer_stage6,z_outer_stage7,delta_outer",
               "outer delta csv header should match schema");
    assertTrue(inner_header ==
                   "i,j,x,y,inside_disk,reconstructed,smooth_valid,metric_domain,paired_seed,interp_allowed,hard_invalid,"
                   "reliable_core,z_inner_stage6,z_inner_stage7,delta_inner",
               "inner delta csv header should match schema");
    assertTrue(thick_header ==
                   "i,j,x,y,inside_disk,reconstructed,smooth_valid,metric_domain,paired_seed,interp_allowed,hard_invalid,"
                   "reliable_core,t_stage6,t_stage7,delta_t",
               "thickness delta csv header should match schema");

    const auto outer_rows = readCsvRows(first.outer_delta_csv_path);
    const auto inner_rows = readCsvRows(first.inner_delta_csv_path);
    const auto thickness_rows = readCsvRows(first.thickness_delta_csv_path);
    const std::size_t expected_row_count = (inputs.stage5.grid.nx * inputs.stage5.grid.ny) + 1;
    assertTrue(outer_rows.size() == expected_row_count, "outer delta csv should include one row per grid node");
    assertTrue(inner_rows.size() == expected_row_count, "inner delta csv should include one row per grid node");
    assertTrue(thickness_rows.size() == expected_row_count, "thickness delta csv should include one row per grid node");

    const std::size_t outer_row = valid_idx + 1;
    const std::size_t inner_row = valid_idx + 1;
    const std::size_t thick_row = valid_idx + 1;
    const double z_outer_stage6 = std::stod(outer_rows[outer_row][12]);
    const double z_outer_stage7 = std::stod(outer_rows[outer_row][13]);
    const double delta_outer = std::stod(outer_rows[outer_row][14]);
    assertTrue(near(z_outer_stage6, inputs.stage6.z_outer_reconstructed[valid_idx], 1e-6),
               "outer stage6 column should match stage6 field");
    assertTrue(near(z_outer_stage7, first.z_outer_smooth[valid_idx], 1e-5), "outer stage7 column should match stage7 field");
    assertTrue(near(delta_outer, first.z_outer_smooth[valid_idx] - inputs.stage6.z_outer_reconstructed[valid_idx], 1e-5),
               "outer delta should equal stage7-stage6");
    const double z_inner_stage6 = std::stod(inner_rows[inner_row][12]);
    const double z_inner_stage7 = std::stod(inner_rows[inner_row][13]);
    const double delta_inner = std::stod(inner_rows[inner_row][14]);
    assertTrue(near(z_inner_stage6, inputs.stage6.z_inner_reconstructed[valid_idx], 1e-6),
               "inner stage6 column should match stage6 field");
    assertTrue(near(z_inner_stage7, first.z_inner_smooth[valid_idx], 1e-5), "inner stage7 column should match stage7 field");
    assertTrue(near(delta_inner, first.z_inner_smooth[valid_idx] - inputs.stage6.z_inner_reconstructed[valid_idx], 1e-5),
               "inner delta should equal stage7-stage6");
    const double t_stage6 = std::stod(thickness_rows[thick_row][12]);
    const double t_stage7 = std::stod(thickness_rows[thick_row][13]);
    const double delta_t = std::stod(thickness_rows[thick_row][14]);
    const double expected_t_stage6 =
        inputs.stage6.z_outer_reconstructed[valid_idx] - inputs.stage6.z_inner_reconstructed[valid_idx];
    const double expected_t_stage7 = first.z_outer_smooth[valid_idx] - first.z_inner_smooth[valid_idx];
    assertTrue(near(t_stage6, expected_t_stage6, 1e-6), "thickness stage6 column should match stage6 field");
    assertTrue(near(t_stage7, expected_t_stage7, 1e-5), "thickness stage7 column should match stage7 field");
    assertTrue(near(delta_t, expected_t_stage7 - expected_t_stage6, 1e-5), "thickness delta should equal stage7-stage6");

    const std::size_t nan_row = nan_idx + 1;
    assertTrue(outer_rows[nan_row][12] == "nan" && outer_rows[nan_row][14] == "nan",
               "outer delta should write nan for invalid nodes");
    assertTrue(inner_rows[nan_row][12] == "nan" && inner_rows[nan_row][14] == "nan",
               "inner delta should write nan for invalid nodes");
    assertTrue(thickness_rows[nan_row][12] == "nan" && thickness_rows[nan_row][14] == "nan",
               "thickness delta should write nan for invalid nodes");

    const std::size_t interp_row = interp_idx + 1;
    assertTrue(outer_rows[outer_row][8] == "1", "paired_seed mask column should match underlying mask");
    assertTrue(outer_rows[interp_row][9] == "1", "interp_allowed mask column should match underlying mask");
    assertTrue(outer_rows[outer_row][7] == "1", "metric_domain mask column should match underlying mask");
    assertTrue(outer_rows[outer_row][10] == "1", "hard_invalid mask column should match underlying mask");

    const auto second = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    std::ifstream first_delta(first.outer_delta_csv_path);
    std::ifstream second_delta(second.outer_delta_csv_path);
    const std::string first_text((std::istreambuf_iterator<char>(first_delta)), std::istreambuf_iterator<char>());
    const std::string second_text((std::istreambuf_iterator<char>(second_delta)), std::istreambuf_iterator<char>());
    assertTrue(first_text == second_text, "Stage 7 s7s6 delta csv should be byte-deterministic");

    removeIfExists(first.outer_smooth_csv_path);
    removeIfExists(first.inner_smooth_csv_path);
    removeIfExists(first.smooth_valid_mask_csv_path);
    removeIfExists(first.metric_domain_mask_csv_path);
    removeIfExists(first.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(first.summary_csv_path);
    removeIfExists(first.outer_delta_csv_path);
    removeIfExists(first.inner_delta_csv_path);
    removeIfExists(first.thickness_delta_csv_path);
}

void testStage7S7S6DeltaCsvsLegacySmooth() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.debug = true;
    config.output_prefix = "stage7_s7s6_smooth";
    config.stage7_export_meshes = false;
    config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::smooth;
    config.stage7_export_s7s6_deltas = true;

    const auto stage7 = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(std::filesystem::exists(stage7.outer_delta_csv_path), "legacy smooth should also emit outer delta csv");
    assertTrue(std::filesystem::exists(stage7.inner_delta_csv_path), "legacy smooth should also emit inner delta csv");
    assertTrue(std::filesystem::exists(stage7.thickness_delta_csv_path),
               "legacy smooth should also emit thickness delta csv");

    removeIfExists(stage7.outer_smooth_csv_path);
    removeIfExists(stage7.inner_smooth_csv_path);
    removeIfExists(stage7.smooth_valid_mask_csv_path);
    removeIfExists(stage7.metric_domain_mask_csv_path);
    removeIfExists(stage7.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(stage7.summary_csv_path);
    removeIfExists(stage7.outer_delta_csv_path);
    removeIfExists(stage7.inner_delta_csv_path);
    removeIfExists(stage7.thickness_delta_csv_path);
}

void testStage7S7S6DeltasDisabledByDefault() {
    SyntheticStage7Inputs inputs = makeSyntheticStage7Inputs();
    FoldPatchAnalysisConfig config;
    config.debug = true;
    config.output_prefix = "stage7_s7s6_disabled";
    config.stage7_export_meshes = false;
    // default should remain disabled unless explicitly requested.
    assertTrue(!config.stage7_export_s7s6_deltas, "s7s6 delta export flag should default to false");

    const auto stage7 = runGeometryAnalysisStage7SurfaceSmoothing(inputs.stage6, inputs.stage5, config, nullptr);
    assertTrue(stage7.outer_delta_csv_path.empty(), "outer delta path should remain empty when disabled");
    assertTrue(stage7.inner_delta_csv_path.empty(), "inner delta path should remain empty when disabled");
    assertTrue(stage7.thickness_delta_csv_path.empty(), "thickness delta path should remain empty when disabled");
    assertTrue(near(stage7.max_abs_outer_delta, 0.0), "max abs outer delta should remain zero when disabled");
    assertTrue(near(stage7.max_abs_inner_delta, 0.0), "max abs inner delta should remain zero when disabled");
    assertTrue(near(stage7.max_abs_thickness_delta, 0.0),
               "max abs thickness delta should remain zero when disabled");
    assertTrue(near(stage7.mean_abs_outer_delta, 0.0), "mean abs outer delta should remain zero when disabled");
    assertTrue(near(stage7.mean_abs_inner_delta, 0.0), "mean abs inner delta should remain zero when disabled");
    assertTrue(near(stage7.mean_abs_thickness_delta, 0.0),
               "mean abs thickness delta should remain zero when disabled");

    removeIfExists(stage7.outer_smooth_csv_path);
    removeIfExists(stage7.inner_smooth_csv_path);
    removeIfExists(stage7.smooth_valid_mask_csv_path);
    removeIfExists(stage7.metric_domain_mask_csv_path);
    removeIfExists(stage7.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(stage7.summary_csv_path);
}

struct SyntheticStage8Inputs {
    GeometryStage5SurfacePrepResult stage5;
    GeometryStage7SmoothedSurfaceResult stage7;
};

SyntheticStage8Inputs makeSyntheticStage8QuadraticInputs() {
    SyntheticStage8Inputs inputs;
    inputs.stage5.success = true;
    inputs.stage5.grid = buildStage4RegularGrid(3.0, 1.0);
    inputs.stage5.node_count = inputs.stage5.grid.nx * inputs.stage5.grid.ny;
    inputs.stage5.reliable_core_mask.assign(inputs.stage5.node_count, 1);
    inputs.stage5.reliable_core_node_count = inputs.stage5.node_count;

    inputs.stage7.success = true;
    inputs.stage7.grid = inputs.stage5.grid;
    inputs.stage7.node_count = inputs.stage5.node_count;
    inputs.stage7.reconstructed_mask.assign(inputs.stage7.node_count, 1);
    inputs.stage7.reliable_core_mask = inputs.stage5.reliable_core_mask;
    inputs.stage7.smooth_valid_mask.assign(inputs.stage7.node_count, 1);
    inputs.stage7.metric_domain_mask.assign(inputs.stage7.node_count, 1);
    inputs.stage7.metric_domain_node_count = inputs.stage7.node_count;
    inputs.stage7.smooth_valid_node_count = inputs.stage7.node_count;
    inputs.stage7.z_outer_smooth.assign(inputs.stage7.node_count, std::numeric_limits<double>::quiet_NaN());
    inputs.stage7.z_inner_smooth.assign(inputs.stage7.node_count, std::numeric_limits<double>::quiet_NaN());

    const double ao = 2.0;
    const double bo = -1.5;
    const double co = 0.5;
    const double do_ = 0.25;
    const double eo = -0.4;
    const double fo = 0.75;
    const double ai = -1.0;
    const double bi = 0.3;
    const double ci = -0.2;
    const double di = 0.1;
    const double ei = 0.15;
    const double fi = -0.05;

    for (std::size_t j = 0; j < inputs.stage7.grid.ny; ++j) {
        for (std::size_t i = 0; i < inputs.stage7.grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, inputs.stage7.grid.nx);
            const double x = inputs.stage7.grid.x_values[i];
            const double y = inputs.stage7.grid.y_values[j];
            inputs.stage7.z_outer_smooth[idx] = ao + (bo * x) + (co * y) + (do_ * x * x) + (eo * x * y) + (fo * y * y);
            inputs.stage7.z_inner_smooth[idx] = ai + (bi * x) + (ci * y) + (di * x * x) + (ei * x * y) + (fi * y * y);
        }
    }
    return inputs;
}

void testStage8ExactQuadraticRecovery() {
    const SyntheticStage8Inputs inputs = makeSyntheticStage8QuadraticInputs();
    FoldPatchAnalysisConfig config;
    config.stage8_export_csv = false;
    config.stage8_fit_radius = 2.5;
    config.stage8_min_points = 9;
    config.stage8_min_directional_span = 0.0;
    const auto stage8 =
        runGeometryAnalysisStage8DerivativeEstimation(inputs.stage7, inputs.stage5, config, nullptr);
    const std::size_t idx = nodeIndex(3, 3, inputs.stage7.grid.nx); // center x=0,y=0
    assertTrue(stage8.derivative_valid_mask[idx] == 1, "quadratic center node should be derivative-valid");
    assertTrue(near(stage8.outer_dz_dx[idx], -1.5, 1e-10), "outer dz/dx should recover exact quadratic slope");
    assertTrue(near(stage8.outer_dz_dy[idx], 0.5, 1e-10), "outer dz/dy should recover exact quadratic slope");
    assertTrue(near(stage8.outer_d2z_dx2[idx], 0.5, 1e-10), "outer d2z/dx2 should recover exact quadratic curvature");
    assertTrue(near(stage8.outer_d2z_dy2[idx], 1.5, 1e-10), "outer d2z/dy2 should recover exact quadratic curvature");
    assertTrue(near(stage8.outer_d2z_dxdy[idx], -0.4, 1e-10), "outer d2z/dxdy should recover exact mixed term");
    assertTrue(near(stage8.inner_dz_dx[idx], 0.3, 1e-10), "inner dz/dx should recover exact quadratic slope");
    assertTrue(near(stage8.inner_dz_dy[idx], -0.2, 1e-10), "inner dz/dy should recover exact quadratic slope");
    assertTrue(near(stage8.inner_d2z_dx2[idx], 0.2, 1e-10), "inner d2z/dx2 should recover exact quadratic curvature");
    assertTrue(near(stage8.inner_d2z_dy2[idx], -0.1, 1e-10), "inner d2z/dy2 should recover exact quadratic curvature");
    assertTrue(near(stage8.inner_d2z_dxdy[idx], 0.15, 1e-10), "inner d2z/dxdy should recover exact mixed term");
    assertTrue(stage8.outer_fit_rms_residual[idx] < 1e-10, "outer fit residual should be near zero for exact quadratic");
    assertTrue(stage8.inner_fit_rms_residual[idx] < 1e-10, "inner fit residual should be near zero for exact quadratic");
}

void testStage8InsufficientNeighborRejection() {
    const SyntheticStage8Inputs inputs = makeSyntheticStage8QuadraticInputs();
    FoldPatchAnalysisConfig config;
    config.stage8_export_csv = false;
    config.stage8_fit_radius = 0.01;
    config.stage8_min_points = 9;
    const auto stage8 =
        runGeometryAnalysisStage8DerivativeEstimation(inputs.stage7, inputs.stage5, config, nullptr);
    const std::size_t idx = nodeIndex(3, 3, inputs.stage7.grid.nx);
    assertTrue(stage8.derivative_valid_mask[idx] == 0, "tiny radius should reject derivative fit");
    assertTrue(stage8.derivative_invalid_insufficient_points_mask[idx] == 1,
               "insufficient-point reason should be flagged");
    assertTrue(std::isnan(stage8.outer_dz_dx[idx]), "invalid derivative output should remain NaN");
}

void testStage8BoundarySkewRejection() {
    SyntheticStage8Inputs inputs = makeSyntheticStage8QuadraticInputs();
    for (std::size_t j = 0; j < inputs.stage7.grid.ny; ++j) {
        for (std::size_t i = 0; i < inputs.stage7.grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, inputs.stage7.grid.nx);
            if (inputs.stage7.grid.x_values[i] < 0.0) {
                inputs.stage7.metric_domain_mask[idx] = 0;
                inputs.stage7.smooth_valid_mask[idx] = 0;
                inputs.stage5.reliable_core_mask[idx] = 0;
            }
        }
    }
    FoldPatchAnalysisConfig config;
    config.stage8_export_csv = false;
    config.stage8_fit_radius = 2.5;
    config.stage8_min_points = 6;
    const auto stage8 =
        runGeometryAnalysisStage8DerivativeEstimation(inputs.stage7, inputs.stage5, config, nullptr);
    const std::size_t idx = nodeIndex(3, 3, inputs.stage7.grid.nx);
    assertTrue(stage8.derivative_valid_mask[idx] == 0, "one-sided neighborhood should be rejected");
    assertTrue(stage8.derivative_invalid_boundary_neighbor_geometry_mask[idx] == 1,
               "boundary-neighbor-geometry reason should be flagged");
}

void testStage8RankDeficientOrPoorConditioningRejection() {
    SyntheticStage8Inputs inputs = makeSyntheticStage8QuadraticInputs();
    for (std::size_t j = 0; j < inputs.stage7.grid.ny; ++j) {
        for (std::size_t i = 0; i < inputs.stage7.grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, inputs.stage7.grid.nx);
            if (j != i) {
                inputs.stage7.metric_domain_mask[idx] = 0;
                inputs.stage7.smooth_valid_mask[idx] = 0;
                inputs.stage5.reliable_core_mask[idx] = 0;
            }
        }
    }
    FoldPatchAnalysisConfig config;
    config.stage8_export_csv = false;
    config.stage8_fit_radius = 3.0;
    config.stage8_min_points = 6;
    config.stage8_require_centered_support = false;
    config.stage8_min_directional_span = 0.5;
    const auto stage8 =
        runGeometryAnalysisStage8DerivativeEstimation(inputs.stage7, inputs.stage5, config, nullptr);
    const std::size_t idx = nodeIndex(3, 3, inputs.stage7.grid.nx);
    assertTrue(stage8.derivative_valid_mask[idx] == 0, "rank-deficient neighborhood should not validate");
    assertTrue(stage8.derivative_invalid_rank_deficient_mask[idx] == 1 ||
                   stage8.derivative_invalid_poor_conditioning_mask[idx] == 1 ||
                   stage8.derivative_invalid_high_residual_mask[idx] == 1 ||
                   stage8.derivative_invalid_boundary_neighbor_geometry_mask[idx] == 1 ||
                   stage8.derivative_invalid_insufficient_points_mask[idx] == 1,
               "pathological support should report a deterministic instability-related rejection reason");
}

void testStage8CsvExportSmoke() {
    const SyntheticStage8Inputs inputs = makeSyntheticStage8QuadraticInputs();
    FoldPatchAnalysisConfig config;
    config.output_prefix = "stage8_smoke";
    config.stage8_export_csv = true;
    config.stage8_fit_radius = 2.5;
    config.stage8_min_points = 9;
    const auto stage8 =
        runGeometryAnalysisStage8DerivativeEstimation(inputs.stage7, inputs.stage5, config, nullptr);
    assertTrue(std::filesystem::exists(stage8.outer_derivatives_csv_path), "outer derivative csv should be exported");
    assertTrue(std::filesystem::exists(stage8.inner_derivatives_csv_path), "inner derivative csv should be exported");
    assertTrue(std::filesystem::exists(stage8.derivative_valid_mask_csv_path), "derivative-valid mask csv should exist");
    assertTrue(std::filesystem::exists(stage8.derivative_failure_reason_csv_path), "failure-reason csv should exist");
    assertTrue(std::filesystem::exists(stage8.derivative_summary_csv_path), "stage8 summary csv should exist");
    assertTrue(stage8.outer_derivatives_csv_path.find("_outer_derivatives.csv") != std::string::npos,
               "outer derivative filename should follow required suffix");
    removeIfExists(stage8.outer_derivatives_csv_path);
    removeIfExists(stage8.inner_derivatives_csv_path);
    removeIfExists(stage8.derivative_valid_mask_csv_path);
    removeIfExists(stage8.derivative_failure_reason_csv_path);
    removeIfExists(stage8.derivative_summary_csv_path);
}

GeometryStage8DerivativeEstimationResult makeSyntheticStage8ResultForStage9(std::size_t nx = 3, std::size_t ny = 3) {
    GeometryStage8DerivativeEstimationResult stage8;
    stage8.success = true;
    stage8.grid.nx = nx;
    stage8.grid.ny = ny;
    stage8.grid.spacing = 1.0;
    stage8.grid.x_min = -1.0;
    stage8.grid.x_max = 1.0;
    stage8.grid.y_min = -1.0;
    stage8.grid.y_max = 1.0;
    stage8.grid.x_values.resize(nx, 0.0);
    stage8.grid.y_values.resize(ny, 0.0);
    for (std::size_t i = 0; i < nx; ++i) {
        stage8.grid.x_values[i] = static_cast<double>(i) - static_cast<double>(nx / 2);
    }
    for (std::size_t j = 0; j < ny; ++j) {
        stage8.grid.y_values[j] = static_cast<double>(j) - static_cast<double>(ny / 2);
    }
    stage8.node_count = nx * ny;
    stage8.reconstructed_mask.assign(stage8.node_count, 1);
    stage8.reliable_core_mask.assign(stage8.node_count, 1);
    stage8.metric_domain_mask.assign(stage8.node_count, 1);
    stage8.derivative_valid_mask.assign(stage8.node_count, 1);
    stage8.derivative_valid_node_count = stage8.node_count;
    stage8.metric_domain_node_count = stage8.node_count;

    stage8.outer_dz_dx.assign(stage8.node_count, 0.0);
    stage8.outer_dz_dy.assign(stage8.node_count, 0.0);
    stage8.outer_d2z_dx2.assign(stage8.node_count, 0.0);
    stage8.outer_d2z_dy2.assign(stage8.node_count, 0.0);
    stage8.outer_d2z_dxdy.assign(stage8.node_count, 0.0);
    stage8.inner_dz_dx.assign(stage8.node_count, 0.0);
    stage8.inner_dz_dy.assign(stage8.node_count, 0.0);
    stage8.inner_d2z_dx2.assign(stage8.node_count, 0.0);
    stage8.inner_d2z_dy2.assign(stage8.node_count, 0.0);
    stage8.inner_d2z_dxdy.assign(stage8.node_count, 0.0);
    return stage8;
}

GeometryStage7SmoothedSurfaceResult makeSyntheticStage7ResultForStage9(const GeometryStage8DerivativeEstimationResult& stage8) {
    GeometryStage7SmoothedSurfaceResult stage7;
    stage7.success = true;
    stage7.grid = stage8.grid;
    stage7.node_count = stage8.node_count;
    stage7.reconstructed_mask = stage8.reconstructed_mask;
    stage7.reliable_core_mask = stage8.reliable_core_mask;
    stage7.metric_domain_mask = stage8.metric_domain_mask;
    stage7.smooth_valid_mask.assign(stage7.node_count, 1);
    stage7.z_outer_smooth.assign(stage7.node_count, 0.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 0.0);
    return stage7;
}

GeometryStage9CurvatureComputationResult makeSyntheticStage9ResultForStage10(
    const GeometryStage7SmoothedSurfaceResult& stage7) {
    GeometryStage9CurvatureComputationResult stage9;
    stage9.success = true;
    stage9.grid = stage7.grid;
    stage9.node_count = stage7.node_count;
    stage9.reconstructed_mask = stage7.reconstructed_mask;
    stage9.reliable_core_mask = stage7.reliable_core_mask;
    stage9.metric_domain_mask = stage7.metric_domain_mask;
    stage9.derivative_valid_mask.assign(stage9.node_count, 1);
    stage9.curvature_valid_mask.assign(stage9.node_count, 1);
    stage9.metric_domain_node_count = stage9.node_count;
    stage9.derivative_valid_node_count = stage9.node_count;
    stage9.curvature_valid_node_count = stage9.node_count;
    return stage9;
}

void testStage10SimpleValidThicknessField() {
    auto stage8 = makeSyntheticStage8ResultForStage9(5, 5);
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "Stage 10 should succeed for positive constant vertical thickness");
    assertTrue(stage10.thickness_valid_node_count == stage10.node_count, "all nodes should be valid");
    for (std::size_t idx = 0; idx < stage10.node_count; ++idx) {
        assertTrue(near(stage10.thickness_vertical[idx], 3.0, 1e-12), "thickness should be outer-inner at each node");
        assertTrue(stage10.thickness_valid_mask[idx] == 1, "all constant-thickness nodes should be valid");
    }
    assertTrue(near(stage10.mean_thickness_vertical, 3.0, 1e-12), "mean thickness should be 3");
    assertTrue(near(stage10.median_thickness_vertical, 3.0, 1e-12), "median thickness should be 3");
    assertTrue(near(stage10.min_thickness_vertical, 3.0, 1e-12), "min thickness should be 3");
    assertTrue(near(stage10.max_thickness_vertical, 3.0, 1e-12), "max thickness should be 3");
    assertTrue(stage10.thickness_method_label == "stage10_vertical_difference_from_stage7_smooth_surfaces",
               "Stage 10 method metadata should match v1 contract");
    assertTrue(stage10.local_thickness_definition == "stage7_z_outer_smooth_minus_stage7_z_inner_smooth",
               "Stage 10 local thickness definition should explicitly reference Stage 7 smooth surfaces");
    assertTrue(stage10.thickness_input_surface_definition == "stage7_smoothed_outer_inner_graph_surfaces",
               "Stage 10 input-surface definition should explicitly reference Stage 7 smooth graph surfaces");
    assertTrue(stage10.thickness_contract_note.find("stage10_v1") != std::string::npos &&
                   stage10.thickness_contract_note.find("stage7_smooth_vertical_separation") != std::string::npos &&
                   stage10.thickness_contract_note.find("repackaged_as_a_dedicated_stage10_field") != std::string::npos,
               "Stage 10 contract note should explicitly describe the Stage 10 v1 repackaged semantics");
}

void testStage10OutsideDomainInvalidation() {
    auto stage8 = makeSyntheticStage8ResultForStage9(5, 5);
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    stage7.metric_domain_mask.assign(stage7.node_count, 0);
    stage7.metric_domain_mask[nodeIndex(1, 1, stage7.grid.nx)] = 1;
    stage7.metric_domain_mask[nodeIndex(2, 1, stage7.grid.nx)] = 1;
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "Stage 10 should still succeed when only a subset is in domain");
    assertTrue(stage10.thickness_attempted_node_count == 2, "only metric-domain nodes should be attempted");
    assertTrue(stage10.thickness_valid_node_count == 2, "only in-domain nodes should validate");
    assertTrue(stage10.thickness_excluded_outside_domain_count == stage10.node_count - 2,
               "outside-domain count should match excluded nodes");
    assertTrue(stage10.thickness_attempted_invalid_node_count == 0,
               "outside-domain exclusions must not increment attempted-invalid count");
    for (std::size_t idx = 0; idx < stage10.node_count; ++idx) {
        if (stage7.metric_domain_mask[idx] == 0) {
            assertTrue(stage10.thickness_attempted_mask[idx] == 0, "outside-domain nodes should not be attempted");
            assertTrue(stage10.thickness_invalid_outside_domain_mask[idx] == 1, "outside-domain reason should be set");
        }
    }
}

void testStage10NegativeOrZeroInvalidation() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    const std::size_t zero_idx = nodeIndex(0, 0, stage7.grid.nx);
    const std::size_t negative_idx = nodeIndex(0, 1, stage7.grid.nx);
    stage7.z_inner_smooth[zero_idx] = stage7.z_outer_smooth[zero_idx];
    stage7.z_inner_smooth[negative_idx] = stage7.z_outer_smooth[negative_idx] + 1.0;
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.thickness_invalid_negative_or_zero_mask[zero_idx] == 1, "zero thickness should be invalid");
    assertTrue(stage10.thickness_invalid_negative_or_zero_mask[negative_idx] == 1, "negative thickness should be invalid");
    assertTrue(stage10.thickness_invalid_negative_or_zero_count == 2, "negative/zero counter should match");
    assertTrue(stage10.thickness_attempted_invalid_node_count == 2,
               "attempted-invalid count should match negative-or-zero invalid nodes");
    assertTrue(stage10.thickness_excluded_outside_domain_count == 0,
               "in-domain failed nodes should not be counted as excluded");
    assertTrue(stage10.thickness_valid_mask[zero_idx] == 0 && stage10.thickness_valid_mask[negative_idx] == 0,
               "nonphysical nodes should not be valid");
    assertTrue(stage10.thickness_valid_node_count == stage10.node_count - 2, "valid count should exclude invalid nodes");
    assertTrue(stage10.min_thickness_vertical > 0.0, "summary stats should only consider positive valid nodes");
}

void testStage10ThresholdInvalidation() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    stage7.z_inner_smooth[nodeIndex(0, 0, stage7.grid.nx)] = 9.2; // 0.8 below min
    stage7.z_inner_smooth[nodeIndex(0, 1, stage7.grid.nx)] = 2.0; // 8.0 above max
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;
    config.stage10_min_thickness = 1.0;
    config.stage10_max_thickness = 5.0;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.thickness_invalid_below_min_threshold_count == 1, "below-min threshold count should be tracked");
    assertTrue(stage10.thickness_invalid_above_max_threshold_count == 1, "above-max threshold count should be tracked");
    assertTrue(stage10.thickness_invalid_below_min_threshold_mask[nodeIndex(0, 0, stage7.grid.nx)] == 1,
               "below-min mask should be set");
    assertTrue(stage10.thickness_invalid_above_max_threshold_mask[nodeIndex(0, 1, stage7.grid.nx)] == 1,
               "above-max mask should be set");
}

void testStage10CurvatureValidRestrictionMode() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    stage9.curvature_valid_mask.assign(stage9.node_count, 0);
    stage9.curvature_valid_mask[nodeIndex(1, 1, stage7.grid.nx)] = 1;
    stage9.curvature_valid_mask[nodeIndex(2, 1, stage7.grid.nx)] = 1;
    stage9.curvature_valid_node_count = 2;
    FoldPatchAnalysisConfig unrestricted;
    unrestricted.stage10_export_csv = false;
    FoldPatchAnalysisConfig restricted = unrestricted;
    restricted.stage10_use_curvature_valid_domain_only = true;

    const auto stage10_unrestricted =
        runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, unrestricted, nullptr);
    const auto stage10_restricted =
        runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, restricted, nullptr);

    assertTrue(stage10_unrestricted.thickness_attempted_node_count == stage7.node_count,
               "unrestricted mode should attempt full metric domain");
    assertTrue(stage10_restricted.thickness_attempted_node_count == 2,
               "restricted mode should attempt only curvature-valid subset");
    assertTrue(stage10_unrestricted.thickness_valid_and_curvature_valid_node_count == 2,
               "unrestricted mode should track valid-curvature intersection cardinality");
    assertTrue(near(stage10_unrestricted.thickness_valid_intersection_fraction_of_metric_domain,
                    2.0 / static_cast<double>(stage7.node_count),
                    1e-12),
               "intersection fraction should be normalized by metric-domain node count");
    assertTrue(stage10_unrestricted.thickness_valid_intersection_fraction_of_metric_domain <= 1.0,
               "intersection fraction should always remain in [0,1]");
    assertTrue(stage10_unrestricted.thickness_domain_definition == "stage7_metric_domain",
               "unrestricted label should match contract");
    assertTrue(stage10_restricted.thickness_domain_definition == "stage7_metric_domain_intersect_stage9_curvature_valid",
               "restricted label should match contract");
}

void testStage10SummaryCsvUsesCorrectedContractFields() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    stage9.curvature_valid_mask.assign(stage9.node_count, 0);
    stage9.curvature_valid_mask[nodeIndex(1, 1, stage7.grid.nx)] = 1;
    stage9.curvature_valid_mask[nodeIndex(2, 1, stage7.grid.nx)] = 1;
    stage9.curvature_valid_node_count = 2;

    FoldPatchAnalysisConfig config;
    config.output_prefix = "stage10_contract_fields_test";
    config.stage10_export_csv = true;
    removeIfExists(config.output_prefix + "_thickness_vertical_summary.csv");

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "Stage 10 should succeed while exporting summary CSV");
    const auto rows = readCsvRows(stage10.thickness_summary_csv_path);
    assertTrue(rows.size() >= 2, "Stage 10 summary CSV should include header and data rows");
    const std::string header_line = [&]() {
        std::ifstream in(stage10.thickness_summary_csv_path);
        std::string line;
        std::getline(in, line);
        return line;
    }();
    assertTrue(header_line.find("thickness_excluded_outside_domain_count") != std::string::npos,
               "summary CSV should include excluded outside-domain count");
    assertTrue(header_line.find("thickness_attempted_invalid_node_count") != std::string::npos,
               "summary CSV should include attempted-invalid count");
    assertTrue(header_line.find("thickness_valid_and_curvature_valid_node_count") != std::string::npos,
               "summary CSV should include valid-curvature intersection count");
    assertTrue(header_line.find("thickness_valid_intersection_fraction_of_metric_domain") != std::string::npos,
               "summary CSV should include corrected intersection fraction field");
    assertTrue(header_line.find("thickness_invalid_node_count") == std::string::npos,
               "summary CSV should remove legacy aggregate invalid field");
    assertTrue(header_line.find("thickness_invalid_outside_domain_count") == std::string::npos,
               "summary CSV should remove legacy outside-domain invalid field");
    assertTrue(header_line.find("thickness_valid_fraction_of_curvature_valid") == std::string::npos,
               "summary CSV should remove legacy invalidly-normalized fraction field");
    removeIfExists(stage10.thickness_vertical_csv_path);
    removeIfExists(stage10.thickness_valid_mask_csv_path);
    removeIfExists(stage10.thickness_invalid_reason_csv_path);
    removeIfExists(stage10.thickness_summary_csv_path);
}

void testStage10Stage8FailureDoesNotBlockVerticalThickness() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    stage8.success = false;
    stage8.derivative_valid_mask.clear();
    auto stage7 = makeSyntheticStage7ResultForStage9(makeSyntheticStage8ResultForStage9());
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 8.5);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "Stage 10 should run even if Stage 8 failed");
    assertTrue(stage10.thickness_valid_node_count == stage7.node_count, "Stage 10 should compute from Stage 7 surfaces");
    const bool has_stage8_notice =
        std::any_of(stage10.messages.begin(), stage10.messages.end(), [](const std::string& msg) {
            return msg.find("Stage 8 derivative estimation unavailable") != std::string::npos;
        });
    assertTrue(has_stage8_notice, "Stage 10 should report Stage 8 unavailability in messages");
}

void testStage10RadialFlatParallelPlanesDiffersFromVerticalOffAxis() {
    auto stage8 = makeSyntheticStage8ResultForStage9(5, 5);
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;
    config.stage10_thickness_method = FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "Stage 10 radial should succeed for flat parallel surfaces");
    const std::size_t off_axis = nodeIndex(3, 2, stage7.grid.nx);
    assertTrue(stage10.thickness_radial[off_axis] > stage10.thickness_vertical[off_axis],
               "off-axis radial thickness should differ from vertical for flat planes");
    assertTrue(stage10.thickness_method_label == "stage10_radial_origin_centered_ray_intersection",
               "radial method label should be explicit");
    assertTrue(near(stage10.mean_thickness_active, stage10.mean_thickness_radial, 1e-12),
               "active mean should map to radial mean in radial mode");
    assertTrue(near(stage10.median_thickness_active, stage10.median_thickness_radial, 1e-12),
               "active median should map to radial median in radial mode");
    assertTrue(std::isfinite(stage10.mean_thickness_active), "active mean should be finite when radial has valid nodes");
}

void testStage10RadialConcentricSphericalCapsApproxConstantGap() {
    auto stage8 = makeSyntheticStage8ResultForStage9(5, 5);
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    const double r_outer = 10.0;
    const double r_inner = 8.0;
    for (std::size_t j = 0; j < stage7.grid.ny; ++j) {
        for (std::size_t i = 0; i < stage7.grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, stage7.grid.nx);
            const double x = stage7.grid.x_values[i];
            const double y = stage7.grid.y_values[j];
            const double r2 = x * x + y * y;
            stage7.z_outer_smooth[idx] = std::sqrt(std::max(0.0, r_outer * r_outer - r2));
            stage7.z_inner_smooth[idx] = std::sqrt(std::max(0.0, r_inner * r_inner - r2));
        }
    }
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;
    config.stage10_thickness_method = FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "Stage 10 radial should succeed for spherical-cap synthetic case");
    assertTrue(std::fabs(stage10.mean_thickness_radial - (r_outer - r_inner)) < 0.1,
               "radial mean thickness should be close to spherical shell gap");
}

void testStage10RadialNoBracketFailureIsTracked() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, -1.0);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;
    config.stage10_thickness_method = FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(!stage10.success, "radial no-bracket everywhere should fail overall");
    assertTrue(stage10.thickness_radial_invalid_no_bracket_count > 0, "no-bracket count should be tracked");
    const bool has_sampling_loss_dominance_note =
        std::any_of(stage10.messages.begin(), stage10.messages.end(), [](const std::string& msg) {
            return msg.find("dominated by inner-surface sampling-domain loss") != std::string::npos;
        });
    assertTrue(!has_sampling_loss_dominance_note,
               "sampling-domain-loss dominance note should not appear when no-bracket failures exist");
}

void testStage10RadialPOutInCsvExport() {
    auto stage8 = makeSyntheticStage8ResultForStage9(5, 5);
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_thickness_method = FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;
    config.stage10_export_csv = true;
    config.output_prefix = "stage10_radial_points_csv";
    removeIfExists(config.output_prefix + "_thickness_radial_P_out_P_in_t.csv");

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "radial Stage 10 should succeed for P_out/P_in export test");
    assertTrue(!stage10.thickness_radial_p_out_in_csv_path.empty(), "radial P_out/P_in CSV path should be populated");
    assertTrue(std::filesystem::exists(stage10.thickness_radial_p_out_in_csv_path),
               "radial P_out/P_in CSV should be created");
    std::ifstream in(stage10.thickness_radial_p_out_in_csv_path);
    std::string header;
    std::getline(in, header);
    assertTrue(header.find("P_out x") != std::string::npos && header.find("P_in z") != std::string::npos &&
                   header.find("t_radial") != std::string::npos,
               "radial P_out/P_in CSV header should include required columns");
    const bool has_export_definition_message =
        std::any_of(stage10.messages.begin(), stage10.messages.end(), [](const std::string& msg) {
            return msg.find("Stage10 radial P_out/P_in export is only defined") != std::string::npos;
        });
    assertTrue(has_export_definition_message, "radial result should describe P_out/P_in export validity domain");

    removeIfExists(stage10.thickness_vertical_csv_path);
    removeIfExists(stage10.thickness_valid_mask_csv_path);
    removeIfExists(stage10.thickness_invalid_reason_csv_path);
    removeIfExists(stage10.thickness_summary_csv_path);
    removeIfExists(stage10.thickness_radial_p_out_in_csv_path);
}


void testStage10RadialInnerDomainLimitationMetadataAndFractions() {
    auto stage8 = makeSyntheticStage8ResultForStage9(5, 5);
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;
    config.stage10_thickness_method = FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;

    auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "radial stage should succeed for baseline planar case");

    for (std::size_t j = 0; j < stage7.grid.ny; ++j) {
        for (std::size_t i = 0; i < stage7.grid.nx; ++i) {
            if (i == 0 || j == 0 || i + 1 == stage7.grid.nx || j + 1 == stage7.grid.ny) {
                const std::size_t idx = nodeIndex(i, j, stage7.grid.nx);
                stage10.thickness_valid_mask[idx] = 0;
                stage10.thickness_radial_valid_mask[idx] = 0;
                stage10.thickness_radial_invalid_outside_inner_domain_mask[idx] = 1;
            }
        }
    }
    stage10.thickness_radial_valid_node_count = 9;
    stage10.thickness_valid_node_count = 9;
    stage10.thickness_radial_invalid_outside_inner_domain_count = 16;
    stage10.thickness_radial_invalid_no_bracket_count = 0;
    stage10.thickness_radial_invalid_root_failure_count = 0;
    const std::size_t metric_domain_count = stage10.node_count;
    stage10.thickness_radial_valid_fraction_of_metric_domain =
        static_cast<double>(stage10.thickness_radial_valid_node_count) / static_cast<double>(metric_domain_count);
    stage10.thickness_radial_invalid_outside_inner_domain_fraction_of_metric_domain =
        static_cast<double>(stage10.thickness_radial_invalid_outside_inner_domain_count) / static_cast<double>(metric_domain_count);

    assertTrue(!stage10.thickness_domain_limitation_note.empty(),
               "radial mode should provide explicit domain limitation note");
    assertTrue(stage10.thickness_domain_limitation_note.find("inner_surface_sampleability") != std::string::npos,
               "limitation note should reference inner-surface sampleability");
    assertTrue(!stage10.radial_invalid_outside_inner_domain_definition.empty(),
               "radial mode should define outside-inner-domain failures");
    assertTrue(stage10.radial_invalid_outside_inner_domain_definition.find("radial_ray") != std::string::npos,
               "outside-inner-domain definition should mention radial ray");
    assertTrue(near(stage10.thickness_radial_valid_fraction_of_metric_domain, 9.0 / 25.0, 1e-12),
               "radial valid fraction should be computed from metric-domain denominator");
    assertTrue(near(stage10.thickness_radial_invalid_outside_inner_domain_fraction_of_metric_domain, 16.0 / 25.0, 1e-12),
               "radial outside-inner-domain fraction should be computed from metric-domain denominator");
}

void testStage10RadialSamplingLossNoteEmissionGuard() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_export_csv = false;
    config.stage10_thickness_method = FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;

    auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "baseline radial run should succeed");
    stage10.thickness_radial_invalid_outside_inner_domain_count = 2;
    stage10.thickness_radial_invalid_no_bracket_count = 1;
    stage10.thickness_radial_invalid_root_failure_count = 0;

    const bool should_emit_sampling_loss_note = stage10.thickness_radial_invalid_outside_inner_domain_count > 0 &&
                                                stage10.thickness_radial_invalid_no_bracket_count == 0 &&
                                                stage10.thickness_radial_invalid_root_failure_count == 0;
    assertTrue(!should_emit_sampling_loss_note,
               "sampling-domain-loss dominance note must not be emitted when no-bracket failures are present");
}

void testStage10RadialSummaryCsvSchemaUsesActiveAndLimitationFields() {
    auto stage8 = makeSyntheticStage8ResultForStage9(5, 5);
    auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    stage7.z_outer_smooth.assign(stage7.node_count, 10.0);
    stage7.z_inner_smooth.assign(stage7.node_count, 7.0);
    const auto stage9 = makeSyntheticStage9ResultForStage10(stage7);
    FoldPatchAnalysisConfig config;
    config.stage10_thickness_method = FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;
    config.stage10_export_csv = true;
    config.output_prefix = "stage10_radial_summary_semantics_test";
    removeIfExists(config.output_prefix + "_thickness_radial_summary.csv");

    const auto stage10 = runGeometryAnalysisStage10ThicknessComputation(stage7, stage8, stage9, config, nullptr);
    assertTrue(stage10.success, "radial Stage 10 should succeed while exporting summary CSV");
    std::ifstream in(stage10.thickness_summary_csv_path);
    std::string header;
    std::getline(in, header);
    assertTrue(header.find("stage10_method") != std::string::npos, "radial summary should include stage10_method");
    assertTrue(header.find("thickness_domain_limitation_note") != std::string::npos,
               "radial summary should include limitation note field");
    assertTrue(header.find("radial_invalid_outside_inner_domain_definition") != std::string::npos,
               "radial summary should include explicit outside-inner-domain definition field");
    assertTrue(header.find("mean_thickness_active") != std::string::npos,
               "radial summary should include active-method mean field");
    assertTrue(header.find("thickness_radial_valid_fraction_of_metric_domain") != std::string::npos,
               "radial summary should include radial valid fraction field");
    assertTrue(header.find("thickness_radial_invalid_outside_inner_domain_fraction_of_metric_domain") != std::string::npos,
               "radial summary should include radial outside-inner-domain fraction field");
    assertTrue(header.find("mean_thickness_vertical") == std::string::npos,
               "radial summary should not use vertical-named fields as active headline metrics");

    removeIfExists(stage10.thickness_vertical_csv_path);
    removeIfExists(stage10.thickness_valid_mask_csv_path);
    removeIfExists(stage10.thickness_invalid_reason_csv_path);
    removeIfExists(stage10.thickness_summary_csv_path);
    removeIfExists(stage10.thickness_radial_p_out_in_csv_path);
}

void testStage9PlaneCurvature() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    const auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    FoldPatchAnalysisConfig config;
    config.stage9_export_csv = false;

    const auto stage9 = runGeometryAnalysisStage9CurvatureComputation(stage7, stage8, config, nullptr);
    assertTrue(stage9.success, "Stage 9 plane computation should succeed");
    assertTrue(stage9.curvature_valid_node_count == stage8.node_count, "all derivative-valid plane nodes should be curvature-valid");
    for (std::size_t idx = 0; idx < stage9.node_count; ++idx) {
        assertTrue(near(stage9.outer_mean_curvature_H[idx], 0.0, 1e-12), "plane outer H should be zero");
        assertTrue(near(stage9.outer_oriented_mean_curvature_H[idx], 0.0, 1e-12), "plane outer oriented H should be zero");
        assertTrue(near(stage9.outer_gaussian_curvature_K[idx], 0.0, 1e-12), "plane outer K should be zero");
        assertTrue(near(stage9.inner_mean_curvature_H[idx], 0.0, 1e-12), "plane inner H should be zero");
        assertTrue(near(stage9.inner_oriented_mean_curvature_H[idx], 0.0, 1e-12), "plane inner oriented H should be zero");
        assertTrue(near(stage9.inner_gaussian_curvature_K[idx], 0.0, 1e-12), "plane inner K should be zero");
        assertTrue(std::fabs(stage9.outer_graph_normal_dot_radial[idx]) <= 1.0 + 1e-12,
                   "plane outer graph_normal_dot_radial should be normalized to [-1,1]");
        assertTrue(std::fabs(stage9.inner_graph_normal_dot_radial[idx]) <= 1.0 + 1e-12,
                   "plane inner graph_normal_dot_radial should be normalized to [-1,1]");
    }
}

void testStage9ParaboloidCurvatureAtOrigin() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    const double alpha = 0.3;
    const double beta = 0.7;
    for (std::size_t j = 0; j < stage8.grid.ny; ++j) {
        for (std::size_t i = 0; i < stage8.grid.nx; ++i) {
            const std::size_t idx = nodeIndex(i, j, stage8.grid.nx);
            const double x = stage8.grid.x_values[i];
            const double y = stage8.grid.y_values[j];
            stage8.outer_dz_dx[idx] = 2.0 * alpha * x;
            stage8.outer_dz_dy[idx] = 2.0 * beta * y;
            stage8.outer_d2z_dx2[idx] = 2.0 * alpha;
            stage8.outer_d2z_dy2[idx] = 2.0 * beta;
            stage8.outer_d2z_dxdy[idx] = 0.0;
            stage8.inner_dz_dx[idx] = stage8.outer_dz_dx[idx];
            stage8.inner_dz_dy[idx] = stage8.outer_dz_dy[idx];
            stage8.inner_d2z_dx2[idx] = stage8.outer_d2z_dx2[idx];
            stage8.inner_d2z_dy2[idx] = stage8.outer_d2z_dy2[idx];
            stage8.inner_d2z_dxdy[idx] = 0.0;
        }
    }
    const auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    FoldPatchAnalysisConfig config;
    config.stage9_export_csv = false;
    const auto stage9 = runGeometryAnalysisStage9CurvatureComputation(stage7, stage8, config, nullptr);

    const std::size_t center = nodeIndex(1, 1, stage8.grid.nx);
    const double expected_H = alpha + beta;
    const double expected_K = 4.0 * alpha * beta;
    assertTrue(near(stage9.outer_mean_curvature_H[center], expected_H, 1e-12), "paraboloid outer H at origin should match analytic");
    assertTrue(near(stage9.outer_gaussian_curvature_K[center], expected_K, 1e-12),
               "paraboloid outer K at origin should match analytic");
    assertTrue(near(stage9.outer_oriented_mean_curvature_H[center], -expected_H, 1e-12),
               "paraboloid outer oriented H should flip when graph normal is not outward");
    assertTrue(stage9.outer_orientation_flip_applied_flag[center] == 1, "outer center should apply orientation flip");
    assertTrue(stage9.outer_outward_normal_alignment_flag[center] == 0, "outer center graph normal should not be outward");
    assertTrue(near(stage9.inner_mean_curvature_H[center], expected_H, 1e-12), "paraboloid inner H at origin should match analytic");
    assertTrue(near(stage9.inner_gaussian_curvature_K[center], expected_K, 1e-12),
               "paraboloid inner K at origin should match analytic");
    assertTrue(near(stage9.inner_oriented_mean_curvature_H[center], -expected_H, 1e-12),
               "paraboloid inner oriented H should flip when graph normal is not outward");
    assertTrue(stage9.inner_orientation_flip_applied_flag[center] == 1, "inner center should apply orientation flip");
    assertTrue(stage9.inner_outward_normal_alignment_flag[center] == 0, "inner center graph normal should not be outward");
}

void testStage9InvalidInputPropagation() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    const std::size_t idx = nodeIndex(0, 0, stage8.grid.nx);
    stage8.outer_dz_dx[idx] = std::numeric_limits<double>::infinity();
    const auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    FoldPatchAnalysisConfig config;
    config.stage9_export_csv = false;
    const auto stage9 = runGeometryAnalysisStage9CurvatureComputation(stage7, stage8, config, nullptr);

    assertTrue(stage9.curvature_valid_mask[idx] == 0, "non-finite derivative input should invalidate curvature node");
    assertTrue(stage9.curvature_invalid_nonfinite_input_mask[idx] == 1, "non-finite input mask should be flagged");
    assertTrue(std::isnan(stage9.outer_mean_curvature_H[idx]), "invalid node outer H should remain NaN");
    assertTrue(std::isnan(stage9.outer_oriented_mean_curvature_H[idx]), "invalid node outer oriented H should remain NaN");
    assertTrue(std::isnan(stage9.outer_gaussian_curvature_K[idx]), "invalid node outer K should remain NaN");
    assertTrue(std::isnan(stage9.inner_mean_curvature_H[idx]), "invalid node inner H should remain NaN");
    assertTrue(std::isnan(stage9.inner_oriented_mean_curvature_H[idx]), "invalid node inner oriented H should remain NaN");
    assertTrue(std::isnan(stage9.inner_gaussian_curvature_K[idx]), "invalid node inner K should remain NaN");
    assertTrue(stage9.curvature_invalid_nonfinite_input_count == 1, "invalid non-finite input count should increment");
}

void testStage9NonFiniteOutputGuard() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    const std::size_t idx = nodeIndex(0, 1, stage8.grid.nx);
    stage8.outer_dz_dx[idx] = 1e308;
    stage8.outer_dz_dy[idx] = 1e308;
    stage8.outer_d2z_dx2[idx] = 1e308;
    stage8.outer_d2z_dy2[idx] = 1e308;
    stage8.outer_d2z_dxdy[idx] = 1e308;
    stage8.inner_dz_dx[idx] = 1e308;
    stage8.inner_dz_dy[idx] = 1e308;
    stage8.inner_d2z_dx2[idx] = 1e308;
    stage8.inner_d2z_dy2[idx] = 1e308;
    stage8.inner_d2z_dxdy[idx] = 1e308;
    const auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    FoldPatchAnalysisConfig config;
    config.stage9_export_csv = false;
    const auto stage9 = runGeometryAnalysisStage9CurvatureComputation(stage7, stage8, config, nullptr);

    assertTrue(stage9.curvature_valid_mask[idx] == 0, "overflow-like derivatives should not produce valid curvature");
    assertTrue(stage9.curvature_invalid_nonfinite_output_mask[idx] == 1, "non-finite output mask should be flagged");
    assertTrue(stage9.curvature_invalid_nonfinite_output_count == 1, "non-finite output counter should increment");
}

void testStage9CsvExportSmoke() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    const auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    FoldPatchAnalysisConfig config;
    config.output_prefix = "stage9_smoke";
    config.stage9_export_csv = true;
    const auto stage9 = runGeometryAnalysisStage9CurvatureComputation(stage7, stage8, config, nullptr);

    assertTrue(std::filesystem::exists(stage9.outer_curvature_csv_path), "outer curvature csv should be exported");
    assertTrue(std::filesystem::exists(stage9.inner_curvature_csv_path), "inner curvature csv should be exported");
    assertTrue(std::filesystem::exists(stage9.curvature_valid_mask_csv_path), "curvature-valid mask csv should be exported");
    assertTrue(std::filesystem::exists(stage9.curvature_summary_csv_path), "curvature summary csv should be exported");
    removeIfExists(stage9.outer_curvature_csv_path);
    removeIfExists(stage9.inner_curvature_csv_path);
    removeIfExists(stage9.curvature_valid_mask_csv_path);
    removeIfExists(stage9.curvature_summary_csv_path);
}

void testStage9QcWarnFlagsAndConfidenceClass() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    const std::size_t center = nodeIndex(1, 1, stage8.grid.nx);
    stage8.outer_d2z_dx2[center] = 10.0;
    stage8.outer_d2z_dy2[center] = 10.0;
    stage8.inner_d2z_dx2[center] = 10.0;
    stage8.inner_d2z_dy2[center] = 10.0;

    const auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    FoldPatchAnalysisConfig config;
    config.stage9_export_csv = false;
    config.stage9_qc_n_tail = 4.0;
    config.stage9_qc_n_spike = 2.0;
    config.stage9_qc_abs_scale_floor = 1e-6;
    const auto stage9 = runGeometryAnalysisStage9CurvatureComputation(stage7, stage8, config, nullptr);

    assertTrue(stage9.outer_K_local_spike_flag[center] == 1, "outer center node should be local-spike flagged");
    assertTrue(stage9.inner_K_local_spike_flag[center] == 1, "inner center node should be local-spike flagged");
    assertTrue(stage9.outer_K_qc_warn_flag[center] == 1, "outer center node should be warned");
    assertTrue(stage9.inner_K_qc_warn_flag[center] == 1, "inner center node should be warned");
    assertTrue(stage9.outer_K_confidence_class[center] == 1, "warned outer node should have confidence class 1");
    assertTrue(stage9.inner_K_confidence_class[center] == 1, "warned inner node should have confidence class 1");
    assertTrue(stage9.outer_K_local_spike_count > 0, "outer local-spike count should increment");
    assertTrue(stage9.inner_K_local_spike_count > 0, "inner local-spike count should increment");

    const std::size_t corner = nodeIndex(0, 0, stage8.grid.nx);
    assertTrue(stage9.outer_K_confidence_class[corner] == 2, "nominal valid outer node should have confidence class 2");
    assertTrue(stage9.inner_K_confidence_class[corner] == 2, "nominal valid inner node should have confidence class 2");
}

void testStage9CsvExportIncludesQcColumns() {
    auto stage8 = makeSyntheticStage8ResultForStage9();
    const auto stage7 = makeSyntheticStage7ResultForStage9(stage8);
    FoldPatchAnalysisConfig config;
    config.output_prefix = "stage9_qc_columns";
    config.stage9_export_csv = true;

    const auto stage9 = runGeometryAnalysisStage9CurvatureComputation(stage7, stage8, config, nullptr);
    const auto outer_rows = readCsvRows(stage9.outer_curvature_csv_path);
    assertTrue(!outer_rows.empty(), "outer curvature csv should include a header");

    const auto& header = outer_rows[0];
    const auto hasColumn = [&header](const std::string& name) {
        return std::find(header.begin(), header.end(), name) != header.end();
    };
    assertTrue(hasColumn("H_raw"), "outer curvature csv should include H_raw");
    assertTrue(hasColumn("H_oriented"), "outer curvature csv should include H_oriented");
    assertTrue(hasColumn("K_raw"), "outer curvature csv should include K_raw");
    assertTrue(hasColumn("graph_normal_dot_radial"), "outer curvature csv should include graph_normal_dot_radial");
    assertTrue(hasColumn("orientation_flip_applied"), "outer curvature csv should include orientation_flip_applied");
    assertTrue(hasColumn("K_global_tail_flag"), "outer curvature csv should include K_global_tail_flag");
    assertTrue(hasColumn("K_local_spike_flag"), "outer curvature csv should include K_local_spike_flag");
    assertTrue(hasColumn("K_qc_warn_flag"), "outer curvature csv should include K_qc_warn_flag");
    assertTrue(hasColumn("K_confidence_class"), "outer curvature csv should include K_confidence_class");

    removeIfExists(stage9.outer_curvature_csv_path);
    removeIfExists(stage9.inner_curvature_csv_path);
    removeIfExists(stage9.curvature_valid_mask_csv_path);
    removeIfExists(stage9.curvature_summary_csv_path);
}

void testStage1ToStage6Integration() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.min_atoms_in_patch = 2;
    config.debug = true;
    config.stage6_export_obj_meshes = true;
    config.output_prefix = "stage6_integration";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1-6 integration should succeed");
    assertTrue(result.stage6_surfaces.success, "Stage 6 should succeed in full pipeline");
    assertTrue(result.stage6_surfaces.reconstructed_node_count > 0, "Stage 6 should reconstruct nodes");
    assertTrue(result.stage6_surfaces.final_valid_analysis_node_count > 0, "Stage 6 should produce final-valid nodes");
    assertTrue(std::isfinite(result.stage6_surfaces.mean_reconstructed_separation), "Stage 6 mean separation should be finite");
    assertTrue(std::filesystem::exists(result.stage6_surfaces.summary_csv_path), "Stage 6 summary csv should exist");
    assertTrue(std::filesystem::exists(result.stage6_surfaces.outer_obj_path), "Stage 6 outer obj should exist");
    assertTrue(std::filesystem::exists(result.stage6_surfaces.inner_obj_path), "Stage 6 inner obj should exist");

    std::filesystem::remove(result.stage2_patch.export_path);
    std::filesystem::remove(result.stage4_raw.outer_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_csv_path);
    std::filesystem::remove(result.stage4_raw.valid_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.outer_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.inner_only_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.negative_thickness_mask_csv_path);
    std::filesystem::remove(result.stage4_raw.summary_csv_path);
    std::filesystem::remove(result.stage4_raw.contact_atoms_pdb_path);
    removeIfExists(result.stage4_raw.stage3_normalized_atoms_csv_path);
    std::filesystem::remove(result.stage5_prep.outer_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.inner_seed_csv_path);
    std::filesystem::remove(result.stage5_prep.paired_seed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.boundary_exclusion_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.interp_allowed_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.hard_invalid_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.reliable_core_mask_csv_path);
    std::filesystem::remove(result.stage5_prep.summary_csv_path);
    std::filesystem::remove(result.stage6_surfaces.outer_reconstructed_csv_path);
    std::filesystem::remove(result.stage6_surfaces.inner_reconstructed_csv_path);
    std::filesystem::remove(result.stage6_surfaces.reconstructed_mask_csv_path);
    std::filesystem::remove(result.stage6_surfaces.final_valid_analysis_mask_csv_path);
    std::filesystem::remove(result.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    std::filesystem::remove(result.stage6_surfaces.summary_csv_path);
    std::filesystem::remove(result.stage6_surfaces.outer_obj_path);
    std::filesystem::remove(result.stage6_surfaces.inner_obj_path);
}

void testStage1ToStage7Integration() {
    Capsid capsid = makeSimpleCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.min_atoms_in_patch = 2;
    config.debug = true;
    config.stage6_export_obj_meshes = false;
    config.stage7_export_meshes = true;
    config.output_prefix = "stage7_integration";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1-7 integration should succeed");
    assertTrue(result.stage7_smooth.success, "Stage 7 should succeed in full pipeline");
    assertTrue(result.stage7_smooth.metric_domain_node_count > 0, "Stage 7 should produce non-empty metric domain");
    assertTrue(std::isfinite(result.stage7_smooth.mean_smooth_separation), "Stage 7 mean smooth separation should be finite");
    assertTrue(std::filesystem::exists(result.stage7_smooth.summary_csv_path), "Stage 7 summary csv should exist");
    assertTrue(std::filesystem::exists(result.stage7_smooth.outer_mesh_path), "Stage 7 outer mesh should exist");

    removeIfExists(result.stage2_patch.export_path);
    removeIfExists(result.stage4_raw.outer_csv_path);
    removeIfExists(result.stage4_raw.inner_csv_path);
    removeIfExists(result.stage4_raw.valid_mask_csv_path);
    removeIfExists(result.stage4_raw.outer_only_mask_csv_path);
    removeIfExists(result.stage4_raw.inner_only_mask_csv_path);
    removeIfExists(result.stage4_raw.negative_thickness_mask_csv_path);
    removeIfExists(result.stage4_raw.summary_csv_path);
    removeIfExists(result.stage4_raw.contact_atoms_pdb_path);
    removeIfExists(result.stage4_raw.stage3_normalized_atoms_csv_path);
    removeIfExists(result.stage5_prep.outer_seed_csv_path);
    removeIfExists(result.stage5_prep.inner_seed_csv_path);
    removeIfExists(result.stage5_prep.paired_seed_mask_csv_path);
    removeIfExists(result.stage5_prep.boundary_exclusion_mask_csv_path);
    removeIfExists(result.stage5_prep.interp_allowed_mask_csv_path);
    removeIfExists(result.stage5_prep.hard_invalid_mask_csv_path);
    removeIfExists(result.stage5_prep.reliable_core_mask_csv_path);
    removeIfExists(result.stage5_prep.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.inner_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.reconstructed_mask_csv_path);
    removeIfExists(result.stage6_surfaces.final_valid_analysis_mask_csv_path);
    removeIfExists(result.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage6_surfaces.summary_csv_path);
    removeIfExists(result.stage7_smooth.outer_smooth_csv_path);
    removeIfExists(result.stage7_smooth.inner_smooth_csv_path);
    removeIfExists(result.stage7_smooth.smooth_valid_mask_csv_path);
    removeIfExists(result.stage7_smooth.metric_domain_mask_csv_path);
    removeIfExists(result.stage7_smooth.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage7_smooth.summary_csv_path);
    removeIfExists(result.stage7_smooth.outer_mesh_path);
    removeIfExists(result.stage7_smooth.inner_mesh_path);
}

void testStage1ToStage9ThinPlateIntegration() {
    Capsid capsid = makeDenseCurvatureCapsid();
    FoldPatchAnalysisConfig config;
    config.enabled = true;
    config.fold_type = 2;
    config.fold_index = 0;
    config.cylinder_radius = 2.0;
    config.grid_spacing = 1.0;
    config.min_atoms_in_patch = 20;
    config.debug = true;
    config.stage6_export_obj_meshes = false;
    config.stage7_export_meshes = true;
    config.stage7_method = FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;
    config.stage8_enabled = true;
    config.stage8_export_csv = true;
    config.stage8_min_points = 6;
    config.stage8_fit_radius = 3.0;
    config.stage8_max_rms_residual = 1e6;
    config.stage8_max_abs_residual = 1e6;
    config.stage8_max_condition_indicator = 1e12;
    config.stage8_require_centered_support = false;
    config.stage8_min_directional_span = 0.1;
    config.stage9_enabled = true;
    config.stage9_export_csv = true;
    config.output_prefix = "stage7_fit_integration";

    const auto result = runFoldPatchGeometryAnalysis(capsid, config, makeParserConfig(), nullptr);
    assertTrue(result.success, "Stage 1-9 integration with thin-plate should succeed");
    assertTrue(result.stage7_smooth.success, "thin-plate Stage 7 should succeed in full pipeline");
    assertTrue(result.stage8_derivatives.success, "Stage 8 should succeed in full pipeline");
    assertTrue(result.stage7_smooth.stage7_method_label == "thin_plate_grid_fit", "Stage 7 method metadata should be populated");
    assertTrue(result.stage7_smooth.stage7_fit_active_node_count > 0, "thin-plate fit domain should be non-empty");
    assertTrue(result.stage8_derivatives.derivative_fit_attempted_node_count > 0,
               "Stage 8 should attempt derivative fits on the thin-plate metric domain");
    assertTrue(std::filesystem::exists(result.stage7_smooth.summary_csv_path), "thin-plate Stage 7 summary should exist");
    assertTrue(std::filesystem::exists(result.stage7_smooth.outer_mesh_path), "thin-plate Stage 7 mesh should exist");
    assertTrue(std::filesystem::exists(result.stage8_derivatives.outer_derivatives_csv_path),
               "Stage 8 outer derivative csv should exist");
    assertTrue(std::filesystem::exists(result.stage8_derivatives.derivative_summary_csv_path),
               "Stage 8 summary csv should exist");
    assertTrue(result.stage9_curvature.success, "Stage 9 should succeed in full pipeline");
    assertTrue(std::filesystem::exists(result.stage9_curvature.outer_curvature_csv_path),
               "Stage 9 outer curvature csv should exist");
    assertTrue(std::filesystem::exists(result.stage9_curvature.curvature_summary_csv_path),
               "Stage 9 summary csv should exist");

    removeIfExists(result.stage2_patch.export_path);
    removeIfExists(result.stage4_raw.outer_csv_path);
    removeIfExists(result.stage4_raw.inner_csv_path);
    removeIfExists(result.stage4_raw.valid_mask_csv_path);
    removeIfExists(result.stage4_raw.outer_only_mask_csv_path);
    removeIfExists(result.stage4_raw.inner_only_mask_csv_path);
    removeIfExists(result.stage4_raw.negative_thickness_mask_csv_path);
    removeIfExists(result.stage4_raw.summary_csv_path);
    removeIfExists(result.stage4_raw.contact_atoms_pdb_path);
    removeIfExists(result.stage4_raw.stage3_normalized_atoms_csv_path);
    removeIfExists(result.stage5_prep.outer_seed_csv_path);
    removeIfExists(result.stage5_prep.inner_seed_csv_path);
    removeIfExists(result.stage5_prep.paired_seed_mask_csv_path);
    removeIfExists(result.stage5_prep.boundary_exclusion_mask_csv_path);
    removeIfExists(result.stage5_prep.interp_allowed_mask_csv_path);
    removeIfExists(result.stage5_prep.hard_invalid_mask_csv_path);
    removeIfExists(result.stage5_prep.reliable_core_mask_csv_path);
    removeIfExists(result.stage5_prep.summary_csv_path);
    removeIfExists(result.stage6_surfaces.outer_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.inner_reconstructed_csv_path);
    removeIfExists(result.stage6_surfaces.reconstructed_mask_csv_path);
    removeIfExists(result.stage6_surfaces.final_valid_analysis_mask_csv_path);
    removeIfExists(result.stage6_surfaces.non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage6_surfaces.summary_csv_path);
    removeIfExists(result.stage7_smooth.outer_smooth_csv_path);
    removeIfExists(result.stage7_smooth.inner_smooth_csv_path);
    removeIfExists(result.stage7_smooth.smooth_valid_mask_csv_path);
    removeIfExists(result.stage7_smooth.metric_domain_mask_csv_path);
    removeIfExists(result.stage7_smooth.smooth_non_crossing_adjustment_mask_csv_path);
    removeIfExists(result.stage7_smooth.summary_csv_path);
    removeIfExists(result.stage7_smooth.outer_mesh_path);
    removeIfExists(result.stage7_smooth.inner_mesh_path);
    removeIfExists(result.stage8_derivatives.outer_derivatives_csv_path);
    removeIfExists(result.stage8_derivatives.inner_derivatives_csv_path);
    removeIfExists(result.stage8_derivatives.derivative_valid_mask_csv_path);
    removeIfExists(result.stage8_derivatives.derivative_failure_reason_csv_path);
    removeIfExists(result.stage8_derivatives.derivative_summary_csv_path);
    removeIfExists(result.stage9_curvature.outer_curvature_csv_path);
    removeIfExists(result.stage9_curvature.inner_curvature_csv_path);
    removeIfExists(result.stage9_curvature.curvature_valid_mask_csv_path);
    removeIfExists(result.stage9_curvature.curvature_summary_csv_path);
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
        testStage4FirstContactEnvelopeCompetingWinners();
        testStage4InnerWinnerCanDisagreeWithLowestCenterZ();
        testStage4FirstContactTieStability();
        testStage4PipelineUsesEnvelopeContactsAndProvenance();
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
        testStage4SummaryCsvIncludesExplicitPatchAndSerialProvenanceColumns();
        testStage5SeedProvenanceSeparatesPatchIndicesFromSerials();
        testStage5BoundaryExclusionBehavior();
        testStage5InterpolationAdmissibilityBehavior();
        testStage5ReliableCoreBehavior();
        testStage5DebugCsvExportAndDeterminism();
        testStage1ToStage5Integration();
        testStage6RequiresSuccessfulStage5();
        testStage6MinSeparationDefaultIsZero();
        testStage6SeedPreservationAndInterpolation();
        testStage6HardInvalidAndMasks();
        testStage6NonCrossingEnforcementAndConvergenceMetadata();
        testStage6DebugCsvArtifactsAndObjExportAndDeterminism();
        testStage6StlExport();
        testStage6CombinedExportDefault();
        testStage1ToStage6Integration();
        testStage7RequiresSuccessfulStage6();
        testStage7SeedPinSwitchAndRoughnessReduction();
        testStage7ReliableCoreDefinesMetricDomain();
        testStage7NonCrossingAndMetadata();
        testStage7MeshExportAndDeterminism();
        testStage7MethodDispatchAndMetadata();
        testStage7FitDomainWeightingMetadata();
        testStage7BoundaryModes();
        testStage7FitVsLegacyAndRoughness();
        testStage7ThinPlateDeterminism();
        testStage7S7S6DeltaCsvsThinPlate();
        testStage7S7S6DeltaCsvsLegacySmooth();
        testStage7S7S6DeltasDisabledByDefault();
        testStage8ExactQuadraticRecovery();
        testStage8InsufficientNeighborRejection();
        testStage8BoundarySkewRejection();
        testStage8RankDeficientOrPoorConditioningRejection();
        testStage8CsvExportSmoke();
        testStage9PlaneCurvature();
        testStage9ParaboloidCurvatureAtOrigin();
        testStage9InvalidInputPropagation();
        testStage9NonFiniteOutputGuard();
        testStage9CsvExportSmoke();
        testStage9QcWarnFlagsAndConfidenceClass();
        testStage9CsvExportIncludesQcColumns();
        testStage10SimpleValidThicknessField();
        testStage10OutsideDomainInvalidation();
        testStage10NegativeOrZeroInvalidation();
        testStage10ThresholdInvalidation();
        testStage10CurvatureValidRestrictionMode();
        testStage10SummaryCsvUsesCorrectedContractFields();
        testStage10Stage8FailureDoesNotBlockVerticalThickness();
        testStage10RadialFlatParallelPlanesDiffersFromVerticalOffAxis();
        testStage10RadialConcentricSphericalCapsApproxConstantGap();
        testStage10RadialNoBracketFailureIsTracked();
        testStage10RadialPOutInCsvExport();
        testStage10RadialInnerDomainLimitationMetadataAndFractions();
        testStage10RadialSamplingLossNoteEmissionGuard();
        testStage10RadialSummaryCsvSchemaUsesActiveAndLimitationFields();
        testStage1ToStage7Integration();
        testStage1ToStage9ThinPlateIntegration();
        std::cout << "All geometry analysis Stage 1/2/3/4/5/6/7/8/9/10 tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Geometry analysis test failure: " << e.what() << '\n';
        return 1;
    }
}
