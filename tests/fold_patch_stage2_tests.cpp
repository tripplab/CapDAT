#include "fold_patch_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <unistd.h>

#include "geometry_symmetry.hpp"
#include "logger.hpp"
#include "pdb_parser.hpp"

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

Capsid buildSyntheticCapsid(const std::vector<std::tuple<double, double, double, std::string, bool>>& atoms) {
    Capsid capsid("synthetic");
    Chain chain(0, 'A');
    Residue residue("GLY", 1, ' ', 'A', 0);

    int serial = 1;
    for (const auto& [x, y, z, element, is_hetatm] : atoms) {
        residue.addAtom(Atom(serial,
                             "CA",
                             ' ',
                             "GLY",
                             'A',
                             1,
                             ' ',
                             x,
                             y,
                             z,
                             1.0,
                             0.0,
                             element,
                             "",
                             is_hetatm));
        ++serial;
    }

    chain.addResidue(std::move(residue));
    capsid.addChain(std::move(chain));
    capsid.finalizeCounts();
    return capsid;
}

Capsid parseRealCapsid() {
    Logger logger;
    logger.setVerbosity(LogLevel::ERROR);
    PdbParser parser(ParserConfig{}, &logger);
    return parser.parseFile(std::string(CAPDAT_SOURCE_DIR) + "/data/1cwp_full.vdb");
}

std::string captureStderrForUnknownElements() {
    const std::string tmp_path = "vdw_unknown_stderr.txt";
    const int saved_stderr = dup(fileno(stderr));
    if (saved_stderr < 0) {
        throw std::runtime_error("Failed to duplicate stderr fd");
    }

    FILE* capture = std::freopen(tmp_path.c_str(), "w", stderr);
    if (capture == nullptr) {
        close(saved_stderr);
        throw std::runtime_error("Failed to redirect stderr for test capture");
    }

    const double r1 = vdwRadius("ZZ");
    const double r2 = vdwRadius("");
    const double r3 = vdwRadius("FE");

    std::fflush(stderr);
    dup2(saved_stderr, fileno(stderr));
    close(saved_stderr);

    assertNear(r1, 1.70, "fallback ZZ", 0.0);
    assertNear(r2, 1.70, "fallback empty", 0.0);
    assertNear(r3, 1.70, "fallback FE", 0.0);

    std::ifstream input(tmp_path);
    std::stringstream buffer;
    buffer << input.rdbuf();
    input.close();

    return buffer.str();
}

void testVdwRadiusKnownElements() {
    assertNear(vdwRadius("H"), 1.20, "vdW H", 0.0);
    assertNear(vdwRadius("c"), 1.70, "vdW C", 0.0);
    assertNear(vdwRadius(" N "), 1.55, "vdW N", 0.0);
    assertNear(vdwRadius("o"), 1.52, "vdW O", 0.0);
    assertNear(vdwRadius("S"), 1.80, "vdW S", 0.0);
    assertNear(vdwRadius("p"), 1.80, "vdW P", 0.0);
}

void testVdwRadiusUnknownElementFallback() {
    const std::string captured = captureStderrForUnknownElements();
    assertTrue(!captured.empty(), "Expected stderr warning output for unknown elements");
}

void testPatchAtomConstruction() {
    const Capsid capsid = buildSyntheticCapsid({{5.0, 0.0, 3.0, "c", false}});

    LocalFrame frame;
    frame.origin = Eigen::Vector3d::Zero();
    frame.axes = Eigen::Matrix3d::Identity();
    frame.fold_name = "2_0";

    FoldPatchAnalysisConfig config;
    config.R = 12.0;
    config.z_margin = 11.0;

    const std::vector<PatchAtom> patch = extractPatch(capsid, frame, config);
    assertTrue(patch.size() == 1, "Expected one selected atom in synthetic patch");

    const PatchAtom& pa = patch.front();
    assertNear(pa.x, 5.0, "patch x", 1e-6);
    assertNear(pa.y, 0.0, "patch y", 1e-6);
    assertNear(pa.z, 3.0, "patch z", 1e-6);
    assertNear(pa.vdw_radius, 1.70, "patch vdW", 1e-6);
}

void testCylinderMembershipPredicate() {
    const Capsid capsid = buildSyntheticCapsid({
        {3.0, 4.0, 2.0, "C", false},
        {9.0, 9.0, 1.0, "C", false},
        {2.0, 2.0, -1.0, "C", false},
        {0.0, 0.0, 12.0, "C", false},
        {0.0, 11.9, 1.0, "C", false},
    });

    LocalFrame frame;
    frame.origin = Eigen::Vector3d::Zero();
    frame.axes = Eigen::Matrix3d::Identity();
    frame.fold_name = "2_0";

    FoldPatchAnalysisConfig config;
    config.R = 12.0;
    config.z_margin = 11.0;

    const std::vector<PatchAtom> patch = extractPatch(capsid, frame, config);
    assertTrue(patch.size() == 3, "Expected exactly three atoms in cylinder");
    assertTrue(patch[0].serial == 1, "First selected atom should be serial 1");
    assertTrue(patch[1].serial == 4, "Second selected atom should be serial 4");
    assertTrue(patch[2].serial == 5, "Third selected atom should be serial 5");
}

void testLocalFrameOrthonormality() {
    const Capsid capsid = parseRealCapsid();
    FoldPatchAnalysisConfig config;
    config.fold_name = "2_0";

    const LocalFrame frame = buildFoldFrame(capsid, config);
    assertTrue(!frame.fold_name.empty(), "buildFoldFrame failed for 2_0");

    double ortho_accum = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double dot_col = 0.0;
            for (int k = 0; k < 3; ++k) {
                dot_col += frame.axes.m[k][i] * frame.axes.m[k][j];
            }
            const double expected = (i == j) ? 1.0 : 0.0;
            const double diff = dot_col - expected;
            ortho_accum += diff * diff;
        }
    }
    assertTrue(std::sqrt(ortho_accum) < 1e-9, "Local frame must be orthonormal");

    const geometry_symmetry::Vector3 expected_axis = geometry_symmetry::foldUnitVector("2_0");
    const double dot_axis = frame.axes.m[0][2] * expected_axis.x +
                            frame.axes.m[1][2] * expected_axis.y +
                            frame.axes.m[2][2] * expected_axis.z;
    const double alignment = 1.0 - std::fabs(dot_axis);
    assertTrue(alignment < 1e-6, "Frame Z axis must align to fold axis");
}

void testIntegrationCoherentExtraction() {
    const Capsid capsid = parseRealCapsid();
    FoldPatchAnalysisConfig config;
    config.min_atoms_in_patch = 1;

    const LocalFrame frame = buildFoldFrame(capsid, config);
    const std::vector<PatchAtom> patch = extractPatch(capsid, frame, config);

    assertTrue(!patch.empty(), "Patch atom extraction should yield at least one atom");

    const auto [zmin_it, zmax_it] = std::minmax_element(
        patch.begin(), patch.end(), [](const PatchAtom& a, const PatchAtom& b) {
            return a.z < b.z;
        });

    assertTrue(zmin_it->z >= 0.0, "Patch z_min must be >= 0");
    assertTrue(zmax_it->z > 0.0, "Patch z_max must be > 0");

    for (const PatchAtom& atom : patch) {
        assertTrue(atom.vdw_radius > 0.0, "All vdW radii must be positive");
    }
}

void testIntegrationPatchPdbExport() {
    const Capsid capsid = parseRealCapsid();
    FoldPatchAnalysisConfig config;
    config.export_patch_pdb = true;
    config.export_grid_pdb = false;
    config.export_surface_obj = false;
    config.export_surface_ply = false;
    config.export_surface_stl = false;
    config.export_fields_csv = false;
    config.export_summary_json = false;
    config.export_summary_csv = false;
    config.min_atoms_in_patch = 1;
    config.out_prefix = "./stage2_patch_test";

    const FoldPatchAnalysisResult result = runFoldPatchAnalysis(capsid, config);
    assertTrue(result.success, "runFoldPatchAnalysis should succeed for real capsid");

    const std::string patch_path = config.out_prefix + "_patch.pdb";
    std::ifstream input(patch_path);
    assertTrue(input.good(), "Patch PDB output file should exist");

    std::string line;
    std::size_t atom_lines = 0;
    std::string last_non_empty;
    while (std::getline(input, line)) {
        if (line.rfind("ATOM", 0) == 0) {
            ++atom_lines;
            assertTrue(line.size() == 80, "ATOM line width must be exactly 80 chars");
        }
        if (!line.empty()) {
            last_non_empty = line;
        }
    }

    assertTrue(atom_lines > 0, "Patch PDB must contain at least one ATOM record");
    assertTrue(last_non_empty == "END", "Patch PDB must terminate with END");
}

} // namespace

int main() {
    try {
        testVdwRadiusKnownElements();
        testVdwRadiusUnknownElementFallback();
        testPatchAtomConstruction();
        testCylinderMembershipPredicate();
        testLocalFrameOrthonormality();
        testIntegrationCoherentExtraction();
        testIntegrationPatchPdbExport();
        std::cout << "All fold patch stage-2 tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Fold patch stage-2 test failure: " << e.what() << '\n';
        return 1;
    }
}
