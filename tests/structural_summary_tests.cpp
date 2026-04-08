#include "capsid.hpp"
#include "logger.hpp"
#include "pdb_parser.hpp"
#include "structural_summary.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

bool nearlyEqual(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

void assertNear(double actual, double expected, const std::string& label, double eps = 1e-9) {
    if (!nearlyEqual(actual, expected, eps)) {
        throw std::runtime_error(label + " mismatch: expected " + std::to_string(expected) +
                                 ", got " + std::to_string(actual));
    }
}

std::string makeRecord(const std::string& record_name,
                       int serial,
                       const std::string& atom_name,
                       const std::string& residue_name,
                       char chain_id,
                       int residue_seq,
                       double x,
                       double y,
                       double z,
                       const std::string& element) {
    char buffer[128];
    std::snprintf(buffer,
                  sizeof(buffer),
                  "%-6s%5d %-4s %-3s %1c%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s",
                  record_name.c_str(),
                  serial,
                  atom_name.c_str(),
                  residue_name.c_str(),
                  chain_id,
                  residue_seq,
                  x,
                  y,
                  z,
                  element.c_str());
    return std::string(buffer);
}

void testConstructedCapsidGeometry() {
    Capsid capsid("constructed");

    Chain chain_a(1, 'A');
    Residue res_a("ALA", 1, ' ', 'A', 1);
    res_a.addAtom(Atom(1, "CA", ' ', "ALA", 'A', 1, ' ', 0.0, 0.0, 0.0, 1.0, 20.0, "C", "", false));
    res_a.addAtom(Atom(2, "CB", ' ', "ALA", 'A', 1, ' ', 2.0, 0.0, 0.0, 1.0, 20.0, "C", "", false));
    chain_a.addResidue(res_a);

    Chain chain_b(2, 'B');
    Residue res_b("GLY", 1, ' ', 'B', 2);
    res_b.addAtom(Atom(3, "CA", ' ', "GLY", 'B', 1, ' ', 0.0, 2.0, 0.0, 1.0, 20.0, "C", "", false));
    res_b.addAtom(Atom(4, "N", ' ', "GLY", 'B', 1, ' ', 0.0, 0.0, 2.0, 1.0, 20.0, "N", "", false));
    chain_b.addResidue(res_b);

    capsid.addChain(chain_a);
    capsid.addChain(chain_b);
    capsid.finalizeCounts();

    const StructuralSummary summary = computeStructuralSummary(capsid);

    assertNear(static_cast<double>(summary.accepted_atom_count), 4.0, "accepted_atom_count");
    assertNear(summary.center_x, 0.5, "center_x");
    assertNear(summary.center_y, 0.5, "center_y");
    assertNear(summary.center_z, 0.5, "center_z");

    assertNear(summary.x_min, 0.0, "x_min");
    assertNear(summary.x_max, 2.0, "x_max");
    assertNear(summary.y_min, 0.0, "y_min");
    assertNear(summary.y_max, 2.0, "y_max");
    assertNear(summary.z_min, 0.0, "z_min");
    assertNear(summary.z_max, 2.0, "z_max");

    assertNear(summary.x_span, 2.0, "x_span");
    assertNear(summary.y_span, 2.0, "y_span");
    assertNear(summary.z_span, 2.0, "z_span");

    assertNear(summary.r_min, std::sqrt(0.75), "r_min", 1e-12);
    assertNear(summary.r_max, std::sqrt(2.75), "r_max", 1e-12);
    assertNear(summary.r_mean, (std::sqrt(0.75) + 3.0 * std::sqrt(2.75)) / 4.0, "r_mean", 1e-12);
    assertNear(summary.shell_thickness_estimate, std::sqrt(2.75) - std::sqrt(0.75), "shell_thickness", 1e-12);

    assertNear(summary.bbox_center_x, 1.0, "bbox_center_x");
    assertNear(summary.bbox_center_y, 1.0, "bbox_center_y");
    assertNear(summary.bbox_center_z, 1.0, "bbox_center_z");

    assertNear(summary.atoms_per_subunit.min, 2.0, "atoms_per_subunit.min");
    assertNear(summary.atoms_per_subunit.max, 2.0, "atoms_per_subunit.max");
    assertNear(summary.atoms_per_subunit.mean, 2.0, "atoms_per_subunit.mean");

    assertNear(summary.residues_per_subunit.min, 1.0, "residues_per_subunit.min");
    assertNear(summary.residues_per_subunit.max, 1.0, "residues_per_subunit.max");
    assertNear(summary.residues_per_subunit.mean, 1.0, "residues_per_subunit.mean");
}

void testParserIntegrationAcceptedAtomsOnly() {
    const std::string path = "tests_tmp_integration.pdb";
    {
        std::ofstream out(path);
        out << makeRecord("ATOM", 1, "N", "ALA", 'A', 1, 1.0, 0.0, 0.0, "N") << '\n';
        out << makeRecord("ATOM", 2, "CA", "ALA", 'A', 1, -1.0, 0.0, 0.0, "C") << '\n';
        out << makeRecord("HETATM", 3, "O", "HOH", 'A', 2, 100.0, 100.0, 100.0, "O") << '\n';
    }

    ParserConfig config;
    config.include_hetatm = false;
    config.protein_only = true;

    Logger logger;
    logger.setVerbosity(LogLevel::ERROR);

    PdbParser parser(config, &logger);
    Capsid capsid = parser.parseFile(path);

    const StructuralSummary summary = computeStructuralSummary(capsid);

    assertNear(static_cast<double>(summary.accepted_atom_count), 2.0, "integration.accepted_atom_count");
    assertNear(summary.center_x, 0.0, "integration.center_x");
    assertNear(summary.center_y, 0.0, "integration.center_y");
    assertNear(summary.center_z, 0.0, "integration.center_z");
    assertNear(summary.r_min, 1.0, "integration.r_min");
    assertNear(summary.r_max, 1.0, "integration.r_max");
    assertNear(summary.shell_thickness_estimate, 0.0, "integration.shell");

    std::remove(path.c_str());
}

void testPolicySensitiveSummaryFollowsAcceptedSet() {
    const std::string path = "tests_tmp_policy.pdb";
    {
        std::ofstream out(path);
        out << makeRecord("ATOM", 1, "N", "ALA", 'A', 1, 0.0, 0.0, 0.0, "N") << '\n';
        out << makeRecord("ATOM", 2, "CA", "ALA", 'A', 1, 2.0, 0.0, 0.0, "C") << '\n';
        out << makeRecord("HETATM", 3, "ZN", "ZN", 'A', 2, 10.0, 0.0, 0.0, "ZN") << '\n';
    }

    Logger logger;
    logger.setVerbosity(LogLevel::ERROR);

    ParserConfig protein_only_config;
    protein_only_config.include_hetatm = true;
    protein_only_config.protein_only = true;

    PdbParser parser_protein_only(protein_only_config, &logger);
    const StructuralSummary summary_protein_only =
        computeStructuralSummary(parser_protein_only.parseFile(path));

    assertNear(static_cast<double>(summary_protein_only.accepted_atom_count), 2.0,
               "policy.protein_only.atom_count");
    assertNear(summary_protein_only.center_x, 1.0, "policy.protein_only.center_x");

    ParserConfig include_non_protein_config;
    include_non_protein_config.include_hetatm = true;
    include_non_protein_config.protein_only = false;

    PdbParser parser_include_non_protein(include_non_protein_config, &logger);
    const StructuralSummary summary_include_non_protein =
        computeStructuralSummary(parser_include_non_protein.parseFile(path));

    assertNear(static_cast<double>(summary_include_non_protein.accepted_atom_count), 3.0,
               "policy.include_non_protein.atom_count");
    assertNear(summary_include_non_protein.center_x, 4.0, "policy.include_non_protein.center_x");

    if (summary_include_non_protein.shell_thickness_estimate <=
        summary_protein_only.shell_thickness_estimate) {
        throw std::runtime_error("policy shell-thickness estimate did not change with acceptance policy");
    }

    std::remove(path.c_str());
}

} // namespace

int main() {
    try {
        testConstructedCapsidGeometry();
        testParserIntegrationAcceptedAtomsOnly();
        testPolicySensitiveSummaryFollowsAcceptedSet();
        std::cout << "All structural summary tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Structural summary test failure: " << e.what() << '\n';
        return 1;
    }
}
