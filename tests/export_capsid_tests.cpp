#include "export_capsid.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

void assertTrue(bool ok, const std::string& msg) {
    if (!ok) throw std::runtime_error(msg);
}

Capsid makeCapsid() {
    Capsid capsid("in-memory");
    Chain chain(1, 'A');
    Residue residue("ALA", 1, ' ', 'A', 1);
    residue.addAtom(Atom(1, "CA", ' ', "ALA", 'A', 1, ' ', 1.0, 2.0, 3.0, 1.0, 20.0, "C", "", false));
    chain.addResidue(residue);
    capsid.addChain(chain);
    capsid.finalizeCounts();
    return capsid;
}

std::string readFile(const std::string& path) {
    std::ifstream in(path);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

void testUntransformedRemark() {
    Capsid capsid = makeCapsid();
    ParserConfig parser_config;
    parser_config.protein_only = true;

    ExportCapsidConfig cfg;
    cfg.output_path = "tmp_export_untransformed.pdb";

    ExportCapsidWriter writer(nullptr);
    (void)writer.write(capsid, cfg, parser_config);

    const std::string text = readFile(cfg.output_path);
    assertTrue(text.find("original parsed input frame") != std::string::npos,
               "untransformed export should describe original frame");
    assertTrue(text.find("reoriented in place") == std::string::npos,
               "untransformed export should not claim reorientation");
    std::remove(cfg.output_path.c_str());
}

void testTransformedRemarkAndCoordinates() {
    Capsid capsid = makeCapsid();

    // Simulate workflow-applied in-place transform and authoritative state.
    capsid.mutableChains()[0].mutableResidues()[0].mutableAtoms()[0].setPosition(9.0, 8.0, 7.0);
    Capsid::OrientationState state;
    state.in_original_parsed_frame = false;
    state.reoriented_in_place = true;
    state.source_mode = Capsid::OrientationSourceMode::custom_vector;
    state.source_description = "1,0,0";
    state.requested_target_axis = 'y';
    state.target_direction = {0.0, 1.0, 0.0};
    state.has_rotation_angle = true;
    state.rotation_angle_radians = 1.57079632679;
    capsid.setOrientationState(state);

    ParserConfig parser_config;
    ExportCapsidConfig cfg;
    cfg.output_path = "tmp_export_transformed.pdb";

    ExportCapsidWriter writer(nullptr);
    (void)writer.write(capsid, cfg, parser_config);

    const std::string text = readFile(cfg.output_path);
    assertTrue(text.find("reoriented in place") != std::string::npos,
               "transformed export should declare in-place reorientation");
    assertTrue(text.find("Target axis: y") != std::string::npos,
               "transformed export should report selected target axis");
    assertTrue(text.find("  9.000   8.000   7.000") != std::string::npos,
               "exported coordinates should match current in-memory coordinates");

    std::remove(cfg.output_path.c_str());
}

} // namespace

int main() {
    try {
        testUntransformedRemark();
        testTransformedRemarkAndCoordinates();
        std::cout << "All export capsid tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Export capsid test failure: " << e.what() << '\n';
        return 1;
    }
}
