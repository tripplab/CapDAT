#include "capsid.hpp"
#include "geometry_symmetry.hpp"
#include "reorientation_workflow.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

bool near(double a, double b, double eps = 1e-8) { return std::fabs(a - b) <= eps; }

void assertTrue(bool ok, const std::string& msg) {
    if (!ok) throw std::runtime_error(msg);
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

void testFoldToXInPlace() {
    Capsid capsid = makeSimpleCapsid();
    const double x_before = capsid.chains()[0].residues()[0].atoms()[0].x();

    ReorientationRequest request;
    request.request_reorientation = true;
    request.source_mode = ReorientationSourceMode::fold;
    request.fold_name = "5_0";
    request.target_axis = 'x';

    const ReorientationResult result = applyReorientationWorkflow(capsid, request, nullptr);
    assertTrue(result.status == ReorientationStatus::success || result.status == ReorientationStatus::identity,
               "fold request should resolve successfully");
    assertTrue(result.requested_target_axis == 'x', "axis x should be preserved");

    const auto& state = capsid.orientationState();
    assertTrue(state.reoriented_in_place, "capsid should be marked reoriented");
    assertTrue(!state.in_original_parsed_frame, "capsid should no longer be in original frame");
    assertTrue(state.requested_target_axis == 'x', "stored requested axis should be x");

    const double x_after = capsid.chains()[0].residues()[0].atoms()[0].x();
    if (result.status == ReorientationStatus::success) {
        assertTrue(!near(x_before, x_after), "coordinates should change for non-identity request");
    }

    const geometry_symmetry::Vector3 rotated = geometry_symmetry::rotateDirection(
        result.rotation_matrix, geometry_symmetry::foldUnitVector("5_0"));
    assertTrue(near(rotated.x, 1.0, 1e-6) && near(rotated.y, 0.0, 1e-6) && near(rotated.z, 0.0, 1e-6),
               "fold should align to +X");
}

void testCustomVectorDefaultAxisZ() {
    Capsid capsid = makeSimpleCapsid();
    ReorientationRequest request;
    request.request_reorientation = true;
    request.source_mode = ReorientationSourceMode::custom_vector;
    request.custom_vector_text = "1,0,0";
    request.target_axis = 'z';

    const ReorientationResult result = applyReorientationWorkflow(capsid, request, nullptr);
    assertTrue(result.requested_target_axis == 'z', "default/selected z axis should be retained");
    assertTrue(capsid.orientationState().requested_target_axis == 'z', "capsid state should store z axis");
}

void testZeroVectorRejected() {
    Capsid capsid = makeSimpleCapsid();
    bool threw = false;

    try {
        ReorientationRequest request;
        request.request_reorientation = true;
        request.source_mode = ReorientationSourceMode::custom_vector;
        request.custom_vector_text = "0,0,0";
        request.target_axis = 'y';
        (void)applyReorientationWorkflow(capsid, request, nullptr);
    } catch (const std::runtime_error&) {
        threw = true;
    }

    assertTrue(threw, "zero vector should be rejected");
    assertTrue(capsid.orientationState().in_original_parsed_frame,
               "failed request must not partially modify orientation state");
}

} // namespace

int main() {
    try {
        testFoldToXInPlace();
        testCustomVectorDefaultAxisZ();
        testZeroVectorRejected();
        std::cout << "All reorientation workflow tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Reorientation workflow test failure: " << e.what() << '\n';
        return 1;
    }
}
