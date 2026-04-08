#include "reorientation_workflow.hpp"

#include <sstream>
#include <stdexcept>

namespace {

geometry_symmetry::Vector3 targetAxisToDirection(char axis) {
    // Axis-selection convention for v01 reorientation:
    // - x maps to +X, y maps to +Y, z maps to +Z.
    // - only lowercase axis tokens are accepted by workflow validation.
    switch (axis) {
        case 'x':
            return {1.0, 0.0, 0.0};
        case 'y':
            return {0.0, 1.0, 0.0};
        case 'z':
            return {0.0, 0.0, 1.0};
        default:
            throw std::runtime_error("Invalid target axis. Allowed values: x, y, z.");
    }
}

std::array<std::array<double, 3>, 3> toCapsidMatrix(const geometry_symmetry::Matrix3& matrix) {
    return matrix.m;
}

std::array<double, 3> toCapsidVector(const geometry_symmetry::Vector3& v) {
    return {v.x, v.y, v.z};
}

Capsid::OrientationSourceMode toCapsidSourceMode(ReorientationSourceMode mode) {
    return mode == ReorientationSourceMode::fold
               ? Capsid::OrientationSourceMode::fold
               : Capsid::OrientationSourceMode::custom_vector;
}

} // namespace

geometry_symmetry::Vector3 parseCliVector(const std::string& vector_text) {
    std::stringstream ss(vector_text);
    std::string token;
    std::array<double, 3> values = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < 3; ++i) {
        if (!std::getline(ss, token, ',')) {
            throw std::runtime_error("Invalid --align-vector format. Expected x,y,z.");
        }
        try {
            values[i] = std::stod(token);
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid --align-vector component: " + token);
        }
    }
    if (std::getline(ss, token, ',')) {
        throw std::runtime_error("Invalid --align-vector format. Expected exactly three components.");
    }
    return {values[0], values[1], values[2]};
}

ReorientationResult applyReorientationWorkflow(Capsid& capsid,
                                               const ReorientationRequest& request,
                                               Logger* logger,
                                               double tolerance) {
    ReorientationResult result;
    result.export_requested = request.request_export;
    result.export_path = request.export_path;

    if (!request.request_reorientation) {
        result.status = ReorientationStatus::identity;
        result.messages.push_back("Reorientation not requested; capsid remains in original parsed frame.");
        return result;
    }

    geometry_symmetry::Vector3 source_direction;
    std::string source_description;

    if (request.source_mode == ReorientationSourceMode::fold) {
        const auto& fold = geometry_symmetry::foldByName(request.fold_name);
        source_direction = fold.unit_vector;
        source_description = fold.name;
        result.messages.push_back("Resolved canonical fold source: " + fold.name);
    } else {
        const geometry_symmetry::Vector3 parsed = parseCliVector(request.custom_vector_text);
        source_direction = geometry_symmetry::normalize(parsed);
        source_description = request.custom_vector_text;
        result.messages.push_back("Resolved custom direction source: " + request.custom_vector_text);
    }

    const geometry_symmetry::Vector3 target_direction = targetAxisToDirection(request.target_axis);
    result.messages.push_back(std::string("Resolved target axis: ") + request.target_axis);

    // Reorientation convention: this workflow requests a pure rotation about the
    // origin that maps source direction to a canonical positive target axis.
    const geometry_symmetry::RotationDefinition rotation =
        geometry_symmetry::alignDirectionToDirection(source_direction, target_direction, tolerance);

    bool modified_coordinates = false;

    if (rotation.status != geometry_symmetry::RotationStatus::identity) {
        for (Chain& chain : capsid.mutableChains()) {
            for (Residue& residue : chain.mutableResidues()) {
                for (Atom& atom : residue.mutableAtoms()) {
                    const geometry_symmetry::Vector3 rotated =
                        geometry_symmetry::rotatePoint(rotation.matrix, {atom.x(), atom.y(), atom.z()});
                    atom.setPosition(rotated.x, rotated.y, rotated.z);
                }
            }
        }
        modified_coordinates = true;
    }

    Capsid::OrientationState state;
    state.in_original_parsed_frame = false;
    state.reoriented_in_place = true;
    state.already_aligned_identity = rotation.status == geometry_symmetry::RotationStatus::identity;
    state.applied_rotation_matrix = toCapsidMatrix(rotation.matrix);
    state.has_rotation_axis = true;
    state.rotation_axis = toCapsidVector(rotation.axis);
    state.has_rotation_angle = true;
    state.rotation_angle_radians = rotation.angle_radians;
    state.source_mode = toCapsidSourceMode(request.source_mode);
    state.source_description = source_description;
    state.source_direction = toCapsidVector(source_direction);
    state.requested_target_axis = request.target_axis;
    state.target_direction = toCapsidVector(target_direction);
    capsid.setOrientationState(state);

    result.resolved_source_description = source_description;
    result.source_direction = source_direction;
    result.requested_target_axis = request.target_axis;
    result.target_direction = target_direction;
    result.rotation_matrix = rotation.matrix;
    result.rotation_axis = rotation.axis;
    result.rotation_angle_radians = rotation.angle_radians;
    result.coordinates_modified_in_place = modified_coordinates;

    if (rotation.status == geometry_symmetry::RotationStatus::identity) {
        result.status = ReorientationStatus::identity;
        result.messages.push_back("Requested alignment already satisfied; identity transform recorded.");
    } else {
        result.status = ReorientationStatus::success;
        result.messages.push_back("Applied in-place coordinate reorientation to current Capsid.");
    }

    if (logger != nullptr) {
        for (const std::string& message : result.messages) {
            logger->info(message);
        }
    }

    return result;
}
