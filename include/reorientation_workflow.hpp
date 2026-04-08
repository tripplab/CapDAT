#ifndef CAPDAT_REORIENTATION_WORKFLOW_HPP
#define CAPDAT_REORIENTATION_WORKFLOW_HPP

#include <array>
#include <string>
#include <vector>

#include "capsid.hpp"
#include "geometry_symmetry.hpp"
#include "logger.hpp"

/**
 * @brief Source selector for user-requested reorientation.
 */
enum class ReorientationSourceMode {
    fold,
    custom_vector
};

/**
 * @brief Request model consumed by the reorientation workflow layer.
 *
 * Main.cpp should only populate this from CLI inputs; all source/axis
 * resolution and validation beyond high-level wiring belongs in the workflow.
 */
struct ReorientationRequest {
    bool request_reorientation = false;
    ReorientationSourceMode source_mode = ReorientationSourceMode::fold;
    std::string fold_name;
    std::string custom_vector_text;
    char target_axis = 'z';
    bool request_export = false;
    std::string export_path;
    bool verbose = false;
};

enum class ReorientationStatus {
    success,
    identity,
    degenerate_input,
    invalid_request
};

/**
 * @brief Structured workflow result for testing and runtime messaging.
 */
struct ReorientationResult {
    ReorientationStatus status = ReorientationStatus::invalid_request;
    std::string resolved_source_description;
    geometry_symmetry::Vector3 source_direction;
    char requested_target_axis = 'z';
    geometry_symmetry::Vector3 target_direction;
    geometry_symmetry::Matrix3 rotation_matrix;
    geometry_symmetry::Vector3 rotation_axis;
    double rotation_angle_radians = 0.0;
    bool coordinates_modified_in_place = false;
    bool export_requested = false;
    std::string export_path;
    std::vector<std::string> messages;
};

/**
 * @brief Apply optional in-place reorientation to the current Capsid.
 */
ReorientationResult applyReorientationWorkflow(Capsid& capsid,
                                               const ReorientationRequest& request,
                                               Logger* logger,
                                               double tolerance = 1e-9);

/**
 * @brief Parse "x,y,z" CLI vector text into a Cartesian vector.
 */
geometry_symmetry::Vector3 parseCliVector(const std::string& vector_text);

#endif // CAPDAT_REORIENTATION_WORKFLOW_HPP
