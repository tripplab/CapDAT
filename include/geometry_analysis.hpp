#ifndef CAPDAT_GEOMETRY_ANALYSIS_HPP
#define CAPDAT_GEOMETRY_ANALYSIS_HPP

#include <string>
#include <vector>

#include "capsid.hpp"
#include "geometry_symmetry.hpp"
#include "logger.hpp"
#include "pdb_parser.hpp"
#include "reorientation_workflow.hpp"

struct FoldPatchAnalysisConfig {
    bool enabled = false;
    int fold_type = 2;
    int fold_index = 0;
    double cylinder_radius = 12.0;
    bool export_rotated_capsid = false;
    std::string output_prefix = "geometry";
};

struct GeometryPreparationResult {
    bool success = false;
    int requested_fold_type = 0;
    int requested_fold_index = 0;
    std::string resolved_fold_name;
    geometry_symmetry::Vector3 resolved_fold_reference_vector;
    geometry_symmetry::Vector3 resolved_fold_unit_vector;
    bool used_identity_rotation = false;
    geometry_symmetry::Matrix3 rotation_matrix;
    geometry_symmetry::Vector3 rotation_axis;
    double rotation_angle_radians = 0.0;
    bool coordinates_modified_in_place = false;
    std::string final_frame_description;
    std::string export_path;
    ReorientationResult reorientation_result;
    std::vector<std::string> messages;
};

struct GeometryAnalysisResult {
    bool success = false;
    GeometryPreparationResult preparation;
    std::vector<std::string> messages;
};

GeometryPreparationResult prepareGeometryAnalysisStage1(Capsid& capsid,
                                                        const FoldPatchAnalysisConfig& config,
                                                        const ParserConfig& parser_config,
                                                        Logger* logger,
                                                        double tolerance = 1e-9);

GeometryAnalysisResult runFoldPatchGeometryAnalysis(Capsid& capsid,
                                                    const FoldPatchAnalysisConfig& config,
                                                    const ParserConfig& parser_config,
                                                    Logger* logger,
                                                    double tolerance = 1e-9);

#endif // CAPDAT_GEOMETRY_ANALYSIS_HPP
