#ifndef CAPDAT_GEOMETRY_ANALYSIS_HPP
#define CAPDAT_GEOMETRY_ANALYSIS_HPP

#include <string>
#include <vector>

#include "atom.hpp"
#include "capsid.hpp"
#include "geometry_symmetry.hpp"
#include "logger.hpp"
#include "pdb_parser.hpp"
#include "reorientation_workflow.hpp"

struct CylinderMembership {
    double radial_xy = 0.0;
    bool in_positive_z = false;
    bool in_cylinder_radius = false;
    bool selected = false;
};

struct PatchAtom {
    geometry_symmetry::Vector3 position;
    std::string element;
    double radial_xy = 0.0;
    bool in_positive_z = false;
    bool in_cylinder_radius = false;
    const Atom* original_atom = nullptr;
};

struct FoldPatchAnalysisConfig {
    bool enabled = false;
    int fold_type = 2;
    int fold_index = 0;
    double cylinder_radius = 12.0;
    std::size_t min_atoms_in_patch = 20;
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

struct GeometryPatchSelectionResult {
    bool success = false;
    std::size_t total_atoms_examined = 0;
    std::size_t selected_atoms_count = 0;
    std::size_t rejected_non_positive_z_count = 0;
    std::size_t rejected_outside_radius_count = 0;
    double cylinder_radius = 0.0;
    std::size_t min_atoms_in_patch = 0;
    std::string export_path;
    std::vector<PatchAtom> patch_atoms;
    std::vector<const Atom*> selected_atom_refs;
    std::vector<std::string> messages;
};

struct GeometryAnalysisResult {
    bool success = false;
    GeometryPreparationResult preparation;
    GeometryPatchSelectionResult stage2_patch;
    std::vector<std::string> messages;
};

CylinderMembership classifyPatchCylinder(const geometry_symmetry::Vector3& position,
                                         double cylinder_radius);

GeometryPreparationResult prepareGeometryAnalysisStage1(Capsid& capsid,
                                                        const FoldPatchAnalysisConfig& config,
                                                        const ParserConfig& parser_config,
                                                        Logger* logger,
                                                        double tolerance = 1e-9);

GeometryPatchSelectionResult runGeometryAnalysisStage2PatchSelection(
    const Capsid& capsid,
    const FoldPatchAnalysisConfig& config,
    const ParserConfig& parser_config,
    const GeometryPreparationResult& stage1_result,
    Logger* logger);

GeometryAnalysisResult runFoldPatchGeometryAnalysis(Capsid& capsid,
                                                    const FoldPatchAnalysisConfig& config,
                                                    const ParserConfig& parser_config,
                                                    Logger* logger,
                                                    double tolerance = 1e-9);

#endif // CAPDAT_GEOMETRY_ANALYSIS_HPP
