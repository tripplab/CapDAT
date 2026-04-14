#ifndef CAPDAT_GEOMETRY_ANALYSIS_HPP
#define CAPDAT_GEOMETRY_ANALYSIS_HPP

#include <cstdint>
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
    /// Rotated working-frame position used by analytical geometry stages.
    geometry_symmetry::Vector3 position;
    /// Normalized element symbol (uppercase, trimmed, parser-safe).
    std::string element;
    /// Van der Waals radius (Å) resolved from normalized element.
    double vdw_radius = 0.0;
    /// Stage 2 cylinder-membership facts preserved for traceability/audits.
    CylinderMembership membership;
    double radial_xy = 0.0;
    bool in_positive_z = false;
    bool in_cylinder_radius = false;
    /// Stable pointer to the original selected atom in the Capsid.
    const Atom* original_atom = nullptr;
};

struct FoldPatchAnalysisConfig {
    enum class Stage7Method : uint8_t { smooth = 0, thin_plate_grid_fit = 1 };
    enum class Stage7BoundaryConditionMode : uint8_t { free = 0, fixed_to_stage6 = 1, soft_to_stage6 = 2 };
    enum class MeshExportFormat : uint8_t { obj = 0, stl = 1 };

    bool enabled = false;
    bool debug = false;
    int fold_type = 2;
    int fold_index = 0;
    double cylinder_radius = 12.0;
    double delta_vdw = 0.0;
    double grid_spacing = 2.0;
    std::size_t min_atoms_in_patch = 20;
    double stage5_boundary_margin = 0.0;
    double stage5_support_radius = 0.0;
    std::size_t stage5_min_support_nodes = 4;
    double stage5_reliable_radius = 0.0;
    double stage6_smoothing_weight = 1.0;
    std::size_t stage6_max_iterations = 500;
    double stage6_convergence_tolerance = 1e-6;
    bool stage6_enforce_non_crossing = true;
    double stage6_min_separation = 0.0;
    bool stage6_export_obj_meshes = true;
    bool stage7_enabled = true;
    Stage7Method stage7_method = Stage7Method::smooth;
    double stage7_smoothing_weight = 0.01;
    std::size_t stage7_max_iterations = 250;
    double stage7_convergence_tolerance = 1e-6;
    bool stage7_preserve_seed_values = false;
    double stage7_lambda = 1.0;
    double stage7_data_weight_seed = 1.0;
    double stage7_data_weight_interp = 0.25;
    bool stage7_use_reliable_core_only_for_fit = false;
    Stage7BoundaryConditionMode stage7_boundary_condition_mode = Stage7BoundaryConditionMode::soft_to_stage6;
    std::size_t stage7_solver_max_iterations = 500;
    double stage7_solver_tolerance = 1e-8;
    bool stage7_export_s7s6_deltas = false;
    bool stage7_enforce_non_crossing = true;
    double stage7_min_separation = 1.0;
    bool stage7_export_meshes = true;
    MeshExportFormat stage6_mesh_export_format = MeshExportFormat::stl;
    bool stage6_split_in_out_meshes = false;
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

struct AnalyticalPatch {
    std::vector<PatchAtom> atoms;
    std::vector<const Atom*> original_atom_refs;
    double cylinder_radius = 0.0;
    std::size_t atom_count = 0;
    std::size_t explicit_vdw_radius_count = 0;
    std::size_t inferred_vdw_radius_count = 0;
    std::size_t fallback_vdw_radius_count = 0;
    std::string export_path;
};

struct GeometryPatchNormalizationResult {
    bool success = false;
    AnalyticalPatch analytical_patch;
    std::vector<std::string> messages;
};

struct Stage4GridDescriptor {
    double x_min = 0.0;
    double x_max = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;
    double spacing = 0.0;
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::vector<double> x_values;
    std::vector<double> y_values;
};

enum class PatchAtomContactRole : uint8_t { none = 0, outer = 1, inner = 2, both = 3 };

struct Stage4RawContactRecord {
    std::size_t i = 0;
    std::size_t j = 0;
    double x = 0.0;
    double y = 0.0;
    double z_outer = 0.0;
    double z_inner = 0.0;
    std::size_t outer_patch_atom_index = 0;
    std::size_t inner_patch_atom_index = 0;
};

struct Stage4NodeFirstContact {
    bool valid = false;
    std::size_t candidate_patch_atom_count = 0;
    /// Winning outer raw envelope contact height (max ray-sphere z_high).
    double z_outer_raw = 0.0;
    /// Winning inner raw envelope contact height (min ray-sphere z_low).
    double z_inner_raw = 0.0;
    /// Patch-atom index of the winning outer z_high contact.
    std::size_t outer_patch_atom_index = 0;
    /// Patch-atom index of the winning inner z_low contact.
    std::size_t inner_patch_atom_index = 0;
};

struct LineSphereIntersection {
    bool intersects = false;
    double z_low = 0.0;
    double z_high = 0.0;
};

struct GeometryStage4RawSheetResult {
    bool success = false;
    Stage4GridDescriptor grid;
    std::vector<double> z_outer_raw;
    std::vector<double> z_inner_raw;
    std::vector<uint8_t> inside_disk_mask;
    std::vector<uint8_t> valid_mask;
    std::vector<std::size_t> candidate_patch_atom_counts;
    std::vector<int> inner_contact_serial_numbers;
    std::vector<int> outer_contact_serial_numbers;
    std::vector<int> outer_contact_patch_atom_indices;
    std::vector<int> inner_contact_patch_atom_indices;
    std::vector<PatchAtomContactRole> atom_roles;
    std::vector<Stage4RawContactRecord> raw_contacts;
    std::size_t contact_search_patch_atom_count = 0;
    std::size_t node_count = 0;
    std::size_t inside_disk_count = 0;
    std::size_t valid_node_count = 0;
    std::size_t invalid_node_count = 0;
    std::size_t outer_only_node_count = 0;
    std::size_t inner_only_node_count = 0;
    std::size_t both_hit_node_count = 0;
    std::size_t zero_thickness_node_count = 0;
    std::size_t negative_thickness_node_count = 0;
    std::size_t unique_outer_contact_atom_count = 0;
    std::size_t unique_inner_contact_atom_count = 0;
    std::size_t unique_both_contact_atom_count = 0;
    std::size_t unique_contact_atom_count = 0;
    std::vector<int> unique_outer_contact_serials;
    std::vector<int> unique_inner_contact_serials;
    std::vector<int> unique_both_contact_serials;
    std::vector<int> unique_contact_serial_union;
    std::size_t unique_outer_contact_serial_count = 0;
    std::size_t unique_inner_contact_serial_count = 0;
    std::size_t unique_both_contact_serial_count = 0;
    std::size_t unique_contact_serial_union_count = 0;
    std::string outer_csv_path;
    std::string inner_csv_path;
    std::string valid_mask_csv_path;
    std::string outer_only_mask_csv_path;
    std::string inner_only_mask_csv_path;
    std::string negative_thickness_mask_csv_path;
    std::string stage3_normalized_atoms_csv_path;
    std::string contact_atoms_pdb_path;
    std::string summary_csv_path;
    std::string stage4_start_timestamp_utc;
    std::string stage4_end_timestamp_utc;
    double stage4_runtime_seconds = 0.0;
    std::vector<std::string> messages;
};

struct GeometryStage5SurfacePrepResult {
    bool success = false;

    Stage4GridDescriptor grid;

    std::vector<double> z_outer_seed;
    std::vector<double> z_inner_seed;

    std::vector<uint8_t> inside_disk_mask;
    std::vector<uint8_t> raw_valid_mask;

    std::vector<uint8_t> outer_seed_mask;
    std::vector<uint8_t> inner_seed_mask;
    std::vector<uint8_t> paired_seed_mask;

    std::vector<uint8_t> boundary_exclusion_mask;
    std::vector<uint8_t> interp_allowed_outer_mask;
    std::vector<uint8_t> interp_allowed_inner_mask;
    std::vector<uint8_t> paired_interp_allowed_mask;
    std::vector<uint8_t> hard_invalid_mask;
    std::vector<uint8_t> reliable_core_mask;

    std::vector<int> outer_contact_serial_numbers;
    std::vector<int> inner_contact_serial_numbers;
    std::vector<int> outer_contact_patch_atom_indices;
    std::vector<int> inner_contact_patch_atom_indices;

    std::vector<int> unique_outer_seed_atom_serials;
    std::vector<int> unique_inner_seed_atom_serials;
    std::vector<std::size_t> unique_outer_seed_patch_atom_indices;
    std::vector<std::size_t> unique_inner_seed_patch_atom_indices;
    std::vector<std::size_t> unique_seed_patch_atom_index_union;

    std::size_t node_count = 0;
    std::size_t inside_disk_count = 0;
    std::size_t raw_valid_node_count = 0;
    std::size_t raw_invalid_node_count = 0;
    std::size_t outer_seed_node_count = 0;
    std::size_t inner_seed_node_count = 0;
    std::size_t paired_seed_node_count = 0;
    std::size_t boundary_excluded_node_count = 0;
    std::size_t interp_allowed_outer_node_count = 0;
    std::size_t interp_allowed_inner_node_count = 0;
    std::size_t paired_interp_allowed_node_count = 0;
    std::size_t hard_invalid_node_count = 0;
    std::size_t reliable_core_node_count = 0;
    std::size_t unique_outer_seed_patch_atom_count = 0;
    std::size_t unique_inner_seed_patch_atom_count = 0;
    std::size_t unique_seed_patch_atom_index_union_count = 0;

    double boundary_margin = 0.0;
    double support_radius = 0.0;
    double reliable_radius = 0.0;

    std::string outer_seed_csv_path;
    std::string inner_seed_csv_path;
    std::string paired_seed_mask_csv_path;
    std::string boundary_exclusion_mask_csv_path;
    std::string interp_allowed_mask_csv_path;
    std::string hard_invalid_mask_csv_path;
    std::string reliable_core_mask_csv_path;
    std::string summary_csv_path;

    std::vector<std::string> messages;
};

struct GeometryStage6SurfaceReconstructionResult {
    bool success = false;

    Stage4GridDescriptor grid;

    std::vector<double> z_outer_reconstructed;
    std::vector<double> z_inner_reconstructed;

    std::vector<uint8_t> inside_disk_mask;
    std::vector<uint8_t> paired_seed_mask;
    std::vector<uint8_t> paired_interp_allowed_mask;
    std::vector<uint8_t> hard_invalid_mask;
    std::vector<uint8_t> reconstructed_mask;
    std::vector<uint8_t> final_valid_analysis_mask;
    std::vector<uint8_t> non_crossing_adjustment_mask;
    std::vector<uint8_t> obj_vertex_mask;

    std::size_t node_count = 0;
    std::size_t reconstructed_node_count = 0;
    std::size_t seed_node_count = 0;
    std::size_t interp_node_count = 0;
    std::size_t final_valid_analysis_node_count = 0;
    std::size_t non_crossing_adjusted_node_count = 0;
    std::size_t unresolved_node_count = 0;

    std::size_t outer_iterations_used = 0;
    std::size_t inner_iterations_used = 0;
    double outer_final_max_update = 0.0;
    double inner_final_max_update = 0.0;
    double min_reconstructed_separation = 0.0;
    double max_reconstructed_separation = 0.0;
    double mean_reconstructed_separation = 0.0;

    std::size_t outer_obj_vertex_count = 0;
    std::size_t inner_obj_vertex_count = 0;
    std::size_t outer_obj_face_count = 0;
    std::size_t inner_obj_face_count = 0;
    std::size_t outer_reconstructed_scalar_node_count = 0;
    std::size_t inner_reconstructed_scalar_node_count = 0;
    std::size_t outer_reconstructed_nodes_not_used_by_any_face = 0;
    std::size_t inner_reconstructed_nodes_not_used_by_any_face = 0;

    std::string outer_reconstructed_csv_path;
    std::string inner_reconstructed_csv_path;
    std::string reconstructed_mask_csv_path;
    std::string final_valid_analysis_mask_csv_path;
    std::string non_crossing_adjustment_mask_csv_path;
    std::string outer_obj_path;
    std::string inner_obj_path;
    std::string summary_csv_path;

    std::vector<std::string> messages;
};

struct GeometryStage7SmoothedSurfaceResult {
    bool success = false;

    Stage4GridDescriptor grid;

    std::vector<double> z_outer_smooth;
    std::vector<double> z_inner_smooth;

    std::vector<uint8_t> reconstructed_mask;
    std::vector<uint8_t> reliable_core_mask;
    std::vector<uint8_t> smooth_valid_mask;
    std::vector<uint8_t> smooth_non_crossing_adjustment_mask;
    std::vector<uint8_t> metric_domain_mask;

    std::size_t node_count = 0;
    std::size_t smooth_valid_node_count = 0;
    std::size_t metric_domain_node_count = 0;
    std::size_t smooth_non_crossing_adjusted_node_count = 0;

    std::size_t outer_iterations_used = 0;
    std::size_t inner_iterations_used = 0;
    double outer_final_max_update = 0.0;
    double inner_final_max_update = 0.0;
    std::string stage7_method_label = "smooth";
    double stage7_lambda = 0.0;
    double stage7_data_weight_seed = 0.0;
    double stage7_data_weight_interp = 0.0;
    std::string stage7_data_weight_policy;
    bool stage7_use_reliable_core_only_for_fit = false;
    std::string stage7_boundary_condition_mode_label;
    std::size_t stage7_fit_active_node_count = 0;
    std::size_t stage7_fit_seed_like_node_count = 0;
    std::size_t stage7_fit_interp_like_node_count = 0;
    double outer_fit_final_residual = 0.0;
    double inner_fit_final_residual = 0.0;
    double outer_fit_max_abs_residual = 0.0;
    double inner_fit_max_abs_residual = 0.0;
    double outer_fit_bending_energy = 0.0;
    double inner_fit_bending_energy = 0.0;
    std::size_t outer_solver_iterations_used = 0;
    std::size_t inner_solver_iterations_used = 0;
    double outer_solver_final_update = 0.0;
    double inner_solver_final_update = 0.0;

    double min_smooth_separation = 0.0;
    double max_smooth_separation = 0.0;
    double mean_smooth_separation = 0.0;
    double max_abs_outer_delta = 0.0;
    double max_abs_inner_delta = 0.0;
    double max_abs_thickness_delta = 0.0;
    double mean_abs_outer_delta = 0.0;
    double mean_abs_inner_delta = 0.0;
    double mean_abs_thickness_delta = 0.0;

    std::string normal_orientation_outer = "toward_positive_z";
    std::string normal_orientation_inner = "toward_negative_z";
    std::string local_thickness_definition = "vertical_z_difference";
    std::string metric_surface_definition = "stage6_reconstructed_then_stage7_smoothed_restricted_to_reliable_core";

    std::string outer_smooth_csv_path;
    std::string inner_smooth_csv_path;
    std::string smooth_valid_mask_csv_path;
    std::string metric_domain_mask_csv_path;
    std::string smooth_non_crossing_adjustment_mask_csv_path;
    std::string outer_delta_csv_path;
    std::string inner_delta_csv_path;
    std::string thickness_delta_csv_path;
    std::string summary_csv_path;

    std::string outer_mesh_path;
    std::string inner_mesh_path;

    std::size_t outer_mesh_vertex_count = 0;
    std::size_t inner_mesh_vertex_count = 0;
    std::size_t outer_mesh_face_count = 0;
    std::size_t inner_mesh_face_count = 0;
    std::size_t outer_mesh_unused_scalar_nodes = 0;
    std::size_t inner_mesh_unused_scalar_nodes = 0;

    std::vector<std::string> messages;
};

struct GeometryAnalysisResult {
    bool success = false;
    GeometryPreparationResult preparation;
    GeometryPatchSelectionResult stage2_patch;
    GeometryPatchNormalizationResult stage3_patch;
    GeometryStage4RawSheetResult stage4_raw;
    GeometryStage5SurfacePrepResult stage5_prep;
    GeometryStage6SurfaceReconstructionResult stage6_surfaces;
    GeometryStage7SmoothedSurfaceResult stage7_smooth;
    std::vector<std::string> messages;
};

CylinderMembership classifyPatchCylinder(const geometry_symmetry::Vector3& position,
                                         double cylinder_radius);
std::string normalizeElementSymbol(const std::string& raw_element);
double vdwRadius(const std::string& normalized_element);
PatchAtom makePatchAtom(const Atom& atom,
                        const geometry_symmetry::Vector3& rotated_position,
                        const CylinderMembership& membership,
                        double delta_vdw);

LineSphereIntersection intersectVerticalLineWithSphere(double x,
                                                       double y,
                                                       const PatchAtom& atom,
                                                       double tolerance = 1e-12);

Stage4NodeFirstContact detectRawFirstContactAtNode(double x,
                                                    double y,
                                                    const std::vector<PatchAtom>& patch_atoms,
                                                    double tie_tolerance = 1e-12);

Stage4GridDescriptor buildStage4RegularGrid(double cylinder_radius,
                                            double spacing,
                                            double tolerance = 1e-12);

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

GeometryPatchNormalizationResult runGeometryAnalysisStage3PatchNormalization(
    const GeometryPatchSelectionResult& stage2_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger);

GeometryStage4RawSheetResult runGeometryAnalysisStage4RawSheetDetection(
    const Capsid& capsid,
    const FoldPatchAnalysisConfig& config,
    const ParserConfig& parser_config,
    const GeometryPatchNormalizationResult& stage3_result,
    Logger* logger,
    double tolerance = 1e-12);

GeometryStage5SurfacePrepResult runGeometryAnalysisStage5SurfacePreparation(
    const GeometryStage4RawSheetResult& stage4_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance = 1e-12);

GeometryStage6SurfaceReconstructionResult runGeometryAnalysisStage6SurfaceReconstruction(
    const GeometryStage5SurfacePrepResult& stage5_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance = 1e-12);

GeometryStage7SmoothedSurfaceResult runGeometryAnalysisStage7SurfaceSmoothing(
    const GeometryStage6SurfaceReconstructionResult& stage6_result,
    const GeometryStage5SurfacePrepResult& stage5_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance = 1e-12);

GeometryAnalysisResult runFoldPatchGeometryAnalysis(Capsid& capsid,
                                                    const FoldPatchAnalysisConfig& config,
                                                    const ParserConfig& parser_config,
                                                    Logger* logger,
                                                    double tolerance = 1e-9);

#endif // CAPDAT_GEOMETRY_ANALYSIS_HPP
