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
    bool enabled = false;
    int fold_type = 2;
    int fold_index = 0;
    double cylinder_radius = 12.0;
    double grid_spacing = 2.0;
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
    double z_outer_raw = 0.0;
    double z_inner_raw = 0.0;
    std::size_t outer_patch_atom_index = 0;
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

struct GeometryAnalysisResult {
    bool success = false;
    GeometryPreparationResult preparation;
    GeometryPatchSelectionResult stage2_patch;
    GeometryPatchNormalizationResult stage3_patch;
    GeometryStage4RawSheetResult stage4_raw;
    std::vector<std::string> messages;
};

CylinderMembership classifyPatchCylinder(const geometry_symmetry::Vector3& position,
                                         double cylinder_radius);
std::string normalizeElementSymbol(const std::string& raw_element);
double vdwRadius(const std::string& normalized_element);
PatchAtom makePatchAtom(const Atom& atom,
                        const geometry_symmetry::Vector3& rotated_position,
                        const CylinderMembership& membership);

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
    Logger* logger);

GeometryStage4RawSheetResult runGeometryAnalysisStage4RawSheetDetection(
    const Capsid& capsid,
    const FoldPatchAnalysisConfig& config,
    const ParserConfig& parser_config,
    const GeometryPatchNormalizationResult& stage3_result,
    Logger* logger,
    double tolerance = 1e-12);

GeometryAnalysisResult runFoldPatchGeometryAnalysis(Capsid& capsid,
                                                    const FoldPatchAnalysisConfig& config,
                                                    const ParserConfig& parser_config,
                                                    Logger* logger,
                                                    double tolerance = 1e-9);

#endif // CAPDAT_GEOMETRY_ANALYSIS_HPP
