#include "geometry_analysis.hpp"

#include "export_capsid.hpp"

#include <cmath>
#include <cctype>
#include <stdexcept>
#include <unordered_map>

namespace {

constexpr double kVdwRadiusFallbackAngstrom = 1.70;

void logMessages(const std::vector<std::string>& messages, Logger* logger) {
    if (logger == nullptr) {
        return;
    }
    for (const std::string& message : messages) {
        logger->info(message);
    }
}

std::string describeFrame(const Capsid::OrientationState& state) {
    if (state.in_original_parsed_frame) {
        return "original_parsed_frame";
    }
    if (state.reoriented_in_place) {
        return "derived_reoriented_frame";
    }
    return "unknown_frame";
}

void validateStage1Config(const FoldPatchAnalysisConfig& config) {
    if (config.fold_type != 2 && config.fold_type != 3 && config.fold_type != 5) {
        throw std::runtime_error("Invalid geometry fold type: expected one of 2, 3, 5");
    }
    if (config.cylinder_radius <= 0.0) {
        throw std::runtime_error("Geometry cylinder radius must be > 0");
    }
}

void validateStage2Config(const FoldPatchAnalysisConfig& config) {
    if (config.cylinder_radius <= 0.0) {
        throw std::runtime_error("Geometry patch selection requires cylinder_radius > 0");
    }
    if (config.min_atoms_in_patch == 0) {
        throw std::runtime_error("Geometry patch selection requires min_atoms_in_patch > 0");
    }
}

std::string trimWhitespace(const std::string& text) {
    std::size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start])) != 0) {
        ++start;
    }

    std::size_t end = text.size();
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
        --end;
    }

    return text.substr(start, end - start);
}

const std::unordered_map<std::string, double>& vdwRadiusLookupTable() {
    static const std::unordered_map<std::string, double> table = {
        {"H", 1.20},
        {"C", 1.70},
        {"N", 1.55},
        {"O", 1.52},
        {"S", 1.80},
        {"P", 1.80},
        {"SE", 1.90},
    };
    return table;
}

} // namespace

CylinderMembership classifyPatchCylinder(const geometry_symmetry::Vector3& position, double cylinder_radius) {
    const double radial_xy = std::sqrt((position.x * position.x) + (position.y * position.y));
    const bool in_positive_z = position.z > 0.0;
    const bool in_cylinder_radius = radial_xy <= cylinder_radius;

    CylinderMembership membership;
    membership.radial_xy = radial_xy;
    membership.in_positive_z = in_positive_z;
    membership.in_cylinder_radius = in_cylinder_radius;
    membership.selected = in_positive_z && in_cylinder_radius;
    return membership;
}

std::string normalizeElementSymbol(const std::string& raw_element) {
    const std::string trimmed = trimWhitespace(raw_element);
    if (trimmed.empty()) {
        return "";
    }

    std::string normalized;
    normalized.reserve(trimmed.size());
    for (const char ch : trimmed) {
        if (std::isalpha(static_cast<unsigned char>(ch)) == 0) {
            continue;
        }
        normalized.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
        if (normalized.size() == 2) {
            break;
        }
    }
    return normalized;
}

double vdwRadius(const std::string& normalized_element) {
    const auto& table = vdwRadiusLookupTable();
    const auto it = table.find(normalized_element);
    if (it != table.end()) {
        return it->second;
    }
    // Stage 3 analytical fallback: unknown/unsupported elements map to 1.70 Å.
    return kVdwRadiusFallbackAngstrom;
}

PatchAtom makePatchAtom(const Atom& atom,
                        const geometry_symmetry::Vector3& rotated_position,
                        const CylinderMembership& membership) {
    PatchAtom patch_atom;
    patch_atom.position = rotated_position;
    patch_atom.element = normalizeElementSymbol(atom.element());
    patch_atom.vdw_radius = vdwRadius(patch_atom.element);
    patch_atom.membership = membership;
    patch_atom.radial_xy = membership.radial_xy;
    patch_atom.in_positive_z = membership.in_positive_z;
    patch_atom.in_cylinder_radius = membership.in_cylinder_radius;
    patch_atom.original_atom = &atom;
    return patch_atom;
}

GeometryPreparationResult prepareGeometryAnalysisStage1(Capsid& capsid,
                                                        const FoldPatchAnalysisConfig& config,
                                                        const ParserConfig& parser_config,
                                                        Logger* logger,
                                                        double tolerance) {
    GeometryPreparationResult result;
    result.requested_fold_type = config.fold_type;
    result.requested_fold_index = config.fold_index;

    validateStage1Config(config);

    const geometry_symmetry::FoldDefinition& fold =
        geometry_symmetry::foldByTypeIndex(config.fold_type, config.fold_index);

    result.resolved_fold_name = fold.name;
    result.resolved_fold_reference_vector = fold.reference_vector;
    result.resolved_fold_unit_vector = fold.unit_vector;

    result.messages.push_back("Geometry Stage 1");
    result.messages.push_back("Geometry Stage 1 request: fold_type=" + std::to_string(config.fold_type) +
                              ", fold_index=" + std::to_string(config.fold_index));
    result.messages.push_back("Geometry Stage 1 resolved canonical fold: " + fold.name);
    result.messages.push_back("Geometry Stage 1 target axis: +Z");

    ReorientationRequest request;
    request.request_reorientation = true;
    request.source_mode = ReorientationSourceMode::fold;
    request.fold_name = fold.name;
    request.target_axis = 'z';
    request.verbose = config.export_rotated_capsid;

    result.reorientation_result = applyReorientationWorkflow(capsid, request, logger, tolerance);
    result.used_identity_rotation = result.reorientation_result.status == ReorientationStatus::identity;
    result.rotation_matrix = result.reorientation_result.rotation_matrix;
    result.rotation_axis = result.reorientation_result.rotation_axis;
    result.rotation_angle_radians = result.reorientation_result.rotation_angle_radians;
    result.coordinates_modified_in_place = result.reorientation_result.coordinates_modified_in_place;

    const auto& orientation = capsid.orientationState();
    if (!orientation.reoriented_in_place || orientation.requested_target_axis != 'z') {
        throw std::runtime_error("Geometry preparation failed during reorientation to +Z");
    }

    result.final_frame_description = describeFrame(orientation);
    result.messages.push_back(std::string("Geometry Stage 1 transform is identity: ") +
                              (result.used_identity_rotation ? "true" : "false"));
    result.messages.push_back(std::string("Geometry Stage 1 coordinates modified in place: ") +
                              (result.coordinates_modified_in_place ? "true" : "false"));
    result.messages.push_back("Geometry Stage 1 final working frame: " + result.final_frame_description);

    if (config.export_rotated_capsid) {
        ExportCapsidConfig writer_config;
        writer_config.output_path = config.output_prefix + "_rotated_capsid.pdb";
        writer_config.emit_header_comments = true;
        writer_config.emit_ter_records = true;
        writer_config.emit_end_record = true;
        writer_config.preserve_atom_serial_numbers = true;
        writer_config.coordinate_records_only = false;

        ExportCapsidWriter writer(logger);
        (void)writer.write(capsid, writer_config, parser_config);
        result.export_path = writer_config.output_path;
        result.messages.push_back("Geometry Stage 1 exported rotated capsid: " + result.export_path);
    }

    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryPatchSelectionResult runGeometryAnalysisStage2PatchSelection(
    const Capsid& capsid,
    const FoldPatchAnalysisConfig& config,
    const ParserConfig& parser_config,
    const GeometryPreparationResult& stage1_result,
    Logger* logger) {
    GeometryPatchSelectionResult result;
    result.cylinder_radius = config.cylinder_radius;
    result.min_atoms_in_patch = config.min_atoms_in_patch;

    if (!stage1_result.success) {
        throw std::runtime_error("Stage 2 cannot run before successful Stage 1 geometry preparation");
    }
    const auto& orientation = capsid.orientationState();
    if (!orientation.reoriented_in_place || orientation.requested_target_axis != 'z') {
        throw std::runtime_error("Stage 2 cannot run before successful Stage 1 geometry preparation");
    }

    validateStage2Config(config);

    result.messages.push_back("Geometry Stage 2");
    result.messages.push_back("Geometry analysis: starting Stage 2 cylindrical patch selection.");
    result.messages.push_back("Geometry Stage 2 cylinder radius: " + std::to_string(config.cylinder_radius));

    for (const Chain& chain : capsid.chains()) {
        for (const Residue& residue : chain.residues()) {
            for (const Atom& atom : residue.atoms()) {
                ++result.total_atoms_examined;
                const geometry_symmetry::Vector3 position{atom.x(), atom.y(), atom.z()};
                const CylinderMembership membership = classifyPatchCylinder(position, config.cylinder_radius);

                if (!membership.in_positive_z) {
                    ++result.rejected_non_positive_z_count;
                }
                if (!membership.in_cylinder_radius) {
                    ++result.rejected_outside_radius_count;
                }

                if (!membership.selected) {
                    continue;
                }

                result.patch_atoms.push_back(makePatchAtom(atom, position, membership));
                result.selected_atom_refs.push_back(&atom);
            }
        }
    }

    result.selected_atoms_count = result.patch_atoms.size();
    if (result.selected_atoms_count < config.min_atoms_in_patch) {
        throw std::runtime_error("Patch selection produced " + std::to_string(result.selected_atoms_count) +
                                 " atoms, below min_atoms_in_patch=" + std::to_string(config.min_atoms_in_patch));
    }

    ExportCapsidConfig writer_config;
    writer_config.output_path = config.output_prefix + "_patch_atoms.pdb";
    writer_config.emit_header_comments = true;
    writer_config.emit_ter_records = true;
    writer_config.emit_end_record = true;
    writer_config.preserve_atom_serial_numbers = true;
    writer_config.coordinate_records_only = false;
    writer_config.atom_subset = &result.selected_atom_refs;

    try {
        ExportCapsidWriter writer(logger);
        (void)writer.write(capsid, writer_config, parser_config);
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to export patch atom subset to PDB");
    }

    result.export_path = writer_config.output_path;
    result.messages.push_back("Geometry Stage 2 examined atoms: " + std::to_string(result.total_atoms_examined));
    result.messages.push_back("Geometry Stage 2 selected atoms: " + std::to_string(result.selected_atoms_count));
    result.messages.push_back("Geometry Stage 2 rejected z<=0: " + std::to_string(result.rejected_non_positive_z_count));
    result.messages.push_back("Geometry Stage 2 rejected radial cutoff: " +
                              std::to_string(result.rejected_outside_radius_count));
    result.messages.push_back("Geometry Stage 2 min_atoms_in_patch: " + std::to_string(result.min_atoms_in_patch));
    result.messages.push_back("Geometry Stage 2 exported patch atoms: " + result.export_path);
    result.messages.push_back("Geometry analysis: completed Stage 2 cylindrical patch selection.");
    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryPatchNormalizationResult runGeometryAnalysisStage3PatchNormalization(
    const GeometryPatchSelectionResult& stage2_result,
    Logger* logger) {
    GeometryPatchNormalizationResult result;
    result.analytical_patch.cylinder_radius = stage2_result.cylinder_radius;
    result.analytical_patch.export_path = stage2_result.export_path;

    if (!stage2_result.success) {
        throw std::runtime_error("Stage 3 cannot run before successful Stage 2 patch selection");
    }
    if (stage2_result.patch_atoms.size() != stage2_result.selected_atom_refs.size()) {
        throw std::runtime_error("Analytical patch build produced inconsistent atom/reference counts");
    }

    result.messages.push_back("Geometry Stage 3");
    result.messages.push_back("Geometry analysis: starting Stage 3 analytical patch normalization.");
    result.messages.push_back("Geometry Stage 3 selected atoms to normalize: " +
                              std::to_string(stage2_result.patch_atoms.size()));

    result.analytical_patch.atoms.reserve(stage2_result.patch_atoms.size());
    result.analytical_patch.original_atom_refs.reserve(stage2_result.selected_atom_refs.size());

    for (std::size_t idx = 0; idx < stage2_result.patch_atoms.size(); ++idx) {
        const PatchAtom& selected = stage2_result.patch_atoms[idx];
        const Atom* original_ref = stage2_result.selected_atom_refs[idx];
        if (original_ref == nullptr || selected.original_atom == nullptr || selected.original_atom != original_ref) {
            throw std::runtime_error("PatchAtom normalization encountered a null original atom reference");
        }

        const PatchAtom normalized = makePatchAtom(*original_ref, selected.position, selected.membership);
        if (vdwRadiusLookupTable().find(normalized.element) == vdwRadiusLookupTable().end()) {
            ++result.analytical_patch.fallback_vdw_radius_count;
        } else {
            ++result.analytical_patch.explicit_vdw_radius_count;
        }

        result.analytical_patch.atoms.push_back(normalized);
        result.analytical_patch.original_atom_refs.push_back(original_ref);
    }

    result.analytical_patch.atom_count = result.analytical_patch.atoms.size();
    if (result.analytical_patch.atom_count != result.analytical_patch.original_atom_refs.size()) {
        throw std::runtime_error("Analytical patch build produced inconsistent atom/reference counts");
    }
    if (result.analytical_patch.atom_count == 0) {
        throw std::runtime_error("Analytical patch is empty after Stage 3 normalization");
    }

    result.messages.push_back("Geometry Stage 3 explicit element vdW assignments: " +
                              std::to_string(result.analytical_patch.explicit_vdw_radius_count));
    result.messages.push_back("Geometry Stage 3 fallback vdW assignments: " +
                              std::to_string(result.analytical_patch.fallback_vdw_radius_count));
    result.messages.push_back("Geometry Stage 3 patch export path (Stage 2 canonical): " +
                              result.analytical_patch.export_path);
    result.messages.push_back("Geometry analysis: completed Stage 3 analytical patch normalization.");
    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryAnalysisResult runFoldPatchGeometryAnalysis(Capsid& capsid,
                                                    const FoldPatchAnalysisConfig& config,
                                                    const ParserConfig& parser_config,
                                                    Logger* logger,
                                                    double tolerance) {
    GeometryAnalysisResult result;
    if (!config.enabled) {
        result.success = true;
        result.messages.push_back("Geometry analysis disabled; skipping Stage 1 preparation.");
        return result;
    }

    result.preparation = prepareGeometryAnalysisStage1(capsid, config, parser_config, logger, tolerance);
    result.stage2_patch =
        runGeometryAnalysisStage2PatchSelection(capsid, config, parser_config, result.preparation, logger);
    result.stage3_patch = runGeometryAnalysisStage3PatchNormalization(result.stage2_patch, logger);
    result.success = result.preparation.success && result.stage2_patch.success && result.stage3_patch.success;
    logMessages(result.messages, logger);

    return result;
}
