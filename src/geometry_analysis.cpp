#include "geometry_analysis.hpp"

#include "export_capsid.hpp"

#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
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

void validateStage4Config(const FoldPatchAnalysisConfig& config) {
    if (config.grid_spacing <= 0.0) {
        throw std::runtime_error("Stage 4 requires grid_spacing > 0");
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

std::size_t stage4NodeIndex(std::size_t i, std::size_t j, std::size_t nx) {
    return (j * nx) + i;
}

bool writeStage4CsvOuter(const GeometryStage4RawSheetResult& result) {
    std::ofstream out(result.outer_csv_path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk,valid,z_outer_raw\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(result.valid_mask[idx])
                << ',';
            if (result.valid_mask[idx] != 0) {
                out << result.z_outer_raw[idx];
            } else {
                out << "nan";
            }
            out << '\n';
        }
    }
    return out.good();
}

bool writeStage4CsvInner(const GeometryStage4RawSheetResult& result) {
    std::ofstream out(result.inner_csv_path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk,valid,z_inner_raw\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(result.valid_mask[idx])
                << ',';
            if (result.valid_mask[idx] != 0) {
                out << result.z_inner_raw[idx];
            } else {
                out << "nan";
            }
            out << '\n';
        }
    }
    return out.good();
}

bool writeStage4CsvValidMask(const GeometryStage4RawSheetResult& result) {
    std::ofstream out(result.valid_mask_csv_path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk,valid\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(result.valid_mask[idx])
                << '\n';
        }
    }
    return out.good();
}

PatchAtomContactRole classifyContactRole(bool used_as_outer, bool used_as_inner) {
    if (used_as_outer && used_as_inner) {
        return PatchAtomContactRole::both;
    }
    if (used_as_outer) {
        return PatchAtomContactRole::outer;
    }
    if (used_as_inner) {
        return PatchAtomContactRole::inner;
    }
    return PatchAtomContactRole::none;
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

LineSphereIntersection intersectVerticalLineWithSphere(double x,
                                                       double y,
                                                       const PatchAtom& atom,
                                                       double tolerance) {
    const double dx = x - atom.position.x;
    const double dy = y - atom.position.y;
    const double d2 = (dx * dx) + (dy * dy);
    const double r2 = atom.vdw_radius * atom.vdw_radius;

    if (d2 > r2 + tolerance) {
        return LineSphereIntersection{};
    }

    const double remaining = std::max(0.0, r2 - d2);
    const double delta = std::sqrt(remaining);
    LineSphereIntersection result;
    result.intersects = true;
    result.z_low = atom.position.z - delta;
    result.z_high = atom.position.z + delta;
    return result;
}

Stage4NodeFirstContact detectRawFirstContactAtNode(double x,
                                                    double y,
                                                    const std::vector<PatchAtom>& patch_atoms,
                                                    double tie_tolerance) {
    Stage4NodeFirstContact result;
    bool any_intersection = false;
    double best_outer = std::numeric_limits<double>::infinity();
    double best_inner = -std::numeric_limits<double>::infinity();

    for (std::size_t idx = 0; idx < patch_atoms.size(); ++idx) {
        const LineSphereIntersection intersection =
            intersectVerticalLineWithSphere(x, y, patch_atoms[idx], tie_tolerance);
        if (!intersection.intersects) {
            continue;
        }

        if (!any_intersection || (intersection.z_high < best_outer - tie_tolerance)) {
            best_outer = intersection.z_high;
            result.outer_patch_atom_index = idx;
        }

        if (!any_intersection || (intersection.z_low > best_inner + tie_tolerance)) {
            best_inner = intersection.z_low;
            result.inner_patch_atom_index = idx;
        }

        any_intersection = true;
    }

    if (!any_intersection) {
        return result;
    }

    result.z_outer_raw = best_outer;
    result.z_inner_raw = best_inner;
    result.valid = std::isfinite(result.z_outer_raw) && std::isfinite(result.z_inner_raw) &&
                   (result.z_inner_raw <= result.z_outer_raw + tie_tolerance);
    return result;
}

Stage4GridDescriptor buildStage4RegularGrid(double cylinder_radius, double spacing, double tolerance) {
    if (cylinder_radius <= 0.0) {
        throw std::runtime_error("Stage 4 requires cylinder_radius > 0");
    }
    if (spacing <= 0.0) {
        throw std::runtime_error("Stage 4 requires grid_spacing > 0");
    }

    Stage4GridDescriptor grid;
    grid.x_min = -cylinder_radius;
    grid.x_max = cylinder_radius;
    grid.y_min = -cylinder_radius;
    grid.y_max = cylinder_radius;
    grid.spacing = spacing;

    const std::size_t nx = static_cast<std::size_t>(std::floor((grid.x_max - grid.x_min) / spacing)) + 1;
    const std::size_t ny = static_cast<std::size_t>(std::floor((grid.y_max - grid.y_min) / spacing)) + 1;
    grid.nx = nx;
    grid.ny = ny;

    grid.x_values.reserve(nx);
    for (std::size_t i = 0; i < nx; ++i) {
        const double x = grid.x_min + (static_cast<double>(i) * spacing);
        const bool at_max = std::fabs(grid.x_max - x) <= tolerance;
        grid.x_values.push_back(at_max ? grid.x_max : x);
    }

    grid.y_values.reserve(ny);
    for (std::size_t j = 0; j < ny; ++j) {
        const double y = grid.y_min + (static_cast<double>(j) * spacing);
        const bool at_max = std::fabs(grid.y_max - y) <= tolerance;
        grid.y_values.push_back(at_max ? grid.y_max : y);
    }

    return grid;
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

GeometryStage4RawSheetResult runGeometryAnalysisStage4RawSheetDetection(
    const Capsid& capsid,
    const FoldPatchAnalysisConfig& config,
    const ParserConfig& parser_config,
    const GeometryPatchNormalizationResult& stage3_result,
    Logger* logger,
    double tolerance) {
    GeometryStage4RawSheetResult result;

    if (!stage3_result.success) {
        throw std::runtime_error("Stage 4 cannot run before successful Stage 3 patch normalization");
    }
    if (stage3_result.analytical_patch.atoms.empty()) {
        throw std::runtime_error("Analytical patch is empty; raw sheet detection cannot proceed");
    }
    validateStage4Config(config);

    result.messages.push_back("Geometry Stage 4");
    result.messages.push_back("Geometry analysis: starting Stage 4 raw outer/inner sheet detection.");
    result.messages.push_back("Geometry Stage 4 grid spacing: " + std::to_string(config.grid_spacing));
    result.messages.push_back("Geometry Stage 4 cylinder radius: " + std::to_string(config.cylinder_radius));

    result.grid = buildStage4RegularGrid(config.cylinder_radius, config.grid_spacing, tolerance);
    result.node_count = result.grid.nx * result.grid.ny;
    result.z_outer_raw.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.z_inner_raw.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.inside_disk_mask.assign(result.node_count, 0);
    result.valid_mask.assign(result.node_count, 0);
    result.atom_roles.assign(stage3_result.analytical_patch.atoms.size(), PatchAtomContactRole::none);

    const double radius2 = config.cylinder_radius * config.cylinder_radius;
    std::vector<bool> used_as_outer(stage3_result.analytical_patch.atoms.size(), false);
    std::vector<bool> used_as_inner(stage3_result.analytical_patch.atoms.size(), false);

    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            const double x = result.grid.x_values[i];
            const double y = result.grid.y_values[j];
            const double r2 = (x * x) + (y * y);
            const bool inside_disk = r2 <= radius2 + tolerance;

            if (!inside_disk) {
                continue;
            }

            result.inside_disk_mask[idx] = 1;
            ++result.inside_disk_count;

            const Stage4NodeFirstContact node =
                detectRawFirstContactAtNode(x, y, stage3_result.analytical_patch.atoms, tolerance);
            if (!node.valid) {
                continue;
            }

            result.valid_mask[idx] = 1;
            result.z_outer_raw[idx] = node.z_outer_raw;
            result.z_inner_raw[idx] = node.z_inner_raw;
            ++result.valid_node_count;

            used_as_outer[node.outer_patch_atom_index] = true;
            used_as_inner[node.inner_patch_atom_index] = true;

            Stage4RawContactRecord record;
            record.i = i;
            record.j = j;
            record.x = x;
            record.y = y;
            record.z_outer = node.z_outer_raw;
            record.z_inner = node.z_inner_raw;
            record.outer_patch_atom_index = node.outer_patch_atom_index;
            record.inner_patch_atom_index = node.inner_patch_atom_index;
            result.raw_contacts.push_back(record);
        }
    }

    if (result.inside_disk_count == 0) {
        throw std::runtime_error("No grid nodes fell inside the analysis disk");
    }

    result.invalid_node_count = result.inside_disk_count - result.valid_node_count;
    if (result.valid_node_count == 0) {
        throw std::runtime_error("Stage 4 produced zero valid raw-contact nodes");
    }

    std::vector<const Atom*> contact_atom_subset;
    contact_atom_subset.reserve(stage3_result.analytical_patch.atoms.size());
    for (std::size_t idx = 0; idx < stage3_result.analytical_patch.atoms.size(); ++idx) {
        result.atom_roles[idx] = classifyContactRole(used_as_outer[idx], used_as_inner[idx]);
        if (used_as_outer[idx]) {
            ++result.unique_outer_contact_atom_count;
        }
        if (used_as_inner[idx]) {
            ++result.unique_inner_contact_atom_count;
        }
        if (used_as_outer[idx] || used_as_inner[idx]) {
            contact_atom_subset.push_back(stage3_result.analytical_patch.atoms[idx].original_atom);
        }
    }

    result.outer_csv_path = config.output_prefix + "_outer_raw.csv";
    result.inner_csv_path = config.output_prefix + "_inner_raw.csv";
    result.valid_mask_csv_path = config.output_prefix + "_valid_mask_raw.csv";
    result.contact_atoms_pdb_path = config.output_prefix + "_raw_contact_points.pdb";

    if (!writeStage4CsvOuter(result)) {
        throw std::runtime_error("Failed to write raw outer grid CSV");
    }
    if (!writeStage4CsvInner(result)) {
        throw std::runtime_error("Failed to write raw inner grid CSV");
    }
    if (!writeStage4CsvValidMask(result)) {
        throw std::runtime_error("Failed to write raw valid-mask grid CSV");
    }

    ExportCapsidConfig writer_config;
    writer_config.output_path = result.contact_atoms_pdb_path;
    writer_config.emit_header_comments = true;
    writer_config.emit_ter_records = true;
    writer_config.emit_end_record = true;
    writer_config.preserve_atom_serial_numbers = true;
    writer_config.coordinate_records_only = false;
    writer_config.atom_subset = &contact_atom_subset;

    try {
        ExportCapsidWriter writer(logger);
        (void)writer.write(capsid, writer_config, parser_config);
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to export raw contact atoms subset");
    }

    result.messages.push_back("Geometry Stage 4 grid dimensions: nx=" + std::to_string(result.grid.nx) +
                              ", ny=" + std::to_string(result.grid.ny));
    result.messages.push_back("Geometry Stage 4 total nodes: " + std::to_string(result.node_count));
    result.messages.push_back("Geometry Stage 4 nodes inside disk: " + std::to_string(result.inside_disk_count));
    result.messages.push_back("Geometry Stage 4 valid nodes: " + std::to_string(result.valid_node_count));
    result.messages.push_back("Geometry Stage 4 invalid nodes: " + std::to_string(result.invalid_node_count));
    result.messages.push_back("Geometry Stage 4 unique outer-contact atoms: " +
                              std::to_string(result.unique_outer_contact_atom_count));
    result.messages.push_back("Geometry Stage 4 unique inner-contact atoms: " +
                              std::to_string(result.unique_inner_contact_atom_count));
    result.messages.push_back("Geometry Stage 4 outer CSV: " + result.outer_csv_path);
    result.messages.push_back("Geometry Stage 4 inner CSV: " + result.inner_csv_path);
    result.messages.push_back("Geometry Stage 4 valid mask CSV: " + result.valid_mask_csv_path);
    result.messages.push_back("Geometry Stage 4 contact atoms PDB: " + result.contact_atoms_pdb_path);
    result.messages.push_back("Geometry analysis: completed Stage 4 raw outer/inner sheet detection.");

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
    result.stage4_raw =
        runGeometryAnalysisStage4RawSheetDetection(capsid, config, parser_config, result.stage3_patch, logger);
    result.success = result.preparation.success && result.stage2_patch.success && result.stage3_patch.success &&
                     result.stage4_raw.success;
    logMessages(result.messages, logger);

    return result;
}
