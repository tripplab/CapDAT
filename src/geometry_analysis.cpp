#include "geometry_analysis.hpp"

#include "export_capsid.hpp"
#include "timer.hpp"

#include <cmath>
#include <filesystem>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace {

constexpr double kVdwRadiusFallbackAngstrom = 1.70;
constexpr std::string_view kWarningMessagePrefix = "[[WARNING]] ";
constexpr std::string_view kNoteMessagePrefix = "[[NOTE]] ";
constexpr std::string_view kInfoMessagePrefix = "[[INFO]] ";


std::string geometryArtifactPath(const FoldPatchAnalysisConfig& config, const std::string& filename) {
    const std::filesystem::path root_dir(config.output_root_dir);
    std::error_code ec;
    std::filesystem::create_directories(root_dir, ec);
    return (root_dir / filename).string();
}

void logMessages(const std::vector<std::string>& messages, Logger* logger) {
    if (logger == nullptr) {
        return;
    }
    for (const std::string& message : messages) {
        if (message.rfind(kWarningMessagePrefix, 0) == 0) {
            logger->warning(message.substr(kWarningMessagePrefix.size()));
            continue;
        }
        if (message.rfind(kNoteMessagePrefix, 0) == 0) {
            logger->note(message.substr(kNoteMessagePrefix.size()));
            continue;
        }
        if (message.rfind(kInfoMessagePrefix, 0) == 0) {
            logger->info(message.substr(kInfoMessagePrefix.size()));
            continue;
        }
        logger->debug(message);
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

std::string jsonEscape(const std::string& input) {
    std::string escaped;
    escaped.reserve(input.size());
    for (const char ch : input) {
        switch (ch) {
        case '\\':
            escaped += "\\\\";
            break;
        case '"':
            escaped += "\\\"";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            escaped.push_back(ch);
            break;
        }
    }
    return escaped;
}

std::string meshFormatLabel(FoldPatchAnalysisConfig::MeshExportFormat format) {
    return format == FoldPatchAnalysisConfig::MeshExportFormat::stl ? "stl" : "obj";
}

std::string formatScientific(double value) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(6) << value;
    return out.str();
}

const char* stage7MethodLabel(FoldPatchAnalysisConfig::Stage7Method method) {
    switch (method) {
    case FoldPatchAnalysisConfig::Stage7Method::smooth:
        return "smooth";
    case FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit:
        return "thin_plate_grid_fit";
    }
    return "unknown";
}

const char* stage7BoundaryConditionModeLabel(FoldPatchAnalysisConfig::Stage7BoundaryConditionMode mode) {
    switch (mode) {
    case FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::free:
        return "free";
    case FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::fixed_to_stage6:
        return "fixed_to_stage6";
    case FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::soft_to_stage6:
        return "soft_to_stage6";
    }
    return "unknown";
}

const char* stage10ThicknessMethodLabel(FoldPatchAnalysisConfig::Stage10ThicknessMethod method) {
    switch (method) {
    case FoldPatchAnalysisConfig::Stage10ThicknessMethod::vertical:
        return "vertical";
    case FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial:
        return "radial";
    }
    return "unknown";
}

bool writeGeometryRunSummaryJson(const FoldPatchAnalysisConfig& config,
                                 const ParserConfig& parser_config,
                                 const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "{\n";
    out << "  \"geometry\": {\n";
    out << "    \"enabled\": " << (config.enabled ? "true" : "false") << ",\n";
    out << "    \"debug\": " << (config.debug ? "true" : "false") << ",\n";
    out << "    \"fold_type\": " << config.fold_type << ",\n";
    out << "    \"fold_index\": " << config.fold_index << ",\n";
    out << "    \"cylinder_radius\": " << config.cylinder_radius << ",\n";
    out << "    \"delta_vdw\": " << config.delta_vdw << ",\n";
    out << "    \"grid_spacing\": " << config.grid_spacing << ",\n";
    out << "    \"min_atoms_in_patch\": " << config.min_atoms_in_patch << ",\n";
    out << "    \"output_root_dir\": \"" << jsonEscape(config.output_root_dir) << "\",\n";
    out << "    \"fold_name\": \"" << jsonEscape(config.fold_name) << "\"\n";
    out << "  },\n";
    out << "  \"stage5\": {\n";
    out << "    \"boundary_margin\": " << config.stage5_boundary_margin << ",\n";
    out << "    \"support_radius\": " << config.stage5_support_radius << ",\n";
    out << "    \"min_support_nodes\": " << config.stage5_min_support_nodes << ",\n";
    out << "    \"reliable_radius\": " << config.stage5_reliable_radius << "\n";
    out << "  },\n";
    out << "  \"stage6\": {\n";
    out << "    \"smoothing_weight\": " << config.stage6_smoothing_weight << ",\n";
    out << "    \"max_iterations\": " << config.stage6_max_iterations << ",\n";
    out << "    \"convergence_tolerance\": " << config.stage6_convergence_tolerance << ",\n";
    out << "    \"enforce_non_crossing\": " << (config.stage6_enforce_non_crossing ? "true" : "false") << ",\n";
    out << "    \"min_separation\": " << config.stage6_min_separation << ",\n";
    out << "    \"export_meshes\": " << (config.stage6_export_obj_meshes ? "true" : "false") << ",\n";
    out << "    \"mesh_export_format\": \"" << meshFormatLabel(config.stage6_mesh_export_format) << "\",\n";
    out << "    \"split_in_out_meshes\": " << (config.stage6_split_in_out_meshes ? "true" : "false") << "\n";
    out << "  },\n";
    out << "  \"stage7\": {\n";
    out << "    \"enabled\": " << (config.stage7_enabled ? "true" : "false") << ",\n";
    out << "    \"method\": \"" << stage7MethodLabel(config.stage7_method) << "\",\n";
    out << "    \"smoothing_weight\": " << config.stage7_smoothing_weight << ",\n";
    out << "    \"max_iterations\": " << config.stage7_max_iterations << ",\n";
    out << "    \"convergence_tolerance\": " << config.stage7_convergence_tolerance << ",\n";
    out << "    \"preserve_seed_values\": " << (config.stage7_preserve_seed_values ? "true" : "false") << ",\n";
    out << "    \"lambda\": " << config.stage7_lambda << ",\n";
    out << "    \"data_weight_seed\": " << config.stage7_data_weight_seed << ",\n";
    out << "    \"data_weight_interp\": " << config.stage7_data_weight_interp << ",\n";
    out << "    \"use_reliable_core_only_for_fit\": "
        << (config.stage7_use_reliable_core_only_for_fit ? "true" : "false") << ",\n";
    out << "    \"boundary_condition_mode\": \"" << stage7BoundaryConditionModeLabel(config.stage7_boundary_condition_mode)
        << "\",\n";
    out << "    \"solver_max_iterations\": " << config.stage7_solver_max_iterations << ",\n";
    out << "    \"solver_tolerance\": " << config.stage7_solver_tolerance << ",\n";
    out << "    \"enforce_non_crossing\": " << (config.stage7_enforce_non_crossing ? "true" : "false") << ",\n";
    out << "    \"min_separation\": " << config.stage7_min_separation << ",\n";
    out << "    \"export_meshes\": " << (config.stage7_export_meshes ? "true" : "false") << "\n";
    out << "  },\n";
    out << "  \"stage8\": {\n";
    out << "    \"enabled\": " << (config.stage8_enabled ? "true" : "false") << ",\n";
    out << "    \"fit_radius\": " << config.stage8_fit_radius << ",\n";
    out << "    \"min_points\": " << config.stage8_min_points << ",\n";
    out << "    \"max_points\": " << config.stage8_max_points << ",\n";
    out << "    \"max_rms_residual\": " << config.stage8_max_rms_residual << ",\n";
    out << "    \"max_abs_residual\": " << config.stage8_max_abs_residual << ",\n";
    out << "    \"max_condition_indicator\": " << config.stage8_max_condition_indicator << ",\n";
    out << "    \"require_centered_support\": " << (config.stage8_require_centered_support ? "true" : "false")
        << ",\n";
    out << "    \"min_directional_span\": " << config.stage8_min_directional_span << ",\n";
    out << "    \"export_csv\": " << (config.stage8_export_csv ? "true" : "false") << "\n";
    out << "  },\n";
    out << "  \"stage9\": {\n";
    out << "    \"enabled\": " << (config.stage9_enabled ? "true" : "false") << ",\n";
    out << "    \"export_csv\": " << (config.stage9_export_csv ? "true" : "false") << ",\n";
    out << "    \"qc_n_tail\": " << config.stage9_qc_n_tail << ",\n";
    out << "    \"qc_n_spike\": " << config.stage9_qc_n_spike << ",\n";
    out << "    \"qc_min_neighbors\": " << config.stage9_qc_min_neighbors << ",\n";
    out << "    \"qc_abs_scale_floor\": " << config.stage9_qc_abs_scale_floor << "\n";
    out << "  },\n";
    out << "  \"stage10\": {\n";
    out << "    \"enabled\": " << (config.stage10_enabled ? "true" : "false") << ",\n";
    out << "    \"thickness_method\": \"" << stage10ThicknessMethodLabel(config.stage10_thickness_method) << "\",\n";
    out << "    \"export_csv\": " << (config.stage10_export_csv ? "true" : "false") << ",\n";
    out << "    \"use_curvature_valid_domain_only\": "
        << (config.stage10_use_curvature_valid_domain_only ? "true" : "false") << ",\n";
    out << "    \"min_thickness\": " << config.stage10_min_thickness << ",\n";
    out << "    \"max_thickness\": " << config.stage10_max_thickness << "\n";
    out << "  },\n";
    out << "  \"parser\": {\n";
    out << "    \"include_hetatm\": " << (parser_config.include_hetatm ? "true" : "false") << ",\n";
    out << "    \"strict_mode\": " << (parser_config.strict_mode ? "true" : "false") << ",\n";
    out << "    \"keep_altloc_all\": " << (parser_config.keep_altloc_all ? "true" : "false") << ",\n";
    out << "    \"ignore_blank_chain_id\": " << (parser_config.ignore_blank_chain_id ? "true" : "false") << ",\n";
    out << "    \"verbose_warnings\": " << (parser_config.verbose_warnings ? "true" : "false") << ",\n";
    out << "    \"protein_only\": " << (parser_config.protein_only ? "true" : "false") << "\n";
    out << "  }\n";
    out << "}\n";
    return out.good();
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

enum class VdwResolutionSource : uint8_t { explicit_element = 0, inferred_from_atom_name = 1, fallback_unknown = 2 };

struct VdwResolution {
    std::string normalized_element;
    VdwResolutionSource source = VdwResolutionSource::fallback_unknown;
};

VdwResolution resolveVdwElement(const Atom& atom);

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
    out << "i,j,x,y,inside_disk,valid,candidate_patch_atom_count,inner_contact_serial_number,"
           "outer_contact_serial_number,inner_contact_patch_atom_index,outer_contact_patch_atom_index\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(result.valid_mask[idx])
                << ',' << result.candidate_patch_atom_counts[idx] << ',';
            if (result.valid_mask[idx] != 0) {
                out << result.inner_contact_serial_numbers[idx] << ',' << result.outer_contact_serial_numbers[idx] << ','
                    << result.inner_contact_patch_atom_indices[idx] << ','
                    << result.outer_contact_patch_atom_indices[idx];
            } else {
                out << "nan,nan,nan,nan";
            }
            out << '\n';
        }
    }
    return out.good();
}

std::string utcTimestampNowIso8601() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &tt);
#else
    gmtime_r(&tt, &tm);
#endif
    std::ostringstream out;
    out << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

bool writeStage4CsvBinaryMask(const GeometryStage4RawSheetResult& result,
                              const std::string& path,
                              const char* header_name,
                              const std::vector<uint8_t>& mask) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk," << header_name << "\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage4SummaryCsv(const GeometryStage4RawSheetResult& result) {
    std::ofstream out(result.summary_csv_path);
    if (!out) {
        return false;
    }
    out << "stage4_start_utc,stage4_end_utc,grid_spacing,cylinder_radius,nx,ny,total_nodes,inside_disk_nodes,"
           "valid_nodes,invalid_nodes,outer_only_nodes,inner_only_nodes,both_hit_nodes,zero_thickness_nodes,"
           "negative_thickness_nodes,unique_outer_contact_patch_atoms,unique_inner_contact_patch_atoms,"
           "unique_both_contact_patch_atoms,unique_contact_patch_atom_union,unique_outer_contact_serials,"
           "unique_inner_contact_serials,unique_both_contact_serials,unique_contact_serial_union,runtime_seconds\n";
    out << result.stage4_start_timestamp_utc << ',' << result.stage4_end_timestamp_utc << ',' << result.grid.spacing
        << ',' << result.grid.x_max << ',' << result.grid.nx << ',' << result.grid.ny << ',' << result.node_count
        << ',' << result.inside_disk_count << ',' << result.valid_node_count << ',' << result.invalid_node_count
        << ',' << result.outer_only_node_count << ',' << result.inner_only_node_count << ','
        << result.both_hit_node_count << ',' << result.zero_thickness_node_count << ','
        << result.negative_thickness_node_count << ',' << result.unique_outer_contact_atom_count << ','
        << result.unique_inner_contact_atom_count << ',' << result.unique_both_contact_atom_count << ','
        << result.unique_contact_atom_count << ',' << result.unique_outer_contact_serial_count << ','
        << result.unique_inner_contact_serial_count << ',' << result.unique_both_contact_serial_count << ','
        << result.unique_contact_serial_union_count << ',' << result.stage4_runtime_seconds << '\n';
    return out.good();
}

double stage5NodeRadius(const Stage4GridDescriptor& grid, std::size_t i, std::size_t j) {
    const double x = grid.x_values[i];
    const double y = grid.y_values[j];
    return std::sqrt((x * x) + (y * y));
}

bool writeStage5SeedCsv(const GeometryStage5SurfacePrepResult& result,
                        const std::string& path,
                        const char* value_header_name,
                        const std::vector<double>& values) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk,raw_valid,paired_seed," << value_header_name << "\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(result.raw_valid_mask[idx])
                << ',' << static_cast<int>(result.paired_seed_mask[idx]) << ',';
            if (result.paired_seed_mask[idx] != 0) {
                out << values[idx];
            } else {
                out << "nan";
            }
            out << '\n';
        }
    }
    return out.good();
}

bool writeStage5MaskCsv(const GeometryStage5SurfacePrepResult& result,
                        const std::string& path,
                        const char* mask_header_name,
                        const std::vector<uint8_t>& mask) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk," << mask_header_name << "\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage5SummaryCsv(const GeometryStage5SurfacePrepResult& result) {
    std::ofstream out(result.summary_csv_path);
    if (!out) {
        return false;
    }
    out << "boundary_margin,support_radius,reliable_radius,total_nodes,inside_disk_nodes,raw_valid_nodes,"
           "raw_invalid_nodes,outer_seed_nodes,inner_seed_nodes,paired_seed_nodes,boundary_excluded_nodes,"
           "interp_allowed_outer_nodes,interp_allowed_inner_nodes,paired_interp_allowed_nodes,hard_invalid_nodes,"
           "reliable_core_nodes,unique_outer_seed_patch_atoms,unique_inner_seed_patch_atoms,"
           "unique_seed_patch_atom_union,unique_outer_seed_serials,unique_inner_seed_serials\n";
    out << result.boundary_margin << ',' << result.support_radius << ',' << result.reliable_radius << ','
        << result.node_count << ',' << result.inside_disk_count << ',' << result.raw_valid_node_count << ','
        << result.raw_invalid_node_count << ',' << result.outer_seed_node_count << ',' << result.inner_seed_node_count
        << ',' << result.paired_seed_node_count << ',' << result.boundary_excluded_node_count << ','
        << result.interp_allowed_outer_node_count << ',' << result.interp_allowed_inner_node_count << ','
        << result.paired_interp_allowed_node_count << ',' << result.hard_invalid_node_count << ','
        << result.reliable_core_node_count << ',' << result.unique_outer_seed_patch_atom_count << ','
        << result.unique_inner_seed_patch_atom_count << ',' << result.unique_seed_patch_atom_index_union_count << ','
        << result.unique_outer_seed_atom_serials.size() << ',' << result.unique_inner_seed_atom_serials.size() << '\n';
    return out.good();
}

struct Stage6FieldSolveResult {
    std::vector<double> values;
    std::vector<uint8_t> solved_interp_mask;
    std::size_t iterations_used = 0;
    double final_max_update = 0.0;
};

struct Stage7FieldSmoothResult {
    std::vector<double> values;
    std::size_t iterations_used = 0;
    double final_max_update = 0.0;
};

struct Stage7FieldFitResult {
    std::vector<double> values;
    std::size_t iterations_used = 0;
    double final_update = 0.0;
    double final_residual = 0.0;
    double max_abs_residual = 0.0;
    double bending_energy = 0.0;
    std::size_t active_node_count = 0;
};

std::vector<std::size_t> stage6FourNeighborIndices(std::size_t i, std::size_t j, std::size_t nx, std::size_t ny) {
    std::vector<std::size_t> neighbors;
    neighbors.reserve(4);
    if (i > 0) {
        neighbors.push_back(stage4NodeIndex(i - 1, j, nx));
    }
    if ((i + 1) < nx) {
        neighbors.push_back(stage4NodeIndex(i + 1, j, nx));
    }
    if (j > 0) {
        neighbors.push_back(stage4NodeIndex(i, j - 1, nx));
    }
    if ((j + 1) < ny) {
        neighbors.push_back(stage4NodeIndex(i, j + 1, nx));
    }
    return neighbors;
}

Stage6FieldSolveResult runStage6FieldReconstruction(const Stage4GridDescriptor& grid,
                                                     const std::vector<double>& seed_values,
                                                     const std::vector<uint8_t>& inside_disk_mask,
                                                     const std::vector<uint8_t>& paired_seed_mask,
                                                     const std::vector<uint8_t>& paired_interp_allowed_mask,
                                                     const std::vector<uint8_t>& hard_invalid_mask,
                                                     const FoldPatchAnalysisConfig& config) {
    const std::size_t node_count = grid.nx * grid.ny;
    Stage6FieldSolveResult solve_result;
    solve_result.values.assign(node_count, std::numeric_limits<double>::quiet_NaN());
    solve_result.solved_interp_mask.assign(node_count, 0);

    // Stage 6 v1 algorithm:
    // 1) keep paired seeds as immutable anchors,
    // 2) initialize interpolation-allowed nodes from local 4-neighbor averages when possible,
    // 3) perform deterministic Jacobi relaxation using only 4-neighbor finite values.
    const double alpha = config.stage6_smoothing_weight / (1.0 + config.stage6_smoothing_weight);
    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (paired_seed_mask[idx] != 0) {
            solve_result.values[idx] = seed_values[idx];
        }
    }

    for (std::size_t j = 0; j < grid.ny; ++j) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
            if (paired_interp_allowed_mask[idx] == 0 || inside_disk_mask[idx] == 0 || hard_invalid_mask[idx] != 0) {
                continue;
            }
            const auto neighbors = stage6FourNeighborIndices(i, j, grid.nx, grid.ny);
            double sum = 0.0;
            std::size_t count = 0;
            for (const std::size_t nidx : neighbors) {
                const double value = solve_result.values[nidx];
                if (std::isfinite(value)) {
                    sum += value;
                    ++count;
                }
            }
            if (count > 0) {
                solve_result.values[idx] = sum / static_cast<double>(count);
            }
        }
    }

    std::vector<double> current = solve_result.values;
    std::vector<double> next = current;
    std::size_t iterations = 0;
    double final_max_update = 0.0;
    for (; iterations < config.stage6_max_iterations; ++iterations) {
        double iteration_max_update = 0.0;
        bool any_updated = false;
        next = current;

        for (std::size_t j = 0; j < grid.ny; ++j) {
            for (std::size_t i = 0; i < grid.nx; ++i) {
                const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
                if (paired_interp_allowed_mask[idx] == 0 || inside_disk_mask[idx] == 0 || hard_invalid_mask[idx] != 0) {
                    continue;
                }

                const auto neighbors = stage6FourNeighborIndices(i, j, grid.nx, grid.ny);
                double sum = 0.0;
                std::size_t count = 0;
                for (const std::size_t nidx : neighbors) {
                    const double value = current[nidx];
                    if (std::isfinite(value)) {
                        sum += value;
                        ++count;
                    }
                }
                if (count == 0) {
                    continue;
                }

                const double average = sum / static_cast<double>(count);
                const double previous = current[idx];
                if (!std::isfinite(previous)) {
                    next[idx] = average;
                    iteration_max_update = std::max(iteration_max_update, std::fabs(average));
                    any_updated = true;
                    continue;
                }
                const double updated = previous + (alpha * (average - previous));
                next[idx] = updated;
                iteration_max_update = std::max(iteration_max_update, std::fabs(updated - previous));
                any_updated = true;
            }
        }

        current.swap(next);
        final_max_update = iteration_max_update;
        if (!any_updated || iteration_max_update < config.stage6_convergence_tolerance) {
            ++iterations;
            break;
        }
    }

    solve_result.values = std::move(current);
    solve_result.iterations_used = iterations;
    solve_result.final_max_update = final_max_update;
    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (paired_interp_allowed_mask[idx] != 0 && std::isfinite(solve_result.values[idx])) {
            solve_result.solved_interp_mask[idx] = 1;
        }
    }
    return solve_result;
}

bool writeStage6FieldCsv(const GeometryStage6SurfaceReconstructionResult& result,
                         const std::string& path,
                         const char* value_header_name,
                         const std::vector<double>& values) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk,paired_seed,paired_interp_allowed,hard_invalid,reconstructed," << value_header_name
        << "\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(result.paired_seed_mask[idx])
                << ',' << static_cast<int>(result.paired_interp_allowed_mask[idx]) << ','
                << static_cast<int>(result.hard_invalid_mask[idx]) << ',' << static_cast<int>(result.reconstructed_mask[idx])
                << ',';
            if (std::isfinite(values[idx])) {
                out << values[idx];
            } else {
                out << "nan";
            }
            out << '\n';
        }
    }
    return out.good();
}

bool writeStage6MaskCsv(const GeometryStage6SurfaceReconstructionResult& result,
                        const std::string& path,
                        const char* header_name,
                        const std::vector<uint8_t>& mask) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,inside_disk," << header_name << "\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.inside_disk_mask[idx]) << ',' << static_cast<int>(mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage6SummaryCsv(const GeometryStage6SurfaceReconstructionResult& result) {
    std::ofstream out(result.summary_csv_path);
    if (!out) {
        return false;
    }
    out << "node_count,seed_nodes,interp_nodes,reconstructed_nodes,final_valid_analysis_nodes,unresolved_nodes,"
           "outer_iterations_used,inner_iterations_used,outer_final_max_update,inner_final_max_update,"
           "non_crossing_adjusted_nodes,min_reconstructed_separation,max_reconstructed_separation,"
           "mean_reconstructed_separation,outer_obj_vertex_count,outer_obj_face_count,inner_obj_vertex_count,"
           "inner_obj_face_count,outer_reconstructed_scalar_node_count,outer_reconstructed_nodes_not_used_by_any_face,"
           "inner_reconstructed_scalar_node_count,inner_reconstructed_nodes_not_used_by_any_face\n";
    out << result.node_count << ',' << result.seed_node_count << ',' << result.interp_node_count << ','
        << result.reconstructed_node_count << ',' << result.final_valid_analysis_node_count << ','
        << result.unresolved_node_count << ',' << result.outer_iterations_used << ',' << result.inner_iterations_used
        << ',' << result.outer_final_max_update << ',' << result.inner_final_max_update << ','
        << result.non_crossing_adjusted_node_count << ',' << result.min_reconstructed_separation << ','
        << result.max_reconstructed_separation << ',' << result.mean_reconstructed_separation << ','
        << result.outer_obj_vertex_count << ',' << result.outer_obj_face_count << ',' << result.inner_obj_vertex_count
        << ',' << result.inner_obj_face_count << ',' << result.outer_reconstructed_scalar_node_count << ','
        << result.outer_reconstructed_nodes_not_used_by_any_face << ','
        << result.inner_reconstructed_scalar_node_count << ','
        << result.inner_reconstructed_nodes_not_used_by_any_face << '\n';
    return out.good();
}

Stage7FieldSmoothResult runStage7FieldSmoothingLegacy(const Stage4GridDescriptor& grid,
                                                      const std::vector<double>& reconstructed_values,
                                                      const std::vector<uint8_t>& reconstructed_mask,
                                                      const std::vector<uint8_t>& paired_seed_mask,
                                                      const FoldPatchAnalysisConfig& config) {
    Stage7FieldSmoothResult smooth_result;
    smooth_result.values = reconstructed_values;
    const double alpha = config.stage7_smoothing_weight / (1.0 + config.stage7_smoothing_weight);

    std::vector<double> current = smooth_result.values;
    std::vector<double> next = current;
    std::size_t iterations = 0;
    double final_max_update = 0.0;
    for (; iterations < config.stage7_max_iterations; ++iterations) {
        double iteration_max_update = 0.0;
        bool any_updated = false;
        next = current;

        for (std::size_t j = 0; j < grid.ny; ++j) {
            for (std::size_t i = 0; i < grid.nx; ++i) {
                const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
                if (reconstructed_mask[idx] == 0 || !std::isfinite(current[idx])) {
                    continue;
                }
                // ON  (pin seed): paired seed nodes remain fixed at their Stage 6 values.
                // OFF (allow move): paired seed nodes are updated like all other reconstructed nodes.
                if (config.stage7_preserve_seed_values && paired_seed_mask[idx] != 0) {
                    continue;
                }
                const auto neighbors = stage6FourNeighborIndices(i, j, grid.nx, grid.ny);
                double sum = 0.0;
                std::size_t count = 0;
                for (const std::size_t nidx : neighbors) {
                    if (reconstructed_mask[nidx] == 0 || !std::isfinite(current[nidx])) {
                        continue;
                    }
                    sum += current[nidx];
                    ++count;
                }
                if (count == 0) {
                    continue;
                }
                const double previous = current[idx];
                const double average = sum / static_cast<double>(count);
                const double updated = previous + (alpha * (average - previous));
                next[idx] = updated;
                iteration_max_update = std::max(iteration_max_update, std::fabs(updated - previous));
                any_updated = true;
            }
        }

        current.swap(next);
        final_max_update = iteration_max_update;
        if (!any_updated || iteration_max_update < config.stage7_convergence_tolerance) {
            ++iterations;
            break;
        }
    }

    smooth_result.values = std::move(current);
    smooth_result.iterations_used = iterations;
    smooth_result.final_max_update = final_max_update;
    return smooth_result;
}

std::vector<uint8_t> buildStage7FitDomainMask(const GeometryStage6SurfaceReconstructionResult& stage6_result,
                                              const GeometryStage5SurfacePrepResult& stage5_result,
                                              const FoldPatchAnalysisConfig& config) {
    std::vector<uint8_t> fit_domain_mask(stage6_result.node_count, 0);
    for (std::size_t idx = 0; idx < stage6_result.node_count; ++idx) {
        if (stage6_result.reconstructed_mask[idx] == 0 || !std::isfinite(stage6_result.z_outer_reconstructed[idx]) ||
            !std::isfinite(stage6_result.z_inner_reconstructed[idx])) {
            continue;
        }
        if (config.stage7_use_reliable_core_only_for_fit && stage5_result.reliable_core_mask[idx] == 0) {
            continue;
        }
        fit_domain_mask[idx] = 1;
    }
    return fit_domain_mask;
}

std::vector<double> buildStage7FidelityWeights(const std::vector<uint8_t>& fit_domain_mask,
                                               const std::vector<uint8_t>& paired_seed_mask,
                                               const FoldPatchAnalysisConfig& config,
                                               std::size_t& seed_like_count,
                                               std::size_t& interp_like_count) {
    std::vector<double> weights(fit_domain_mask.size(), 0.0);
    seed_like_count = 0;
    interp_like_count = 0;
    for (std::size_t idx = 0; idx < fit_domain_mask.size(); ++idx) {
        if (fit_domain_mask[idx] == 0) {
            continue;
        }
        if (paired_seed_mask[idx] != 0) {
            weights[idx] = config.stage7_data_weight_seed;
            ++seed_like_count;
        } else {
            weights[idx] = config.stage7_data_weight_interp;
            ++interp_like_count;
        }
    }
    return weights;
}

std::vector<uint8_t> buildStage7BoundaryMask(const Stage4GridDescriptor& grid, const std::vector<uint8_t>& fit_domain_mask) {
    std::vector<uint8_t> boundary_mask(grid.nx * grid.ny, 0);
    for (std::size_t j = 0; j < grid.ny; ++j) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
            if (fit_domain_mask[idx] == 0) {
                continue;
            }
            bool is_boundary = false;
            const auto neighbors = stage6FourNeighborIndices(i, j, grid.nx, grid.ny);
            for (const std::size_t nidx : neighbors) {
                if (fit_domain_mask[nidx] == 0) {
                    is_boundary = true;
                    break;
                }
            }
            if (i == 0 || j == 0 || i + 1 == grid.nx || j + 1 == grid.ny) {
                is_boundary = true;
            }
            boundary_mask[idx] = is_boundary ? 1 : 0;
        }
    }
    return boundary_mask;
}

void applyStage7Laplacian(const Stage4GridDescriptor& grid,
                          const std::vector<uint8_t>& fit_domain_mask,
                          const std::vector<double>& values,
                          std::vector<double>& out) {
    out.assign(grid.nx * grid.ny, 0.0);
    for (std::size_t j = 0; j < grid.ny; ++j) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
            if (fit_domain_mask[idx] == 0) {
                continue;
            }
            const auto neighbors = stage6FourNeighborIndices(i, j, grid.nx, grid.ny);
            double sum_neighbors = 0.0;
            std::size_t count = 0;
            for (const std::size_t nidx : neighbors) {
                if (fit_domain_mask[nidx] == 0) {
                    continue;
                }
                sum_neighbors += values[nidx];
                ++count;
            }
            out[idx] = (static_cast<double>(count) * values[idx]) - sum_neighbors;
        }
    }
}

double computeDiscreteStage7BendingEnergy(const Stage4GridDescriptor& grid,
                                          const std::vector<uint8_t>& fit_domain_mask,
                                          const std::vector<double>& values) {
    std::vector<double> laplacian;
    applyStage7Laplacian(grid, fit_domain_mask, values, laplacian);
    double energy = 0.0;
    for (std::size_t idx = 0; idx < laplacian.size(); ++idx) {
        if (fit_domain_mask[idx] == 0) {
            continue;
        }
        energy += laplacian[idx] * laplacian[idx];
    }
    return energy;
}

void applyStage7System(const Stage4GridDescriptor& grid,
                       const std::vector<uint8_t>& fit_domain_mask,
                       const std::vector<double>& fidelity_weights,
                       double lambda,
                       const std::vector<double>& x,
                       std::vector<double>& out) {
    std::vector<double> lx;
    std::vector<double> ltlx;
    applyStage7Laplacian(grid, fit_domain_mask, x, lx);
    applyStage7Laplacian(grid, fit_domain_mask, lx, ltlx);
    out.assign(grid.nx * grid.ny, 0.0);
    for (std::size_t idx = 0; idx < out.size(); ++idx) {
        if (fit_domain_mask[idx] == 0) {
            continue;
        }
        out[idx] = (fidelity_weights[idx] * x[idx]) + (lambda * ltlx[idx]);
    }
}

Stage7FieldFitResult runStage7FieldThinPlateGridFit(const Stage4GridDescriptor& grid,
                                                     const std::vector<double>& reconstructed_values,
                                                     const std::vector<uint8_t>& fit_domain_mask,
                                                     const std::vector<uint8_t>& boundary_mask,
                                                     const std::vector<double>& fidelity_weights,
                                                     const FoldPatchAnalysisConfig& config) {
    const std::size_t node_count = grid.nx * grid.ny;
    Stage7FieldFitResult fit_result;
    fit_result.values = reconstructed_values;
    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (fit_domain_mask[idx] != 0) {
            ++fit_result.active_node_count;
        }
    }
    const double soft_boundary_factor = 5.0;
    std::vector<double> effective_weights = fidelity_weights;
    std::vector<uint8_t> fixed_mask(node_count, 0);
    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (fit_domain_mask[idx] == 0) {
            continue;
        }
        if (config.stage7_boundary_condition_mode == FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::fixed_to_stage6 &&
            boundary_mask[idx] != 0) {
            fixed_mask[idx] = 1;
            effective_weights[idx] = std::max(effective_weights[idx], 1.0);
        } else if (config.stage7_boundary_condition_mode ==
                       FoldPatchAnalysisConfig::Stage7BoundaryConditionMode::soft_to_stage6 &&
                   boundary_mask[idx] != 0) {
            effective_weights[idx] *= soft_boundary_factor;
        }
    }

    std::vector<double> b(node_count, 0.0);
    std::vector<double> x = fit_result.values;
    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (fit_domain_mask[idx] == 0) {
            x[idx] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        if (!std::isfinite(x[idx])) {
            x[idx] = reconstructed_values[idx];
        }
        if (fixed_mask[idx] == 0) {
            b[idx] = effective_weights[idx] * reconstructed_values[idx];
        }
    }

    std::vector<double> ax(node_count, 0.0);
    applyStage7System(grid, fit_domain_mask, effective_weights, config.stage7_lambda, x, ax);
    std::vector<double> r(node_count, 0.0);
    std::vector<double> p(node_count, 0.0);
    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (fit_domain_mask[idx] == 0 || fixed_mask[idx] != 0) {
            continue;
        }
        r[idx] = b[idx] - ax[idx];
        p[idx] = r[idx];
    }

    auto maskedDot = [&](const std::vector<double>& a, const std::vector<double>& bvec) {
        double dot = 0.0;
        for (std::size_t idx = 0; idx < node_count; ++idx) {
            if (fit_domain_mask[idx] == 0 || fixed_mask[idx] != 0) {
                continue;
            }
            dot += a[idx] * bvec[idx];
        }
        return dot;
    };

    double rr = maskedDot(r, r);
    std::size_t iterations = 0;
    double final_update = 0.0;
    for (; iterations < config.stage7_solver_max_iterations; ++iterations) {
        if (rr <= (config.stage7_solver_tolerance * config.stage7_solver_tolerance)) {
            ++iterations;
            break;
        }
        std::vector<double> ap(node_count, 0.0);
        applyStage7System(grid, fit_domain_mask, effective_weights, config.stage7_lambda, p, ap);
        const double pap = maskedDot(p, ap);
        if (pap <= 0.0) {
            break;
        }
        const double alpha = rr / pap;
        double iteration_max_update = 0.0;
        for (std::size_t idx = 0; idx < node_count; ++idx) {
            if (fit_domain_mask[idx] == 0 || fixed_mask[idx] != 0) {
                continue;
            }
            const double delta = alpha * p[idx];
            x[idx] += delta;
            r[idx] -= alpha * ap[idx];
            iteration_max_update = std::max(iteration_max_update, std::fabs(delta));
        }
        final_update = iteration_max_update;
        const double rr_new = maskedDot(r, r);
        if (rr_new <= (config.stage7_solver_tolerance * config.stage7_solver_tolerance)) {
            rr = rr_new;
            ++iterations;
            break;
        }
        const double beta = rr_new / rr;
        rr = rr_new;
        for (std::size_t idx = 0; idx < node_count; ++idx) {
            if (fit_domain_mask[idx] == 0 || fixed_mask[idx] != 0) {
                continue;
            }
            p[idx] = r[idx] + (beta * p[idx]);
        }
    }

    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (fit_domain_mask[idx] == 0) {
            x[idx] = std::numeric_limits<double>::quiet_NaN();
        } else if (fixed_mask[idx] != 0) {
            x[idx] = reconstructed_values[idx];
        }
    }
    std::vector<double> residual_system(node_count, 0.0);
    applyStage7System(grid, fit_domain_mask, effective_weights, config.stage7_lambda, x, residual_system);
    double residual_sq_sum = 0.0;
    double max_abs_residual = 0.0;
    std::size_t residual_count = 0;
    for (std::size_t idx = 0; idx < node_count; ++idx) {
        if (fit_domain_mask[idx] == 0) {
            continue;
        }
        const double rhs = fixed_mask[idx] != 0 ? reconstructed_values[idx] : b[idx];
        const double residual = residual_system[idx] - rhs;
        residual_sq_sum += residual * residual;
        max_abs_residual = std::max(max_abs_residual, std::fabs(residual));
        ++residual_count;
    }
    fit_result.values = std::move(x);
    fit_result.iterations_used = iterations;
    fit_result.final_update = final_update;
    fit_result.final_residual =
        residual_count > 0 ? std::sqrt(residual_sq_sum / static_cast<double>(residual_count)) : 0.0;
    fit_result.max_abs_residual = max_abs_residual;
    fit_result.bending_energy = computeDiscreteStage7BendingEnergy(grid, fit_domain_mask, fit_result.values);
    return fit_result;
}

bool writeStage7FieldCsv(const GeometryStage7SmoothedSurfaceResult& result,
                         const std::string& path,
                         const char* value_header_name,
                         const std::vector<double>& values) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,reconstructed,reliable_core,metric_domain," << value_header_name << "\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.reconstructed_mask[idx]) << ','
                << static_cast<int>(result.reliable_core_mask[idx]) << ','
                << static_cast<int>(result.metric_domain_mask[idx]) << ',';
            if (std::isfinite(values[idx])) {
                out << values[idx];
            } else {
                out << "nan";
            }
            out << '\n';
        }
    }
    return out.good();
}

bool writeStage7MaskCsv(const GeometryStage7SmoothedSurfaceResult& result,
                        const std::string& path,
                        const char* mask_name,
                        const std::vector<uint8_t>& mask) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y," << mask_name << "\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(mask[idx]) << '\n';
        }
    }
    return out.good();
}

enum class Stage7DeltaCsvKind : uint8_t { outer = 0, inner = 1, thickness = 2 };

bool writeStage7DeltaCsv(const Stage4GridDescriptor& grid,
                         const std::vector<uint8_t>& inside_disk_mask,
                         const GeometryStage5SurfacePrepResult& stage5_result,
                         const GeometryStage6SurfaceReconstructionResult& stage6_result,
                         const GeometryStage7SmoothedSurfaceResult& stage7_result,
                         Stage7DeltaCsvKind kind,
                         const std::string& path,
                         double& max_abs_delta,
                         double& mean_abs_delta) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    if (kind == Stage7DeltaCsvKind::outer) {
        out << "i,j,x,y,inside_disk,reconstructed,smooth_valid,metric_domain,paired_seed,interp_allowed,hard_invalid,"
               "reliable_core,z_outer_stage6,z_outer_stage7,delta_outer\n";
    } else if (kind == Stage7DeltaCsvKind::inner) {
        out << "i,j,x,y,inside_disk,reconstructed,smooth_valid,metric_domain,paired_seed,interp_allowed,hard_invalid,"
               "reliable_core,z_inner_stage6,z_inner_stage7,delta_inner\n";
    } else {
        out << "i,j,x,y,inside_disk,reconstructed,smooth_valid,metric_domain,paired_seed,interp_allowed,hard_invalid,"
               "reliable_core,t_stage6,t_stage7,delta_t\n";
    }

    max_abs_delta = 0.0;
    double sum_abs_delta = 0.0;
    std::size_t finite_delta_count = 0;
    for (std::size_t j = 0; j < grid.ny; ++j) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
            out << i << ',' << j << ',' << grid.x_values[i] << ',' << grid.y_values[j] << ','
                << static_cast<int>(inside_disk_mask[idx]) << ','
                << static_cast<int>(stage6_result.reconstructed_mask[idx]) << ','
                << static_cast<int>(stage7_result.smooth_valid_mask[idx]) << ','
                << static_cast<int>(stage7_result.metric_domain_mask[idx]) << ','
                << static_cast<int>(stage5_result.paired_seed_mask[idx]) << ','
                << static_cast<int>(stage5_result.paired_interp_allowed_mask[idx]) << ','
                << static_cast<int>(stage5_result.hard_invalid_mask[idx]) << ','
                << static_cast<int>(stage5_result.reliable_core_mask[idx]) << ',';

            double stage6_value = std::numeric_limits<double>::quiet_NaN();
            double stage7_value = std::numeric_limits<double>::quiet_NaN();
            double delta_value = std::numeric_limits<double>::quiet_NaN();
            if (kind == Stage7DeltaCsvKind::outer) {
                stage6_value = stage6_result.z_outer_reconstructed[idx];
                stage7_value = stage7_result.z_outer_smooth[idx];
            } else if (kind == Stage7DeltaCsvKind::inner) {
                stage6_value = stage6_result.z_inner_reconstructed[idx];
                stage7_value = stage7_result.z_inner_smooth[idx];
            } else {
                const double z_outer_stage6 = stage6_result.z_outer_reconstructed[idx];
                const double z_inner_stage6 = stage6_result.z_inner_reconstructed[idx];
                const double z_outer_stage7 = stage7_result.z_outer_smooth[idx];
                const double z_inner_stage7 = stage7_result.z_inner_smooth[idx];
                if (std::isfinite(z_outer_stage6) && std::isfinite(z_inner_stage6)) {
                    stage6_value = z_outer_stage6 - z_inner_stage6;
                }
                if (std::isfinite(z_outer_stage7) && std::isfinite(z_inner_stage7)) {
                    stage7_value = z_outer_stage7 - z_inner_stage7;
                }
            }
            if (std::isfinite(stage6_value) && std::isfinite(stage7_value)) {
                delta_value = stage7_value - stage6_value;
            }

            if (std::isfinite(stage6_value)) {
                out << stage6_value;
            } else {
                out << "nan";
            }
            out << ',';
            if (std::isfinite(stage7_value)) {
                out << stage7_value;
            } else {
                out << "nan";
            }
            out << ',';
            if (std::isfinite(delta_value)) {
                out << delta_value;
                const double abs_delta = std::fabs(delta_value);
                max_abs_delta = std::max(max_abs_delta, abs_delta);
                sum_abs_delta += abs_delta;
                ++finite_delta_count;
            } else {
                out << "nan";
            }
            out << '\n';
        }
    }
    mean_abs_delta = finite_delta_count > 0 ? (sum_abs_delta / static_cast<double>(finite_delta_count)) : 0.0;
    return out.good();
}

bool writeStage7SummaryCsv(const GeometryStage7SmoothedSurfaceResult& result, const FoldPatchAnalysisConfig& config) {
    std::ofstream out(result.summary_csv_path);
    if (!out) {
        return false;
    }
    out << "stage7_enabled,stage7_method,smoothing_weight,max_iterations,convergence_tolerance,preserve_seed_values,"
           "stage7_lambda,stage7_data_weight_seed,stage7_data_weight_interp,stage7_data_weight_policy,"
           "stage7_use_reliable_core_only_for_fit,stage7_boundary_condition_mode,stage7_fit_active_node_count,"
           "stage7_fit_seed_like_node_count,stage7_fit_interp_like_node_count,stage7_solver_max_iterations,"
           "stage7_solver_tolerance,enforce_non_crossing,min_separation,node_count,smooth_valid_nodes,metric_domain_nodes,"
           "non_crossing_adjusted_nodes,outer_iterations_used,inner_iterations_used,outer_final_max_update,"
           "inner_final_max_update,outer_fit_final_residual,inner_fit_final_residual,outer_fit_max_abs_residual,"
           "inner_fit_max_abs_residual,outer_fit_bending_energy,inner_fit_bending_energy,outer_solver_iterations_used,"
           "inner_solver_iterations_used,outer_solver_final_update,inner_solver_final_update,min_smooth_separation,max_smooth_separation,mean_smooth_separation,"
           "max_abs_outer_delta,mean_abs_outer_delta,max_abs_inner_delta,mean_abs_inner_delta,max_abs_thickness_delta,mean_abs_thickness_delta,"
           "normal_orientation_outer,normal_orientation_inner,local_thickness_definition,metric_surface_definition,"
           "outer_mesh_vertex_count,outer_mesh_face_count,inner_mesh_vertex_count,inner_mesh_face_count,"
           "outer_mesh_unused_scalar_nodes,inner_mesh_unused_scalar_nodes\n";
    out << (config.stage7_enabled ? 1 : 0) << ',' << result.stage7_method_label << ',' << config.stage7_smoothing_weight
        << ','
        << config.stage7_max_iterations << ',' << config.stage7_convergence_tolerance << ','
        << (config.stage7_preserve_seed_values ? 1 : 0) << ',' << result.stage7_lambda << ','
        << result.stage7_data_weight_seed << ',' << result.stage7_data_weight_interp << ','
        << result.stage7_data_weight_policy << ',' << (result.stage7_use_reliable_core_only_for_fit ? 1 : 0) << ','
        << result.stage7_boundary_condition_mode_label << ',' << result.stage7_fit_active_node_count << ','
        << result.stage7_fit_seed_like_node_count << ',' << result.stage7_fit_interp_like_node_count << ','
        << config.stage7_solver_max_iterations << ',' << config.stage7_solver_tolerance << ','
        << (config.stage7_enforce_non_crossing ? 1 : 0) << ','
        << config.stage7_min_separation << ',' << result.node_count << ',' << result.smooth_valid_node_count << ','
        << result.metric_domain_node_count << ',' << result.smooth_non_crossing_adjusted_node_count << ','
        << result.outer_iterations_used << ',' << result.inner_iterations_used << ',' << result.outer_final_max_update
        << ',' << result.inner_final_max_update << ',' << result.outer_fit_final_residual << ','
        << result.inner_fit_final_residual << ',' << result.outer_fit_max_abs_residual << ','
        << result.inner_fit_max_abs_residual << ',' << result.outer_fit_bending_energy << ','
        << result.inner_fit_bending_energy << ',' << result.outer_solver_iterations_used << ','
        << result.inner_solver_iterations_used << ',' << result.outer_solver_final_update << ','
        << result.inner_solver_final_update << ',' << result.min_smooth_separation << ','
        << result.max_smooth_separation << ',' << result.mean_smooth_separation << ','
        << result.max_abs_outer_delta << ',' << result.mean_abs_outer_delta << ',' << result.max_abs_inner_delta << ','
        << result.mean_abs_inner_delta << ',' << result.max_abs_thickness_delta << ','
        << result.mean_abs_thickness_delta << ',' << result.normal_orientation_outer
        << ',' << result.normal_orientation_inner << ',' << result.local_thickness_definition << ','
        << result.metric_surface_definition << ',' << result.outer_mesh_vertex_count << ','
        << result.outer_mesh_face_count << ',' << result.inner_mesh_vertex_count << ',' << result.inner_mesh_face_count
        << ',' << result.outer_mesh_unused_scalar_nodes << ',' << result.inner_mesh_unused_scalar_nodes << '\n';
    return out.good();
}

std::size_t countStage6ReconstructedScalarNodes(const Stage4GridDescriptor& grid,
                                                const std::vector<uint8_t>& obj_vertex_mask,
                                                const std::vector<double>& values) {
    std::size_t count = 0;
    for (std::size_t idx = 0; idx < (grid.nx * grid.ny); ++idx) {
        if (obj_vertex_mask[idx] != 0 && std::isfinite(values[idx])) {
            ++count;
        }
    }
    return count;
}

std::size_t countStage6ReconstructedNodesUsedByFaces(const Stage4GridDescriptor& grid,
                                                     const std::vector<uint8_t>& obj_vertex_mask,
                                                     const std::vector<double>& values) {
    std::vector<uint8_t> used(grid.nx * grid.ny, 0);
    for (std::size_t j = 0; (j + 1) < grid.ny; ++j) {
        for (std::size_t i = 0; (i + 1) < grid.nx; ++i) {
            const std::size_t idx00 = stage4NodeIndex(i, j, grid.nx);
            const std::size_t idx10 = stage4NodeIndex(i + 1, j, grid.nx);
            const std::size_t idx01 = stage4NodeIndex(i, j + 1, grid.nx);
            const std::size_t idx11 = stage4NodeIndex(i + 1, j + 1, grid.nx);
            if (obj_vertex_mask[idx00] == 0 || obj_vertex_mask[idx10] == 0 || obj_vertex_mask[idx01] == 0 ||
                obj_vertex_mask[idx11] == 0 || !std::isfinite(values[idx00]) || !std::isfinite(values[idx10]) ||
                !std::isfinite(values[idx01]) || !std::isfinite(values[idx11])) {
                continue;
            }
            used[idx00] = 1;
            used[idx10] = 1;
            used[idx01] = 1;
            used[idx11] = 1;
        }
    }
    std::size_t count = 0;
    for (const uint8_t flag : used) {
        if (flag != 0) {
            ++count;
        }
    }
    return count;
}

struct Stage6ObjExportResult {
    std::size_t vertex_count = 0;
    std::size_t face_count = 0;
};

struct Stage6StlExportResult {
    std::size_t vertex_count = 0;
    std::size_t face_count = 0;
};

struct Stage6DualMeshExportResult {
    std::size_t outer_vertex_count = 0;
    std::size_t outer_face_count = 0;
    std::size_t inner_vertex_count = 0;
    std::size_t inner_face_count = 0;
};

Stage6ObjExportResult writeStage6ObjMesh(const Stage4GridDescriptor& grid,
                                         const std::vector<uint8_t>& obj_vertex_mask,
                                         const std::vector<double>& values,
                                         const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open Stage 6 OBJ path for writing: " + path);
    }

    Stage6ObjExportResult export_result;
    std::vector<int> vertex_indices(grid.nx * grid.ny, -1);
    int next_vertex_index = 1;
    for (std::size_t j = 0; j < grid.ny; ++j) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
            if (obj_vertex_mask[idx] == 0 || !std::isfinite(values[idx])) {
                continue;
            }
            out << "v " << grid.x_values[i] << ' ' << grid.y_values[j] << ' ' << values[idx] << '\n';
            vertex_indices[idx] = next_vertex_index;
            ++next_vertex_index;
            ++export_result.vertex_count;
        }
    }

    for (std::size_t j = 0; (j + 1) < grid.ny; ++j) {
        for (std::size_t i = 0; (i + 1) < grid.nx; ++i) {
            const std::size_t idx00 = stage4NodeIndex(i, j, grid.nx);
            const std::size_t idx10 = stage4NodeIndex(i + 1, j, grid.nx);
            const std::size_t idx01 = stage4NodeIndex(i, j + 1, grid.nx);
            const std::size_t idx11 = stage4NodeIndex(i + 1, j + 1, grid.nx);
            if (obj_vertex_mask[idx00] == 0 || obj_vertex_mask[idx10] == 0 || obj_vertex_mask[idx01] == 0 ||
                obj_vertex_mask[idx11] == 0 || !std::isfinite(values[idx00]) || !std::isfinite(values[idx10]) ||
                !std::isfinite(values[idx01]) || !std::isfinite(values[idx11])) {
                continue;
            }
            const int v00 = vertex_indices[idx00];
            const int v10 = vertex_indices[idx10];
            const int v01 = vertex_indices[idx01];
            const int v11 = vertex_indices[idx11];
            out << "f " << v00 << ' ' << v10 << ' ' << v01 << '\n';
            out << "f " << v10 << ' ' << v11 << ' ' << v01 << '\n';
            export_result.face_count += 2;
        }
    }
    if (!out.good()) {
        throw std::runtime_error("Failed while writing Stage 6 OBJ mesh: " + path);
    }
    return export_result;
}

Stage6StlExportResult writeStage6StlMesh(const Stage4GridDescriptor& grid,
                                         const std::vector<uint8_t>& obj_vertex_mask,
                                         const std::vector<double>& values,
                                         const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open Stage 6 STL path for writing: " + path);
    }

    const auto writeFacet = [&](double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2) {
        const double ux = x1 - x0;
        const double uy = y1 - y0;
        const double uz = z1 - z0;
        const double vx = x2 - x0;
        const double vy = y2 - y0;
        const double vz = z2 - z0;
        double nx = (uy * vz) - (uz * vy);
        double ny = (uz * vx) - (ux * vz);
        double nz = (ux * vy) - (uy * vx);
        const double norm = std::sqrt((nx * nx) + (ny * ny) + (nz * nz));
        if (norm > 0.0) {
            nx /= norm;
            ny /= norm;
            nz /= norm;
        }
        out << "  facet normal " << nx << ' ' << ny << ' ' << nz << '\n'
            << "    outer loop\n"
            << "      vertex " << x0 << ' ' << y0 << ' ' << z0 << '\n'
            << "      vertex " << x1 << ' ' << y1 << ' ' << z1 << '\n'
            << "      vertex " << x2 << ' ' << y2 << ' ' << z2 << '\n'
            << "    endloop\n"
            << "  endfacet\n";
    };

    Stage6StlExportResult export_result;
    std::vector<uint8_t> vertex_used(grid.nx * grid.ny, 0);
    out << "solid capdat_stage6_mesh\n";
    for (std::size_t j = 0; (j + 1) < grid.ny; ++j) {
        for (std::size_t i = 0; (i + 1) < grid.nx; ++i) {
            const std::size_t idx00 = stage4NodeIndex(i, j, grid.nx);
            const std::size_t idx10 = stage4NodeIndex(i + 1, j, grid.nx);
            const std::size_t idx01 = stage4NodeIndex(i, j + 1, grid.nx);
            const std::size_t idx11 = stage4NodeIndex(i + 1, j + 1, grid.nx);
            if (obj_vertex_mask[idx00] == 0 || obj_vertex_mask[idx10] == 0 || obj_vertex_mask[idx01] == 0 ||
                obj_vertex_mask[idx11] == 0 || !std::isfinite(values[idx00]) || !std::isfinite(values[idx10]) ||
                !std::isfinite(values[idx01]) || !std::isfinite(values[idx11])) {
                continue;
            }

            writeFacet(grid.x_values[i], grid.y_values[j], values[idx00], grid.x_values[i + 1], grid.y_values[j], values[idx10],
                       grid.x_values[i], grid.y_values[j + 1], values[idx01]);
            writeFacet(grid.x_values[i + 1], grid.y_values[j], values[idx10], grid.x_values[i + 1], grid.y_values[j + 1], values[idx11],
                       grid.x_values[i], grid.y_values[j + 1], values[idx01]);
            export_result.face_count += 2;

            vertex_used[idx00] = 1;
            vertex_used[idx10] = 1;
            vertex_used[idx01] = 1;
            vertex_used[idx11] = 1;
        }
    }
    out << "endsolid capdat_stage6_mesh\n";

    for (const uint8_t used : vertex_used) {
        if (used != 0) {
            ++export_result.vertex_count;
        }
    }
    if (!out.good()) {
        throw std::runtime_error("Failed while writing Stage 6 STL mesh: " + path);
    }
    return export_result;
}

Stage6DualMeshExportResult writeStage6ObjMeshesCombined(const Stage4GridDescriptor& grid,
                                                        const std::vector<uint8_t>& obj_vertex_mask,
                                                        const std::vector<double>& outer_values,
                                                        const std::vector<double>& inner_values,
                                                        const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open Stage 6 combined OBJ path for writing: " + path);
    }

    Stage6DualMeshExportResult export_result;
    int global_vertex_index = 1;
    const auto writeSurface = [&](const std::vector<double>& values,
                                  const std::string& object_name,
                                  std::size_t& vertex_count,
                                  std::size_t& face_count) {
        out << "o " << object_name << '\n';
        std::vector<int> vertex_indices(grid.nx * grid.ny, -1);
        for (std::size_t j = 0; j < grid.ny; ++j) {
            for (std::size_t i = 0; i < grid.nx; ++i) {
                const std::size_t idx = stage4NodeIndex(i, j, grid.nx);
                if (obj_vertex_mask[idx] == 0 || !std::isfinite(values[idx])) {
                    continue;
                }
                out << "v " << grid.x_values[i] << ' ' << grid.y_values[j] << ' ' << values[idx] << '\n';
                vertex_indices[idx] = global_vertex_index;
                ++global_vertex_index;
                ++vertex_count;
            }
        }

        for (std::size_t j = 0; (j + 1) < grid.ny; ++j) {
            for (std::size_t i = 0; (i + 1) < grid.nx; ++i) {
                const std::size_t idx00 = stage4NodeIndex(i, j, grid.nx);
                const std::size_t idx10 = stage4NodeIndex(i + 1, j, grid.nx);
                const std::size_t idx01 = stage4NodeIndex(i, j + 1, grid.nx);
                const std::size_t idx11 = stage4NodeIndex(i + 1, j + 1, grid.nx);
                if (obj_vertex_mask[idx00] == 0 || obj_vertex_mask[idx10] == 0 || obj_vertex_mask[idx01] == 0 ||
                    obj_vertex_mask[idx11] == 0 || !std::isfinite(values[idx00]) || !std::isfinite(values[idx10]) ||
                    !std::isfinite(values[idx01]) || !std::isfinite(values[idx11])) {
                    continue;
                }
                const int v00 = vertex_indices[idx00];
                const int v10 = vertex_indices[idx10];
                const int v01 = vertex_indices[idx01];
                const int v11 = vertex_indices[idx11];
                out << "f " << v00 << ' ' << v10 << ' ' << v01 << '\n';
                out << "f " << v10 << ' ' << v11 << ' ' << v01 << '\n';
                face_count += 2;
            }
        }
    };

    writeSurface(outer_values, "outer_surface", export_result.outer_vertex_count, export_result.outer_face_count);
    writeSurface(inner_values, "inner_surface", export_result.inner_vertex_count, export_result.inner_face_count);

    if (!out.good()) {
        throw std::runtime_error("Failed while writing combined Stage 6 OBJ mesh: " + path);
    }
    return export_result;
}

Stage6DualMeshExportResult writeStage6StlMeshesCombined(const Stage4GridDescriptor& grid,
                                                        const std::vector<uint8_t>& obj_vertex_mask,
                                                        const std::vector<double>& outer_values,
                                                        const std::vector<double>& inner_values,
                                                        const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open Stage 6 combined STL path for writing: " + path);
    }

    out << "solid capdat_mesh\n";
    Stage6DualMeshExportResult export_result;

    const auto writeSurface = [&](const std::vector<double>& values, std::size_t& vertex_count, std::size_t& face_count) {
        std::vector<uint8_t> vertex_used(grid.nx * grid.ny, 0);
        for (std::size_t j = 0; (j + 1) < grid.ny; ++j) {
            for (std::size_t i = 0; (i + 1) < grid.nx; ++i) {
                const std::size_t idx00 = stage4NodeIndex(i, j, grid.nx);
                const std::size_t idx10 = stage4NodeIndex(i + 1, j, grid.nx);
                const std::size_t idx01 = stage4NodeIndex(i, j + 1, grid.nx);
                const std::size_t idx11 = stage4NodeIndex(i + 1, j + 1, grid.nx);
                if (obj_vertex_mask[idx00] == 0 || obj_vertex_mask[idx10] == 0 || obj_vertex_mask[idx01] == 0 ||
                    obj_vertex_mask[idx11] == 0 || !std::isfinite(values[idx00]) || !std::isfinite(values[idx10]) ||
                    !std::isfinite(values[idx01]) || !std::isfinite(values[idx11])) {
                    continue;
                }

                const double x00 = grid.x_values[i];
                const double y00 = grid.y_values[j];
                const double x10 = grid.x_values[i + 1];
                const double y10 = grid.y_values[j];
                const double x01 = grid.x_values[i];
                const double y01 = grid.y_values[j + 1];
                const double x11 = grid.x_values[i + 1];
                const double y11 = grid.y_values[j + 1];

                const auto writeFacet = [&](double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2) {
                    const double ux = x1 - x0;
                    const double uy = y1 - y0;
                    const double uz = z1 - z0;
                    const double vx = x2 - x0;
                    const double vy = y2 - y0;
                    const double vz = z2 - z0;
                    double nx = (uy * vz) - (uz * vy);
                    double ny = (uz * vx) - (ux * vz);
                    double nz = (ux * vy) - (uy * vx);
                    const double norm = std::sqrt((nx * nx) + (ny * ny) + (nz * nz));
                    if (norm > 0.0) {
                        nx /= norm;
                        ny /= norm;
                        nz /= norm;
                    }
                    out << "  facet normal " << nx << ' ' << ny << ' ' << nz << '\n'
                        << "    outer loop\n"
                        << "      vertex " << x0 << ' ' << y0 << ' ' << z0 << '\n'
                        << "      vertex " << x1 << ' ' << y1 << ' ' << z1 << '\n'
                        << "      vertex " << x2 << ' ' << y2 << ' ' << z2 << '\n'
                        << "    endloop\n"
                        << "  endfacet\n";
                };

                writeFacet(x00, y00, values[idx00], x10, y10, values[idx10], x01, y01, values[idx01]);
                writeFacet(x10, y10, values[idx10], x11, y11, values[idx11], x01, y01, values[idx01]);
                face_count += 2;

                vertex_used[idx00] = 1;
                vertex_used[idx10] = 1;
                vertex_used[idx01] = 1;
                vertex_used[idx11] = 1;
            }
        }

        for (const uint8_t used : vertex_used) {
            if (used != 0) {
                ++vertex_count;
            }
        }
    };

    writeSurface(outer_values, export_result.outer_vertex_count, export_result.outer_face_count);
    writeSurface(inner_values, export_result.inner_vertex_count, export_result.inner_face_count);
    out << "endsolid capdat_mesh\n";

    if (!out.good()) {
        throw std::runtime_error("Failed while writing combined Stage 6 STL mesh: " + path);
    }
    return export_result;
}

const char* toVdwSourceLabel(VdwResolutionSource source) {
    switch (source) {
    case VdwResolutionSource::explicit_element:
        return "explicit_element";
    case VdwResolutionSource::inferred_from_atom_name:
        return "inferred_from_atom_name";
    case VdwResolutionSource::fallback_unknown:
    default:
        return "fallback_unknown";
    }
}

bool writeStage3NormalizedAtomsCsv(const std::string& path, const GeometryPatchNormalizationResult& stage3_result) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }

    out << "patch_atom_index,serial_number,atom_name,residue_name,chain_id,residue_seq,insertion_code,alt_loc,"
           "is_hetatm,raw_element,assigned_element,vdw_assignment_source,vdw_radius,radial_xy,in_positive_z,"
           "in_cylinder_radius,x,y,z\n";

    for (std::size_t idx = 0; idx < stage3_result.analytical_patch.atoms.size(); ++idx) {
        const PatchAtom& atom = stage3_result.analytical_patch.atoms[idx];
        const Atom* original = atom.original_atom;
        if (original == nullptr) {
            return false;
        }
        const VdwResolution vdw_resolution = resolveVdwElement(*original);
        const char chain_id = original->chainId() == '\0' ? ' ' : original->chainId();
        const char insertion_code = original->insertionCode() == '\0' ? ' ' : original->insertionCode();
        const char alt_loc = original->altLoc() == '\0' ? ' ' : original->altLoc();

        out << idx << ',' << original->serial() << ',' << original->name() << ',' << original->residueName() << ','
            << chain_id << ',' << original->residueSeq() << ',' << insertion_code << ',' << alt_loc << ','
            << (original->isHetatm() ? 1 : 0) << ',' << original->element() << ',' << atom.element << ','
            << toVdwSourceLabel(vdw_resolution.source) << ',' << atom.vdw_radius << ',' << atom.radial_xy << ','
            << (atom.in_positive_z ? 1 : 0) << ',' << (atom.in_cylinder_radius ? 1 : 0) << ',' << atom.position.x
            << ',' << atom.position.y << ',' << atom.position.z << '\n';
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

std::string inferElementSymbolFromAtomName(const std::string& atom_name) {
    std::string letters;
    letters.reserve(atom_name.size());
    for (const char ch : atom_name) {
        if (std::isalpha(static_cast<unsigned char>(ch)) != 0) {
            letters.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
        }
    }
    if (letters.empty()) {
        return "";
    }

    if (letters.size() >= 2 && letters[0] == 'S' && letters[1] == 'E') {
        return "SE";
    }
    return std::string(1, letters[0]);
}

VdwResolution resolveVdwElement(const Atom& atom) {
    const auto& table = vdwRadiusLookupTable();
    const std::string explicit_element = normalizeElementSymbol(atom.element());
    if (!explicit_element.empty() && table.find(explicit_element) != table.end()) {
        return VdwResolution{explicit_element, VdwResolutionSource::explicit_element};
    }

    const std::string inferred_element = inferElementSymbolFromAtomName(atom.name());
    if (!inferred_element.empty() && table.find(inferred_element) != table.end()) {
        return VdwResolution{inferred_element, VdwResolutionSource::inferred_from_atom_name};
    }

    return VdwResolution{"", VdwResolutionSource::fallback_unknown};
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
                        const CylinderMembership& membership,
                        double delta_vdw) {
    const VdwResolution resolved = resolveVdwElement(atom);

    PatchAtom patch_atom;
    patch_atom.position = rotated_position;
    patch_atom.element = resolved.normalized_element;
    patch_atom.vdw_radius = vdwRadius(patch_atom.element) + delta_vdw;
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
    double best_outer_contact_height = -std::numeric_limits<double>::infinity();
    double best_inner_contact_height = std::numeric_limits<double>::infinity();
    std::size_t best_outer_patch_atom_index = 0;
    std::size_t best_inner_patch_atom_index = 0;

    for (std::size_t idx = 0; idx < patch_atoms.size(); ++idx) {
        const LineSphereIntersection intersection =
            intersectVerticalLineWithSphere(x, y, patch_atoms[idx], tie_tolerance);
        if (!intersection.intersects) {
            continue;
        }
        ++result.candidate_patch_atom_count;

        if (!any_intersection || (intersection.z_high > best_outer_contact_height + tie_tolerance)) {
            best_outer_contact_height = intersection.z_high;
            best_outer_patch_atom_index = idx;
        }

        if (!any_intersection || (intersection.z_low < best_inner_contact_height - tie_tolerance)) {
            best_inner_contact_height = intersection.z_low;
            best_inner_patch_atom_index = idx;
        }

        any_intersection = true;
    }

    if (!any_intersection) {
        return result;
    }

    result.z_outer_raw = best_outer_contact_height;
    result.z_inner_raw = best_inner_contact_height;
    result.outer_patch_atom_index = best_outer_patch_atom_index;
    result.inner_patch_atom_index = best_inner_patch_atom_index;
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
    std::string final_working_frame_message = result.final_frame_description;
    if (result.used_identity_rotation && !result.coordinates_modified_in_place &&
        result.final_frame_description == "derived_reoriented_frame") {
        final_working_frame_message += " (already aligned to +Z; workflow recorded without coordinate changes)";
    }
    result.messages.push_back("Geometry Stage 1 final working frame: " + final_working_frame_message);

    if (config.export_rotated_capsid) {
        ExportCapsidConfig writer_config;
        writer_config.output_path = geometryArtifactPath(config, "_rotated_capsid.pdb");
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

                result.patch_atoms.push_back(makePatchAtom(atom, position, membership, config.delta_vdw));
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
    writer_config.output_path = geometryArtifactPath(config, "_patch_atoms.pdb");
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
    const FoldPatchAnalysisConfig& config,
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
    result.messages.push_back("Geometry Stage 3 vdW delta offset: " + std::to_string(config.delta_vdw));
    result.messages.push_back("Geometry Stage 3 selected atoms to normalize: " +
                              std::to_string(stage2_result.patch_atoms.size()));

    result.analytical_patch.atoms.reserve(stage2_result.patch_atoms.size());
    result.analytical_patch.original_atom_refs.reserve(stage2_result.selected_atom_refs.size());
    geometry_symmetry::Vector3 min_position{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
    };
    geometry_symmetry::Vector3 max_position{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };

    for (std::size_t idx = 0; idx < stage2_result.patch_atoms.size(); ++idx) {
        const PatchAtom& selected = stage2_result.patch_atoms[idx];
        const Atom* original_ref = stage2_result.selected_atom_refs[idx];
        if (original_ref == nullptr || selected.original_atom == nullptr || selected.original_atom != original_ref) {
            throw std::runtime_error("PatchAtom normalization encountered a null original atom reference");
        }

        const PatchAtom normalized =
            makePatchAtom(*original_ref, selected.position, selected.membership, config.delta_vdw);
        const VdwResolution vdw_resolution = resolveVdwElement(*original_ref);
        if (vdw_resolution.source == VdwResolutionSource::explicit_element) {
            ++result.analytical_patch.explicit_vdw_radius_count;
        } else if (vdw_resolution.source == VdwResolutionSource::inferred_from_atom_name) {
            ++result.analytical_patch.inferred_vdw_radius_count;
        } else {
            ++result.analytical_patch.fallback_vdw_radius_count;
        }

        result.analytical_patch.atoms.push_back(normalized);
        result.analytical_patch.original_atom_refs.push_back(original_ref);
        min_position.x = std::min(min_position.x, normalized.position.x);
        min_position.y = std::min(min_position.y, normalized.position.y);
        min_position.z = std::min(min_position.z, normalized.position.z);
        max_position.x = std::max(max_position.x, normalized.position.x);
        max_position.y = std::max(max_position.y, normalized.position.y);
        max_position.z = std::max(max_position.z, normalized.position.z);
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
    result.messages.push_back("Geometry Stage 3 inferred-from-name vdW assignments: " +
                              std::to_string(result.analytical_patch.inferred_vdw_radius_count));
    if (result.analytical_patch.fallback_vdw_radius_count == 0) {
        result.messages.push_back("Geometry Stage 3 fallback vdW assignments: " +
                                  std::to_string(result.analytical_patch.fallback_vdw_radius_count));
    } else {
        result.messages.push_back(std::string(kWarningMessagePrefix) +
                                  "Geometry Stage 3 fallback vdW assignments: " +
                                  std::to_string(result.analytical_patch.fallback_vdw_radius_count));
    }
    result.messages.push_back("Geometry Stage 3 patch atom bounds x:[" + std::to_string(min_position.x) + ", " +
                              std::to_string(max_position.x) + "] y:[" + std::to_string(min_position.y) + ", " +
                              std::to_string(max_position.y) + "] z:[" + std::to_string(min_position.z) + ", " +
                              std::to_string(max_position.z) + "]");
    result.messages.push_back("Geometry Stage 3 patch export path (Stage 2 canonical): " +
                              result.analytical_patch.export_path);
    if (result.analytical_patch.fallback_vdw_radius_count > 0) {
        result.messages.push_back(
            std::string(kNoteMessagePrefix) +
            "Geometry Stage 3 note: fallback vdW radii were applied where explicit/inferred element resolution "
            "was unavailable.");
    }
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
    Timer stage4_timer;
    stage4_timer.start();
    result.stage4_start_timestamp_utc = utcTimestampNowIso8601();

    if (!stage3_result.success) {
        throw std::runtime_error("Stage 4 cannot run before successful Stage 3 patch normalization");
    }
    if (stage3_result.analytical_patch.atoms.empty()) {
        throw std::runtime_error("Analytical patch is empty; raw sheet detection cannot proceed");
    }
    validateStage4Config(config);

    result.messages.push_back("Geometry Stage 4");
    result.messages.push_back(
        "Geometry analysis: starting Stage 4 raw outer/inner ray-sphere first-contact envelope detection.");
    result.messages.push_back("Geometry Stage 4 grid spacing: " + std::to_string(config.grid_spacing));
    result.messages.push_back("Geometry Stage 4 cylinder radius: " + std::to_string(config.cylinder_radius));
    if (config.debug) {
        result.stage3_normalized_atoms_csv_path = geometryArtifactPath(config, "_stage3_normalized_atoms.csv");
        if (!writeStage3NormalizedAtomsCsv(result.stage3_normalized_atoms_csv_path, stage3_result)) {
            throw std::runtime_error("Failed to write Stage 3 normalized atoms CSV for Stage 4 input traceability");
        }
        result.messages.push_back("Geometry Stage 4 Stage 3 normalized atom input CSV: " +
                                  result.stage3_normalized_atoms_csv_path);
        result.messages.push_back("Geometry Stage 4 Stage 3 normalized atom rows: " +
                                  std::to_string(stage3_result.analytical_patch.atoms.size()));
    }

    result.grid = buildStage4RegularGrid(config.cylinder_radius, config.grid_spacing, tolerance);
    result.contact_search_patch_atom_count = stage3_result.analytical_patch.atoms.size();
    result.node_count = result.grid.nx * result.grid.ny;
    result.z_outer_raw.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.z_inner_raw.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.inside_disk_mask.assign(result.node_count, 0);
    result.valid_mask.assign(result.node_count, 0);
    result.candidate_patch_atom_counts.assign(result.node_count, 0);
    result.inner_contact_serial_numbers.assign(result.node_count, 0);
    result.outer_contact_serial_numbers.assign(result.node_count, 0);
    result.outer_contact_patch_atom_indices.assign(result.node_count, -1);
    result.inner_contact_patch_atom_indices.assign(result.node_count, -1);
    result.atom_roles.assign(stage3_result.analytical_patch.atoms.size(), PatchAtomContactRole::none);
    std::vector<uint8_t> outer_only_mask(result.node_count, 0);
    std::vector<uint8_t> inner_only_mask(result.node_count, 0);
    std::vector<uint8_t> negative_thickness_mask(result.node_count, 0);

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
                result.candidate_patch_atom_counts[idx] = node.candidate_patch_atom_count;
                continue;
            }

            result.valid_mask[idx] = 1;
            result.candidate_patch_atom_counts[idx] = node.candidate_patch_atom_count;
            result.z_outer_raw[idx] = node.z_outer_raw;
            result.z_inner_raw[idx] = node.z_inner_raw;
            const Atom* outer_atom = stage3_result.analytical_patch.atoms[node.outer_patch_atom_index].original_atom;
            const Atom* inner_atom = stage3_result.analytical_patch.atoms[node.inner_patch_atom_index].original_atom;
            if (outer_atom == nullptr || inner_atom == nullptr) {
                throw std::runtime_error("Stage 4 encountered null original atom while writing node contacts");
            }
            result.outer_contact_serial_numbers[idx] = outer_atom->serial();
            result.inner_contact_serial_numbers[idx] = inner_atom->serial();
            result.outer_contact_patch_atom_indices[idx] = static_cast<int>(node.outer_patch_atom_index);
            result.inner_contact_patch_atom_indices[idx] = static_cast<int>(node.inner_patch_atom_index);
            ++result.valid_node_count;
            ++result.both_hit_node_count;

            const double thickness = node.z_outer_raw - node.z_inner_raw;
            if (std::fabs(thickness) <= tolerance) {
                ++result.zero_thickness_node_count;
            }
            if (thickness < -tolerance) {
                ++result.negative_thickness_node_count;
                negative_thickness_mask[idx] = 1;
            }

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
    geometry_symmetry::Vector3 contact_min_position{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    geometry_symmetry::Vector3 contact_max_position{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()};
    for (std::size_t idx = 0; idx < stage3_result.analytical_patch.atoms.size(); ++idx) {
        result.atom_roles[idx] = classifyContactRole(used_as_outer[idx], used_as_inner[idx]);
        if (used_as_outer[idx] && used_as_inner[idx]) {
            ++result.unique_both_contact_atom_count;
        }
        if (used_as_outer[idx]) {
            ++result.unique_outer_contact_atom_count;
        }
        if (used_as_inner[idx]) {
            ++result.unique_inner_contact_atom_count;
        }
        if (used_as_outer[idx] || used_as_inner[idx]) {
            ++result.unique_contact_atom_count;
            const geometry_symmetry::Vector3& position = stage3_result.analytical_patch.atoms[idx].position;
            contact_min_position.x = std::min(contact_min_position.x, position.x);
            contact_min_position.y = std::min(contact_min_position.y, position.y);
            contact_min_position.z = std::min(contact_min_position.z, position.z);
            contact_max_position.x = std::max(contact_max_position.x, position.x);
            contact_max_position.y = std::max(contact_max_position.y, position.y);
            contact_max_position.z = std::max(contact_max_position.z, position.z);
            contact_atom_subset.push_back(stage3_result.analytical_patch.atoms[idx].original_atom);
        }
    }

    for (std::size_t idx = 0; idx < stage3_result.analytical_patch.atoms.size(); ++idx) {
        const bool is_outer = used_as_outer[idx];
        const bool is_inner = used_as_inner[idx];
        if (!is_outer && !is_inner) {
            continue;
        }
        const Atom* original_atom = stage3_result.analytical_patch.atoms[idx].original_atom;
        if (original_atom == nullptr) {
            throw std::runtime_error("Stage 4 encountered null original atom while aggregating serial provenance");
        }
        const int serial = original_atom->serial();
        if (is_outer) {
            result.unique_outer_contact_serials.push_back(serial);
        }
        if (is_inner) {
            result.unique_inner_contact_serials.push_back(serial);
        }
        if (is_outer && is_inner) {
            result.unique_both_contact_serials.push_back(serial);
        }
        result.unique_contact_serial_union.push_back(serial);
    }
    std::sort(result.unique_outer_contact_serials.begin(), result.unique_outer_contact_serials.end());
    result.unique_outer_contact_serials.erase(
        std::unique(result.unique_outer_contact_serials.begin(), result.unique_outer_contact_serials.end()),
        result.unique_outer_contact_serials.end());
    std::sort(result.unique_inner_contact_serials.begin(), result.unique_inner_contact_serials.end());
    result.unique_inner_contact_serials.erase(
        std::unique(result.unique_inner_contact_serials.begin(), result.unique_inner_contact_serials.end()),
        result.unique_inner_contact_serials.end());
    std::sort(result.unique_both_contact_serials.begin(), result.unique_both_contact_serials.end());
    result.unique_both_contact_serials.erase(
        std::unique(result.unique_both_contact_serials.begin(), result.unique_both_contact_serials.end()),
        result.unique_both_contact_serials.end());
    std::sort(result.unique_contact_serial_union.begin(), result.unique_contact_serial_union.end());
    result.unique_contact_serial_union.erase(
        std::unique(result.unique_contact_serial_union.begin(), result.unique_contact_serial_union.end()),
        result.unique_contact_serial_union.end());
    result.unique_outer_contact_serial_count = result.unique_outer_contact_serials.size();
    result.unique_inner_contact_serial_count = result.unique_inner_contact_serials.size();
    result.unique_both_contact_serial_count = result.unique_both_contact_serials.size();
    result.unique_contact_serial_union_count = result.unique_contact_serial_union.size();

    if (config.debug) {
        result.outer_csv_path = geometryArtifactPath(config, "_outer_raw.csv");
        result.inner_csv_path = geometryArtifactPath(config, "_inner_raw.csv");
        result.valid_mask_csv_path = geometryArtifactPath(config, "_valid_mask_raw.csv");
        result.outer_only_mask_csv_path = geometryArtifactPath(config, "_outer_only_mask_raw.csv");
        result.inner_only_mask_csv_path = geometryArtifactPath(config, "_inner_only_mask_raw.csv");
        result.negative_thickness_mask_csv_path = geometryArtifactPath(config, "_negative_thickness_mask_raw.csv");
        result.summary_csv_path = geometryArtifactPath(config, "_stage4_summary.csv");
    }
    result.contact_atoms_pdb_path = geometryArtifactPath(config, "_raw_contact_atoms.pdb");

    if (config.debug) {
        if (!writeStage4CsvOuter(result)) {
            throw std::runtime_error("Failed to write raw outer grid CSV");
        }
        if (!writeStage4CsvInner(result)) {
            throw std::runtime_error("Failed to write raw inner grid CSV");
        }
        if (!writeStage4CsvValidMask(result)) {
            throw std::runtime_error("Failed to write raw valid-mask grid CSV");
        }
        if (!writeStage4CsvBinaryMask(result, result.outer_only_mask_csv_path, "outer_only", outer_only_mask)) {
            throw std::runtime_error("Failed to write raw outer-only mask CSV");
        }
        if (!writeStage4CsvBinaryMask(result, result.inner_only_mask_csv_path, "inner_only", inner_only_mask)) {
            throw std::runtime_error("Failed to write raw inner-only mask CSV");
        }
        if (!writeStage4CsvBinaryMask(
                result, result.negative_thickness_mask_csv_path, "negative_thickness", negative_thickness_mask)) {
            throw std::runtime_error("Failed to write raw negative-thickness mask CSV");
        }
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

    stage4_timer.stop();
    result.stage4_end_timestamp_utc = utcTimestampNowIso8601();
    result.stage4_runtime_seconds = stage4_timer.elapsedSeconds();
    if (config.debug && !writeStage4SummaryCsv(result)) {
        throw std::runtime_error("Failed to write Stage 4 summary CSV");
    }

    result.messages.push_back("Geometry Stage 4 grid dimensions: nx=" + std::to_string(result.grid.nx) +
                              ", ny=" + std::to_string(result.grid.ny));
    result.messages.push_back("Geometry Stage 4 total nodes: " + std::to_string(result.node_count));
    result.messages.push_back("Geometry Stage 4 nodes inside disk: " + std::to_string(result.inside_disk_count));
    result.messages.push_back("Geometry Stage 4 valid nodes: " + std::to_string(result.valid_node_count));
    result.messages.push_back("Geometry Stage 4 invalid nodes: " + std::to_string(result.invalid_node_count));
    result.messages.push_back("Geometry Stage 4 outer-only nodes: " + std::to_string(result.outer_only_node_count));
    result.messages.push_back("Geometry Stage 4 inner-only nodes: " + std::to_string(result.inner_only_node_count));
    result.messages.push_back("Geometry Stage 4 both-hit nodes: " + std::to_string(result.both_hit_node_count));
    if (result.zero_thickness_node_count == 0) {
        result.messages.push_back("Geometry Stage 4 zero-thickness nodes: " +
                                  std::to_string(result.zero_thickness_node_count));
    } else {
        result.messages.push_back(std::string(kWarningMessagePrefix) +
                                  "Geometry Stage 4 zero-thickness nodes: " +
                                  std::to_string(result.zero_thickness_node_count));
    }
    if (result.negative_thickness_node_count == 0) {
        result.messages.push_back("Geometry Stage 4 negative-thickness nodes: " +
                                  std::to_string(result.negative_thickness_node_count));
    } else {
        result.messages.push_back(std::string(kWarningMessagePrefix) +
                                  "Geometry Stage 4 negative-thickness nodes: " +
                                  std::to_string(result.negative_thickness_node_count));
    }
    result.messages.push_back("Geometry Stage 4 patch atoms used for contact search: " +
                              std::to_string(result.contact_search_patch_atom_count));
    result.messages.push_back("Geometry Stage 4 unique outer-contact patch atoms (index-based): " +
                              std::to_string(result.unique_outer_contact_atom_count));
    result.messages.push_back("Geometry Stage 4 unique inner-contact patch atoms (index-based): " +
                              std::to_string(result.unique_inner_contact_atom_count));
    result.messages.push_back("Geometry Stage 4 unique both-contact patch atoms (index-based): " +
                              std::to_string(result.unique_both_contact_atom_count));
    result.messages.push_back("Geometry Stage 4 unique contact patch atoms union (index-based): " +
                              std::to_string(result.unique_contact_atom_count));
    result.messages.push_back("Geometry Stage 4 unique outer-contact atom serials: " +
                              std::to_string(result.unique_outer_contact_serial_count));
    result.messages.push_back("Geometry Stage 4 unique inner-contact atom serials: " +
                              std::to_string(result.unique_inner_contact_serial_count));
    result.messages.push_back("Geometry Stage 4 unique both-contact atom serials: " +
                              std::to_string(result.unique_both_contact_serial_count));
    result.messages.push_back("Geometry Stage 4 unique contact atom serial union: " +
                              std::to_string(result.unique_contact_serial_union_count));
    result.messages.push_back("Geometry Stage 4 unique contact atom bounds x:[" +
                              std::to_string(contact_min_position.x) + ", " +
                              std::to_string(contact_max_position.x) + "] y:[" +
                              std::to_string(contact_min_position.y) + ", " +
                              std::to_string(contact_max_position.y) + "] z:[" +
                              std::to_string(contact_min_position.z) + ", " +
                              std::to_string(contact_max_position.z) + "]");
    if (config.debug) {
        result.messages.push_back("Geometry Stage 4 outer CSV: " + result.outer_csv_path);
        result.messages.push_back("Geometry Stage 4 inner CSV: " + result.inner_csv_path);
        result.messages.push_back("Geometry Stage 4 valid mask CSV: " + result.valid_mask_csv_path);
        result.messages.push_back("Geometry Stage 4 outer-only mask CSV: " + result.outer_only_mask_csv_path);
        result.messages.push_back("Geometry Stage 4 inner-only mask CSV: " + result.inner_only_mask_csv_path);
        result.messages.push_back("Geometry Stage 4 negative-thickness mask CSV: " +
                                  result.negative_thickness_mask_csv_path);
        result.messages.push_back("Geometry Stage 4 summary CSV: " + result.summary_csv_path);
    }
    result.messages.push_back("Geometry Stage 4 contact atoms PDB: " + result.contact_atoms_pdb_path);
    result.messages.push_back("Geometry Stage 4 start timestamp (UTC): " + result.stage4_start_timestamp_utc);
    result.messages.push_back("Geometry Stage 4 end timestamp (UTC): " + result.stage4_end_timestamp_utc);
    result.messages.push_back("Geometry Stage 4 runtime seconds: " + std::to_string(result.stage4_runtime_seconds));
    if (config.debug &&
        (result.zero_thickness_node_count > 0 || result.negative_thickness_node_count > 0 ||
         result.unique_both_contact_atom_count > 0)) {
        result.messages.push_back(
            std::string(kNoteMessagePrefix) +
            "Geometry Stage 4 note: non-zero zero-thickness/negative-thickness nodes or both-contact atoms "
            "indicate potential raw-sheet overlap/thickness anomalies.");
    }
    result.messages.push_back("Geometry analysis: completed Stage 4 raw outer/inner sheet detection.");

    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryStage5SurfacePrepResult runGeometryAnalysisStage5SurfacePreparation(
    const GeometryStage4RawSheetResult& stage4_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance) {
    GeometryStage5SurfacePrepResult result;
    if (!stage4_result.success) {
        throw std::runtime_error("Stage 5 cannot run before successful Stage 4 raw sheet detection");
    }
    if (stage4_result.grid.nx == 0 || stage4_result.grid.ny == 0) {
        throw std::runtime_error("Stage 5 requires a non-empty Stage 4 grid");
    }
    if (stage4_result.node_count == 0) {
        throw std::runtime_error("Stage 5 requires Stage 4 node_count > 0");
    }
    if (stage4_result.inside_disk_count == 0) {
        throw std::runtime_error("Stage 5 requires Stage 4 inside-disk node_count > 0");
    }
    if (stage4_result.valid_node_count == 0) {
        throw std::runtime_error("Stage 5 requires Stage 4 valid_node_count > 0");
    }

    result.messages.push_back("Geometry Stage 5");
    result.messages.push_back("Geometry analysis: starting Stage 5 surface preparation.");

    result.grid = stage4_result.grid;
    result.node_count = stage4_result.node_count;
    result.inside_disk_mask = stage4_result.inside_disk_mask;
    result.raw_valid_mask = stage4_result.valid_mask;
    result.outer_contact_serial_numbers = stage4_result.outer_contact_serial_numbers;
    result.inner_contact_serial_numbers = stage4_result.inner_contact_serial_numbers;
    result.outer_contact_patch_atom_indices = stage4_result.outer_contact_patch_atom_indices;
    result.inner_contact_patch_atom_indices = stage4_result.inner_contact_patch_atom_indices;
    result.inside_disk_count = stage4_result.inside_disk_count;
    result.raw_valid_node_count = stage4_result.valid_node_count;
    result.raw_invalid_node_count = result.inside_disk_count - result.raw_valid_node_count;

    result.z_outer_seed.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.z_inner_seed.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.outer_seed_mask.assign(result.node_count, 0);
    result.inner_seed_mask.assign(result.node_count, 0);
    result.paired_seed_mask.assign(result.node_count, 0);
    result.boundary_exclusion_mask.assign(result.node_count, 0);
    result.interp_allowed_outer_mask.assign(result.node_count, 0);
    result.interp_allowed_inner_mask.assign(result.node_count, 0);
    result.paired_interp_allowed_mask.assign(result.node_count, 0);
    result.hard_invalid_mask.assign(result.node_count, 0);
    result.reliable_core_mask.assign(result.node_count, 0);

    const double derived_boundary_margin =
        config.stage5_boundary_margin > 0.0 ? config.stage5_boundary_margin : (2.0 * result.grid.spacing);
    const double derived_support_radius =
        config.stage5_support_radius > 0.0 ? config.stage5_support_radius : (2.5 * result.grid.spacing);
    const double derived_reliable_radius =
        config.stage5_reliable_radius > 0.0
            ? config.stage5_reliable_radius
            : std::max(0.0, config.cylinder_radius - derived_boundary_margin);

    if (derived_boundary_margin < 0.0) {
        throw std::runtime_error("Stage 5 requires boundary_margin >= 0");
    }
    if (derived_support_radius <= 0.0) {
        throw std::runtime_error("Stage 5 requires support_radius > 0");
    }
    if (derived_reliable_radius < 0.0) {
        throw std::runtime_error("Stage 5 requires reliable_radius >= 0");
    }
    if (config.stage5_min_support_nodes == 0) {
        throw std::runtime_error("Stage 5 requires stage5_min_support_nodes > 0");
    }

    result.boundary_margin = derived_boundary_margin;
    result.support_radius = derived_support_radius;
    result.reliable_radius = std::min(config.cylinder_radius, derived_reliable_radius);

    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            if (result.inside_disk_mask[idx] == 0) {
                continue;
            }
            if (result.raw_valid_mask[idx] == 0) {
                continue;
            }
            result.z_outer_seed[idx] = stage4_result.z_outer_raw[idx];
            result.z_inner_seed[idx] = stage4_result.z_inner_raw[idx];
            result.outer_seed_mask[idx] = 1;
            result.inner_seed_mask[idx] = 1;
            result.paired_seed_mask[idx] = 1;
            ++result.outer_seed_node_count;
            ++result.inner_seed_node_count;
            ++result.paired_seed_node_count;
        }
    }

    const double boundary_limit = config.cylinder_radius - result.boundary_margin;
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            if (result.inside_disk_mask[idx] == 0) {
                continue;
            }
            const double radial = stage5NodeRadius(result.grid, i, j);
            if (radial > boundary_limit + tolerance) {
                result.boundary_exclusion_mask[idx] = 1;
                ++result.boundary_excluded_node_count;
            }
        }
    }

    // Stage 5 seed-atom traceability is defined over the full paired-seed domain.
    // Do not filter by boundary/reliability/interpolation masks here: every
    // Stage 4-valid node promoted to paired_seed contributes provenance serials.
    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.paired_seed_mask[idx] == 0) {
            continue;
        }
        const int outer_patch_index = result.outer_contact_patch_atom_indices[idx];
        const int inner_patch_index = result.inner_contact_patch_atom_indices[idx];
        if (outer_patch_index >= 0) {
            result.unique_outer_seed_patch_atom_indices.push_back(static_cast<std::size_t>(outer_patch_index));
        }
        if (inner_patch_index >= 0) {
            result.unique_inner_seed_patch_atom_indices.push_back(static_cast<std::size_t>(inner_patch_index));
        }

        if (outer_patch_index >= 0 || inner_patch_index >= 0) {
            if (outer_patch_index >= 0) {
                result.unique_seed_patch_atom_index_union.push_back(static_cast<std::size_t>(outer_patch_index));
            }
            if (inner_patch_index >= 0) {
                result.unique_seed_patch_atom_index_union.push_back(static_cast<std::size_t>(inner_patch_index));
            }
        } else {
            // Fallback path for synthetic tests or legacy Stage 4 payloads that
            // do not provide per-node patch-atom indices.
            if (result.outer_contact_serial_numbers[idx] != 0 || result.inner_contact_serial_numbers[idx] != 0) {
                result.unique_seed_patch_atom_index_union.push_back(idx);
            }
        }

        result.unique_outer_seed_atom_serials.push_back(result.outer_contact_serial_numbers[idx]);
        result.unique_inner_seed_atom_serials.push_back(result.inner_contact_serial_numbers[idx]);
    }
    std::sort(result.unique_outer_seed_patch_atom_indices.begin(), result.unique_outer_seed_patch_atom_indices.end());
    result.unique_outer_seed_patch_atom_indices.erase(std::unique(result.unique_outer_seed_patch_atom_indices.begin(),
                                                                  result.unique_outer_seed_patch_atom_indices.end()),
                                                      result.unique_outer_seed_patch_atom_indices.end());
    std::sort(result.unique_inner_seed_patch_atom_indices.begin(), result.unique_inner_seed_patch_atom_indices.end());
    result.unique_inner_seed_patch_atom_indices.erase(std::unique(result.unique_inner_seed_patch_atom_indices.begin(),
                                                                  result.unique_inner_seed_patch_atom_indices.end()),
                                                      result.unique_inner_seed_patch_atom_indices.end());
    std::sort(result.unique_seed_patch_atom_index_union.begin(), result.unique_seed_patch_atom_index_union.end());
    result.unique_seed_patch_atom_index_union.erase(std::unique(result.unique_seed_patch_atom_index_union.begin(),
                                                                result.unique_seed_patch_atom_index_union.end()),
                                                    result.unique_seed_patch_atom_index_union.end());
    std::sort(result.unique_outer_seed_atom_serials.begin(), result.unique_outer_seed_atom_serials.end());
    result.unique_outer_seed_atom_serials.erase(
        std::unique(result.unique_outer_seed_atom_serials.begin(), result.unique_outer_seed_atom_serials.end()),
        result.unique_outer_seed_atom_serials.end());
    std::sort(result.unique_inner_seed_atom_serials.begin(), result.unique_inner_seed_atom_serials.end());
    result.unique_inner_seed_atom_serials.erase(
        std::unique(result.unique_inner_seed_atom_serials.begin(), result.unique_inner_seed_atom_serials.end()),
        result.unique_inner_seed_atom_serials.end());
    result.unique_outer_seed_patch_atom_count = result.unique_outer_seed_patch_atom_indices.size();
    result.unique_inner_seed_patch_atom_count = result.unique_inner_seed_patch_atom_indices.size();
    result.unique_seed_patch_atom_index_union_count = result.unique_seed_patch_atom_index_union.size();

    const double support_radius2 = result.support_radius * result.support_radius;
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            if (result.inside_disk_mask[idx] == 0 || result.paired_seed_mask[idx] != 0) {
                continue;
            }
            if (result.boundary_exclusion_mask[idx] != 0) {
                result.hard_invalid_mask[idx] = 1;
                ++result.hard_invalid_node_count;
                continue;
            }

            std::size_t paired_neighbors = 0;
            std::size_t outer_neighbors = 0;
            std::size_t inner_neighbors = 0;
            const double x0 = result.grid.x_values[i];
            const double y0 = result.grid.y_values[j];
            for (std::size_t nj = 0; nj < result.grid.ny; ++nj) {
                const double dy = result.grid.y_values[nj] - y0;
                for (std::size_t ni = 0; ni < result.grid.nx; ++ni) {
                    const double dx = result.grid.x_values[ni] - x0;
                    const double d2 = (dx * dx) + (dy * dy);
                    if (d2 > support_radius2 + tolerance) {
                        continue;
                    }
                    const std::size_t nidx = stage4NodeIndex(ni, nj, result.grid.nx);
                    if (result.paired_seed_mask[nidx] != 0) {
                        ++paired_neighbors;
                    }
                    if (result.outer_seed_mask[nidx] != 0) {
                        ++outer_neighbors;
                    }
                    if (result.inner_seed_mask[nidx] != 0) {
                        ++inner_neighbors;
                    }
                }
            }

            const bool allowed = paired_neighbors >= config.stage5_min_support_nodes &&
                                 outer_neighbors >= config.stage5_min_support_nodes &&
                                 inner_neighbors >= config.stage5_min_support_nodes;
            if (allowed) {
                result.interp_allowed_outer_mask[idx] = 1;
                result.interp_allowed_inner_mask[idx] = 1;
                result.paired_interp_allowed_mask[idx] = 1;
                ++result.interp_allowed_outer_node_count;
                ++result.interp_allowed_inner_node_count;
                ++result.paired_interp_allowed_node_count;
            } else {
                result.hard_invalid_mask[idx] = 1;
                ++result.hard_invalid_node_count;
            }
        }
    }

    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            if (result.inside_disk_mask[idx] == 0 || result.hard_invalid_mask[idx] != 0) {
                continue;
            }
            const double radial = stage5NodeRadius(result.grid, i, j);
            if (radial <= result.reliable_radius + tolerance) {
                result.reliable_core_mask[idx] = 1;
                ++result.reliable_core_node_count;
            }
        }
    }

    if (config.debug) {
        result.outer_seed_csv_path = geometryArtifactPath(config, "_outer_seed.csv");
        result.inner_seed_csv_path = geometryArtifactPath(config, "_inner_seed.csv");
        result.paired_seed_mask_csv_path = geometryArtifactPath(config, "_paired_seed_mask.csv");
        result.boundary_exclusion_mask_csv_path = geometryArtifactPath(config, "_boundary_exclusion_mask.csv");
        result.interp_allowed_mask_csv_path = geometryArtifactPath(config, "_interp_allowed_mask.csv");
        result.hard_invalid_mask_csv_path = geometryArtifactPath(config, "_hard_invalid_mask.csv");
        result.reliable_core_mask_csv_path = geometryArtifactPath(config, "_reliable_core_mask.csv");
        result.summary_csv_path = geometryArtifactPath(config, "_stage5_summary.csv");

        if (!writeStage5SeedCsv(result, result.outer_seed_csv_path, "z_outer_seed", result.z_outer_seed)) {
            throw std::runtime_error("Failed to write Stage 5 outer seed CSV");
        }
        if (!writeStage5SeedCsv(result, result.inner_seed_csv_path, "z_inner_seed", result.z_inner_seed)) {
            throw std::runtime_error("Failed to write Stage 5 inner seed CSV");
        }
        if (!writeStage5MaskCsv(result, result.paired_seed_mask_csv_path, "paired_seed", result.paired_seed_mask)) {
            throw std::runtime_error("Failed to write Stage 5 paired-seed mask CSV");
        }
        if (!writeStage5MaskCsv(
                result, result.boundary_exclusion_mask_csv_path, "boundary_excluded", result.boundary_exclusion_mask)) {
            throw std::runtime_error("Failed to write Stage 5 boundary-exclusion mask CSV");
        }
        if (!writeStage5MaskCsv(
                result, result.interp_allowed_mask_csv_path, "paired_interp_allowed", result.paired_interp_allowed_mask)) {
            throw std::runtime_error("Failed to write Stage 5 interpolation-allowed mask CSV");
        }
        if (!writeStage5MaskCsv(result, result.hard_invalid_mask_csv_path, "hard_invalid", result.hard_invalid_mask)) {
            throw std::runtime_error("Failed to write Stage 5 hard-invalid mask CSV");
        }
        if (!writeStage5MaskCsv(result, result.reliable_core_mask_csv_path, "reliable_core", result.reliable_core_mask)) {
            throw std::runtime_error("Failed to write Stage 5 reliable-core mask CSV");
        }
        if (!writeStage5SummaryCsv(result)) {
            throw std::runtime_error("Failed to write Stage 5 summary CSV");
        }
    }

    result.messages.push_back("Geometry Stage 5 boundary margin: " + std::to_string(result.boundary_margin));
    result.messages.push_back("Geometry Stage 5 support radius: " + std::to_string(result.support_radius));
    result.messages.push_back("Geometry Stage 5 reliable radius: " + std::to_string(result.reliable_radius));
    result.messages.push_back("Geometry Stage 5 raw valid nodes: " + std::to_string(result.raw_valid_node_count));
    result.messages.push_back("Geometry Stage 5 paired seed nodes: " + std::to_string(result.paired_seed_node_count));
    result.messages.push_back("Geometry Stage 5 boundary-excluded nodes: " +
                              std::to_string(result.boundary_excluded_node_count));
    result.messages.push_back("Geometry Stage 5 interpolation-allowed nodes: " +
                              std::to_string(result.paired_interp_allowed_node_count));
    result.messages.push_back("Geometry Stage 5 hard-invalid nodes: " + std::to_string(result.hard_invalid_node_count));
    result.messages.push_back("Geometry Stage 5 reliable core nodes: " +
                              std::to_string(result.reliable_core_node_count));
    result.messages.push_back("Geometry Stage 5 unique outer seed patch atoms (index-based): " +
                              std::to_string(result.unique_outer_seed_patch_atom_count));
    result.messages.push_back("Geometry Stage 5 unique inner seed patch atoms (index-based): " +
                              std::to_string(result.unique_inner_seed_patch_atom_count));
    result.messages.push_back("Geometry Stage 5 unique seed patch atom union (index-based): " +
                              std::to_string(result.unique_seed_patch_atom_index_union_count));
    result.messages.push_back("Geometry Stage 5 unique outer seed atom serials: " +
                              std::to_string(result.unique_outer_seed_atom_serials.size()));
    result.messages.push_back("Geometry Stage 5 unique inner seed atom serials: " +
                              std::to_string(result.unique_inner_seed_atom_serials.size()));
    if (config.debug) {
        result.messages.push_back("Geometry Stage 5 outer seed CSV: " + result.outer_seed_csv_path);
        result.messages.push_back("Geometry Stage 5 inner seed CSV: " + result.inner_seed_csv_path);
        result.messages.push_back("Geometry Stage 5 paired-seed mask CSV: " + result.paired_seed_mask_csv_path);
        result.messages.push_back("Geometry Stage 5 boundary-exclusion mask CSV: " +
                                  result.boundary_exclusion_mask_csv_path);
        result.messages.push_back("Geometry Stage 5 interpolation-allowed mask CSV: " +
                                  result.interp_allowed_mask_csv_path);
        result.messages.push_back("Geometry Stage 5 hard-invalid mask CSV: " + result.hard_invalid_mask_csv_path);
        result.messages.push_back("Geometry Stage 5 reliable-core mask CSV: " + result.reliable_core_mask_csv_path);
        result.messages.push_back("Geometry Stage 5 summary CSV: " + result.summary_csv_path);
    }
    result.messages.push_back("Geometry analysis: completed Stage 5 surface preparation.");

    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryStage6SurfaceReconstructionResult runGeometryAnalysisStage6SurfaceReconstruction(
    const GeometryStage5SurfacePrepResult& stage5_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance) {
    (void)tolerance;
    GeometryStage6SurfaceReconstructionResult result;
    if (!stage5_result.success) {
        throw std::runtime_error("Stage 6 cannot run before successful Stage 5 surface preparation");
    }
    if (stage5_result.grid.nx == 0 || stage5_result.grid.ny == 0 || stage5_result.node_count == 0) {
        throw std::runtime_error("Stage 6 requires a non-empty Stage 5 grid");
    }
    if (stage5_result.paired_seed_node_count == 0) {
        throw std::runtime_error("Stage 6 requires Stage 5 paired_seed_node_count > 0");
    }
    if (stage5_result.paired_seed_node_count + stage5_result.paired_interp_allowed_node_count == 0) {
        throw std::runtime_error("Stage 6 requires Stage 5 seeds and/or interpolation-allowed nodes");
    }
    if (config.stage6_max_iterations == 0) {
        throw std::runtime_error("Stage 6 requires stage6_max_iterations > 0");
    }
    if (config.stage6_convergence_tolerance <= 0.0) {
        throw std::runtime_error("Stage 6 requires stage6_convergence_tolerance > 0");
    }
    if (config.stage6_smoothing_weight <= 0.0) {
        throw std::runtime_error("Stage 6 requires stage6_smoothing_weight > 0");
    }

    result.messages.push_back("Geometry Stage 6");
    result.messages.push_back("Geometry analysis: starting Stage 6 surface reconstruction.");

    result.grid = stage5_result.grid;
    result.node_count = stage5_result.node_count;
    result.seed_node_count = stage5_result.paired_seed_node_count;
    result.interp_node_count = stage5_result.paired_interp_allowed_node_count;
    result.inside_disk_mask = stage5_result.inside_disk_mask;
    result.paired_seed_mask = stage5_result.paired_seed_mask;
    result.paired_interp_allowed_mask = stage5_result.paired_interp_allowed_mask;
    result.hard_invalid_mask = stage5_result.hard_invalid_mask;
    result.z_outer_reconstructed.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.z_inner_reconstructed.assign(result.node_count, std::numeric_limits<double>::quiet_NaN());
    result.reconstructed_mask.assign(result.node_count, 0);
    result.final_valid_analysis_mask.assign(result.node_count, 0);
    result.non_crossing_adjustment_mask.assign(result.node_count, 0);
    result.obj_vertex_mask.assign(result.node_count, 0);

    const Stage6FieldSolveResult outer = runStage6FieldReconstruction(result.grid,
                                                                       stage5_result.z_outer_seed,
                                                                       stage5_result.inside_disk_mask,
                                                                       stage5_result.paired_seed_mask,
                                                                       stage5_result.paired_interp_allowed_mask,
                                                                       stage5_result.hard_invalid_mask,
                                                                       config);
    const Stage6FieldSolveResult inner = runStage6FieldReconstruction(result.grid,
                                                                       stage5_result.z_inner_seed,
                                                                       stage5_result.inside_disk_mask,
                                                                       stage5_result.paired_seed_mask,
                                                                       stage5_result.paired_interp_allowed_mask,
                                                                       stage5_result.hard_invalid_mask,
                                                                       config);
    result.z_outer_reconstructed = outer.values;
    result.z_inner_reconstructed = inner.values;
    result.outer_iterations_used = outer.iterations_used;
    result.inner_iterations_used = inner.iterations_used;
    result.outer_final_max_update = outer.final_max_update;
    result.inner_final_max_update = inner.final_max_update;

    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        const bool is_seed = result.paired_seed_mask[idx] != 0;
        const bool is_interp = result.paired_interp_allowed_mask[idx] != 0;
        if (!is_seed && !is_interp) {
            result.z_outer_reconstructed[idx] = std::numeric_limits<double>::quiet_NaN();
            result.z_inner_reconstructed[idx] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        if (result.inside_disk_mask[idx] == 0 || result.hard_invalid_mask[idx] != 0) {
            result.z_outer_reconstructed[idx] = std::numeric_limits<double>::quiet_NaN();
            result.z_inner_reconstructed[idx] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        if (config.stage6_enforce_non_crossing && std::isfinite(result.z_outer_reconstructed[idx]) &&
            std::isfinite(result.z_inner_reconstructed[idx]) &&
            result.z_outer_reconstructed[idx] < (result.z_inner_reconstructed[idx] + config.stage6_min_separation)) {
            const double mid = 0.5 * (result.z_outer_reconstructed[idx] + result.z_inner_reconstructed[idx]);
            result.z_inner_reconstructed[idx] = mid - (0.5 * config.stage6_min_separation);
            result.z_outer_reconstructed[idx] = mid + (0.5 * config.stage6_min_separation);
            result.non_crossing_adjustment_mask[idx] = 1;
            ++result.non_crossing_adjusted_node_count;
        }
    }

    double separation_sum = 0.0;
    std::size_t separation_count = 0;
    bool first_separation = true;
    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        const bool in_domain = result.paired_seed_mask[idx] != 0 || result.paired_interp_allowed_mask[idx] != 0;
        const bool finite_pair =
            std::isfinite(result.z_outer_reconstructed[idx]) && std::isfinite(result.z_inner_reconstructed[idx]);
        if (in_domain && finite_pair) {
            result.reconstructed_mask[idx] = 1;
            ++result.reconstructed_node_count;
            const double separation = result.z_outer_reconstructed[idx] - result.z_inner_reconstructed[idx];
            separation_sum += separation;
            if (first_separation) {
                result.min_reconstructed_separation = separation;
                result.max_reconstructed_separation = separation;
                first_separation = false;
            } else {
                result.min_reconstructed_separation = std::min(result.min_reconstructed_separation, separation);
                result.max_reconstructed_separation = std::max(result.max_reconstructed_separation, separation);
            }
            ++separation_count;
        }
        if (result.reconstructed_mask[idx] != 0 && result.inside_disk_mask[idx] != 0 && result.hard_invalid_mask[idx] == 0) {
            result.final_valid_analysis_mask[idx] = 1;
            result.obj_vertex_mask[idx] = 1;
            ++result.final_valid_analysis_node_count;
        }
    }
    result.unresolved_node_count = result.interp_node_count > 0 ? (result.interp_node_count - std::min(result.interp_node_count,
                                                                                                       result.reconstructed_node_count > result.seed_node_count
                                                                                                           ? result.reconstructed_node_count - result.seed_node_count
                                                                                                           : static_cast<std::size_t>(0)))
                                                                 : 0;

    if (result.reconstructed_node_count == 0 || separation_count == 0) {
        throw std::runtime_error("Stage 6 produced zero reconstructed nodes");
    }
    result.mean_reconstructed_separation = separation_sum / static_cast<double>(separation_count);

    if (config.debug) {
        result.outer_reconstructed_csv_path = geometryArtifactPath(config, "_outer_reconstructed.csv");
        result.inner_reconstructed_csv_path = geometryArtifactPath(config, "_inner_reconstructed.csv");
        result.reconstructed_mask_csv_path = geometryArtifactPath(config, "_reconstructed_mask.csv");
        result.final_valid_analysis_mask_csv_path = geometryArtifactPath(config, "_final_valid_analysis_mask.csv");
        result.non_crossing_adjustment_mask_csv_path = geometryArtifactPath(config, "_non_crossing_adjustment_mask.csv");
        result.summary_csv_path = geometryArtifactPath(config, "_stage6_summary.csv");

        if (!writeStage6FieldCsv(
                result, result.outer_reconstructed_csv_path, "z_outer_reconstructed", result.z_outer_reconstructed)) {
            throw std::runtime_error("Failed to write Stage 6 outer reconstructed CSV");
        }
        if (!writeStage6FieldCsv(
                result, result.inner_reconstructed_csv_path, "z_inner_reconstructed", result.z_inner_reconstructed)) {
            throw std::runtime_error("Failed to write Stage 6 inner reconstructed CSV");
        }
        if (!writeStage6MaskCsv(result, result.reconstructed_mask_csv_path, "reconstructed", result.reconstructed_mask)) {
            throw std::runtime_error("Failed to write Stage 6 reconstructed mask CSV");
        }
        if (!writeStage6MaskCsv(result,
                                result.final_valid_analysis_mask_csv_path,
                                "final_valid_analysis",
                                result.final_valid_analysis_mask)) {
            throw std::runtime_error("Failed to write Stage 6 final-valid-analysis mask CSV");
        }
        if (!writeStage6MaskCsv(result,
                                result.non_crossing_adjustment_mask_csv_path,
                                "non_crossing_adjusted",
                                result.non_crossing_adjustment_mask)) {
            throw std::runtime_error("Failed to write Stage 6 non-crossing-adjustment mask CSV");
        }
        if (!writeStage6SummaryCsv(result)) {
            throw std::runtime_error("Failed to write Stage 6 summary CSV");
        }
    }

    if (config.stage6_export_obj_meshes) {
        const bool use_stl = config.stage6_mesh_export_format == FoldPatchAnalysisConfig::MeshExportFormat::stl;
        const std::string extension = use_stl ? ".stl" : ".obj";
        if (config.stage6_split_in_out_meshes) {
            result.outer_obj_path = geometryArtifactPath(config, "_outer_surface") + extension;
            result.inner_obj_path = geometryArtifactPath(config, "_inner_surface") + extension;
            if (use_stl) {
                const Stage6StlExportResult outer_mesh =
                    writeStage6StlMesh(result.grid, result.obj_vertex_mask, result.z_outer_reconstructed, result.outer_obj_path);
                const Stage6StlExportResult inner_mesh =
                    writeStage6StlMesh(result.grid, result.obj_vertex_mask, result.z_inner_reconstructed, result.inner_obj_path);
                result.outer_obj_vertex_count = outer_mesh.vertex_count;
                result.outer_obj_face_count = outer_mesh.face_count;
                result.inner_obj_vertex_count = inner_mesh.vertex_count;
                result.inner_obj_face_count = inner_mesh.face_count;
            } else {
                const Stage6ObjExportResult outer_obj =
                    writeStage6ObjMesh(result.grid, result.obj_vertex_mask, result.z_outer_reconstructed, result.outer_obj_path);
                const Stage6ObjExportResult inner_obj =
                    writeStage6ObjMesh(result.grid, result.obj_vertex_mask, result.z_inner_reconstructed, result.inner_obj_path);
                result.outer_obj_vertex_count = outer_obj.vertex_count;
                result.outer_obj_face_count = outer_obj.face_count;
                result.inner_obj_vertex_count = inner_obj.vertex_count;
                result.inner_obj_face_count = inner_obj.face_count;
            }
        } else {
            result.outer_obj_path = geometryArtifactPath(config, "_surface") + extension;
            result.inner_obj_path = result.outer_obj_path;
            if (use_stl) {
                const Stage6DualMeshExportResult combined =
                    writeStage6StlMeshesCombined(result.grid,
                                                 result.obj_vertex_mask,
                                                 result.z_outer_reconstructed,
                                                 result.z_inner_reconstructed,
                                                 result.outer_obj_path);
                result.outer_obj_vertex_count = combined.outer_vertex_count;
                result.outer_obj_face_count = combined.outer_face_count;
                result.inner_obj_vertex_count = combined.inner_vertex_count;
                result.inner_obj_face_count = combined.inner_face_count;
            } else {
                const Stage6DualMeshExportResult combined =
                    writeStage6ObjMeshesCombined(result.grid,
                                                 result.obj_vertex_mask,
                                                 result.z_outer_reconstructed,
                                                 result.z_inner_reconstructed,
                                                 result.outer_obj_path);
                result.outer_obj_vertex_count = combined.outer_vertex_count;
                result.outer_obj_face_count = combined.outer_face_count;
                result.inner_obj_vertex_count = combined.inner_vertex_count;
                result.inner_obj_face_count = combined.inner_face_count;
            }
        }

        result.outer_reconstructed_scalar_node_count =
            countStage6ReconstructedScalarNodes(result.grid, result.obj_vertex_mask, result.z_outer_reconstructed);
        result.inner_reconstructed_scalar_node_count =
            countStage6ReconstructedScalarNodes(result.grid, result.obj_vertex_mask, result.z_inner_reconstructed);

        const std::size_t outer_used_by_faces =
            countStage6ReconstructedNodesUsedByFaces(result.grid, result.obj_vertex_mask, result.z_outer_reconstructed);
        const std::size_t inner_used_by_faces =
            countStage6ReconstructedNodesUsedByFaces(result.grid, result.obj_vertex_mask, result.z_inner_reconstructed);
        result.outer_reconstructed_nodes_not_used_by_any_face =
            result.outer_reconstructed_scalar_node_count > outer_used_by_faces
                ? result.outer_reconstructed_scalar_node_count - outer_used_by_faces
                : 0;
        result.inner_reconstructed_nodes_not_used_by_any_face =
            result.inner_reconstructed_scalar_node_count > inner_used_by_faces
                ? result.inner_reconstructed_scalar_node_count - inner_used_by_faces
                : 0;
    }

    result.messages.push_back("Geometry Stage 6 enforce non-crossing: " +
                              std::to_string(config.stage6_enforce_non_crossing ? 1 : 0));
    result.messages.push_back("Geometry Stage 6 minimum separation: " + std::to_string(config.stage6_min_separation));
    result.messages.push_back("Geometry Stage 6 seed nodes: " + std::to_string(result.seed_node_count));
    result.messages.push_back("Geometry Stage 6 interpolation nodes: " + std::to_string(result.interp_node_count));
    result.messages.push_back("Geometry Stage 6 reconstructed nodes: " + std::to_string(result.reconstructed_node_count));
    result.messages.push_back("Geometry Stage 6 unresolved nodes: " + std::to_string(result.unresolved_node_count));
    result.messages.push_back("Geometry Stage 6 outer iterations used: " + std::to_string(result.outer_iterations_used));
    result.messages.push_back("Geometry Stage 6 inner iterations used: " + std::to_string(result.inner_iterations_used));
    result.messages.push_back("Geometry Stage 6 outer final max update: " + std::to_string(result.outer_final_max_update));
    result.messages.push_back("Geometry Stage 6 inner final max update: " + std::to_string(result.inner_final_max_update));
    result.messages.push_back("Geometry Stage 6 non-crossing adjusted nodes: " +
                              std::to_string(result.non_crossing_adjusted_node_count));
    result.messages.push_back("Geometry Stage 6 min reconstructed separation: " +
                              std::to_string(result.min_reconstructed_separation));
    result.messages.push_back("Geometry Stage 6 max reconstructed separation: " +
                              std::to_string(result.max_reconstructed_separation));
    result.messages.push_back("Geometry Stage 6 mean reconstructed separation: " +
                              std::to_string(result.mean_reconstructed_separation));
    if (config.debug) {
        result.messages.push_back("Geometry Stage 6 outer reconstructed CSV: " + result.outer_reconstructed_csv_path);
        result.messages.push_back("Geometry Stage 6 inner reconstructed CSV: " + result.inner_reconstructed_csv_path);
        result.messages.push_back("Geometry Stage 6 reconstructed mask CSV: " + result.reconstructed_mask_csv_path);
        result.messages.push_back("Geometry Stage 6 final-valid-analysis mask CSV: " +
                                  result.final_valid_analysis_mask_csv_path);
        result.messages.push_back("Geometry Stage 6 non-crossing-adjustment mask CSV: " +
                                  result.non_crossing_adjustment_mask_csv_path);
        result.messages.push_back("Geometry Stage 6 summary CSV: " + result.summary_csv_path);
    }
    if (config.stage6_export_obj_meshes) {
        const std::string mesh_label =
            config.stage6_mesh_export_format == FoldPatchAnalysisConfig::MeshExportFormat::stl ? "STL" : "OBJ";
        if (config.stage6_split_in_out_meshes) {
            result.messages.push_back("Geometry Stage 6 outer " + mesh_label + ": " + result.outer_obj_path +
                                      " (reconstructed scalar nodes=" +
                                      std::to_string(result.outer_reconstructed_scalar_node_count) +
                                      ", mesh vertices emitted=" + std::to_string(result.outer_obj_vertex_count) +
                                      ", mesh faces emitted=" + std::to_string(result.outer_obj_face_count) +
                                      ", reconstructed nodes not used by any face=" +
                                      std::to_string(result.outer_reconstructed_nodes_not_used_by_any_face) + ")");
            result.messages.push_back("Geometry Stage 6 inner " + mesh_label + ": " + result.inner_obj_path +
                                      " (reconstructed scalar nodes=" +
                                      std::to_string(result.inner_reconstructed_scalar_node_count) +
                                      ", mesh vertices emitted=" + std::to_string(result.inner_obj_vertex_count) +
                                      ", mesh faces emitted=" + std::to_string(result.inner_obj_face_count) +
                                      ", reconstructed nodes not used by any face=" +
                                      std::to_string(result.inner_reconstructed_nodes_not_used_by_any_face) + ")");
        } else {
            result.messages.push_back("Geometry Stage 6 " + mesh_label + ": " + result.outer_obj_path +
                                      " (outer reconstructed scalar nodes=" +
                                      std::to_string(result.outer_reconstructed_scalar_node_count) +
                                      ", outer mesh vertices emitted=" + std::to_string(result.outer_obj_vertex_count) +
                                      ", outer mesh faces emitted=" + std::to_string(result.outer_obj_face_count) +
                                      ", outer reconstructed nodes not used by any face=" +
                                      std::to_string(result.outer_reconstructed_nodes_not_used_by_any_face) +
                                      ", inner reconstructed scalar nodes=" +
                                      std::to_string(result.inner_reconstructed_scalar_node_count) +
                                      ", inner mesh vertices emitted=" + std::to_string(result.inner_obj_vertex_count) +
                                      ", inner mesh faces emitted=" + std::to_string(result.inner_obj_face_count) +
                                      ", inner reconstructed nodes not used by any face=" +
                                      std::to_string(result.inner_reconstructed_nodes_not_used_by_any_face) + ")");
        }
    }
    result.messages.push_back("Geometry analysis: completed Stage 6 smooth surface reconstruction.");

    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryStage7SmoothedSurfaceResult runGeometryAnalysisStage7SurfaceSmoothing(
    const GeometryStage6SurfaceReconstructionResult& stage6_result,
    const GeometryStage5SurfacePrepResult& stage5_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance) {
    (void)tolerance;
    GeometryStage7SmoothedSurfaceResult result;
    if (!stage6_result.success) {
        throw std::runtime_error("Stage 7 cannot run before successful Stage 6 surface reconstruction");
    }
    if (stage6_result.node_count == 0 || stage6_result.grid.nx == 0 || stage6_result.grid.ny == 0) {
        throw std::runtime_error("Stage 7 requires a non-empty Stage 6 grid");
    }
    if (stage6_result.reconstructed_node_count == 0) {
        throw std::runtime_error("Stage 7 requires Stage 6 reconstructed nodes");
    }
    if (stage5_result.reliable_core_node_count == 0) {
        throw std::runtime_error("Stage 7 requires Stage 5 reliable-core nodes");
    }
    if (config.stage7_max_iterations == 0) {
        throw std::runtime_error("Stage 7 requires stage7_max_iterations > 0");
    }
    if (config.stage7_convergence_tolerance <= 0.0) {
        throw std::runtime_error("Stage 7 requires stage7_convergence_tolerance > 0");
    }
    if (config.stage7_smoothing_weight <= 0.0) {
        throw std::runtime_error("Stage 7 requires stage7_smoothing_weight > 0");
    }
    if (config.stage7_lambda <= 0.0) {
        throw std::runtime_error("Stage 7 requires stage7_lambda > 0");
    }
    if (config.stage7_data_weight_seed < 0.0 || config.stage7_data_weight_interp < 0.0) {
        throw std::runtime_error("Stage 7 requires non-negative fidelity weights");
    }
    if (config.stage7_data_weight_seed <= 0.0 && config.stage7_data_weight_interp <= 0.0) {
        throw std::runtime_error("Stage 7 requires at least one positive fidelity weight");
    }
    if (config.stage7_solver_max_iterations == 0) {
        throw std::runtime_error("Stage 7 requires stage7_solver_max_iterations > 0");
    }
    if (config.stage7_solver_tolerance <= 0.0) {
        throw std::runtime_error("Stage 7 requires stage7_solver_tolerance > 0");
    }

    result.messages.push_back("Geometry Stage 7");
    result.messages.push_back("Geometry analysis: starting Stage 7 surface smoothing / regularization.");
    result.grid = stage6_result.grid;
    result.node_count = stage6_result.node_count;
    result.reconstructed_mask = stage6_result.reconstructed_mask;
    result.reliable_core_mask = stage5_result.reliable_core_mask;
    result.smooth_valid_mask.assign(result.node_count, 0);
    result.metric_domain_mask.assign(result.node_count, 0);
    result.smooth_non_crossing_adjustment_mask.assign(result.node_count, 0);
    result.stage7_method_label = stage7MethodLabel(config.stage7_method);
    result.stage7_lambda = config.stage7_lambda;
    result.stage7_data_weight_seed = config.stage7_data_weight_seed;
    result.stage7_data_weight_interp = config.stage7_data_weight_interp;
    result.stage7_data_weight_policy =
        "paired_seed -> stage7_data_weight_seed ; reconstructed_nonseed -> stage7_data_weight_interp";
    result.stage7_use_reliable_core_only_for_fit = config.stage7_use_reliable_core_only_for_fit;
    result.stage7_boundary_condition_mode_label = stage7BoundaryConditionModeLabel(config.stage7_boundary_condition_mode);

    if (config.stage7_method == FoldPatchAnalysisConfig::Stage7Method::smooth) {
        const Stage7FieldSmoothResult outer = runStage7FieldSmoothingLegacy(result.grid,
                                                                             stage6_result.z_outer_reconstructed,
                                                                             stage6_result.reconstructed_mask,
                                                                             stage5_result.paired_seed_mask,
                                                                             config);
        const Stage7FieldSmoothResult inner = runStage7FieldSmoothingLegacy(result.grid,
                                                                             stage6_result.z_inner_reconstructed,
                                                                             stage6_result.reconstructed_mask,
                                                                             stage5_result.paired_seed_mask,
                                                                             config);
        result.z_outer_smooth = outer.values;
        result.z_inner_smooth = inner.values;
        result.outer_iterations_used = outer.iterations_used;
        result.inner_iterations_used = inner.iterations_used;
        result.outer_final_max_update = outer.final_max_update;
        result.inner_final_max_update = inner.final_max_update;
    } else {
        const std::vector<uint8_t> fit_domain_mask = buildStage7FitDomainMask(stage6_result, stage5_result, config);
        result.stage7_fit_active_node_count = 0;
        for (const uint8_t v : fit_domain_mask) {
            if (v != 0) {
                ++result.stage7_fit_active_node_count;
            }
        }
        if (result.stage7_fit_active_node_count == 0) {
            throw std::runtime_error("Stage 7 thin-plate fit has zero active fit-domain nodes");
        }
        const std::vector<uint8_t> boundary_mask = buildStage7BoundaryMask(result.grid, fit_domain_mask);
        const std::vector<double> fidelity_weights = buildStage7FidelityWeights(
            fit_domain_mask, stage5_result.paired_seed_mask, config, result.stage7_fit_seed_like_node_count,
            result.stage7_fit_interp_like_node_count);

        const Stage7FieldFitResult outer = runStage7FieldThinPlateGridFit(result.grid,
                                                                           stage6_result.z_outer_reconstructed,
                                                                           fit_domain_mask,
                                                                           boundary_mask,
                                                                           fidelity_weights,
                                                                           config);
        const Stage7FieldFitResult inner = runStage7FieldThinPlateGridFit(result.grid,
                                                                           stage6_result.z_inner_reconstructed,
                                                                           fit_domain_mask,
                                                                           boundary_mask,
                                                                           fidelity_weights,
                                                                           config);
        result.z_outer_smooth = outer.values;
        result.z_inner_smooth = inner.values;
        result.outer_iterations_used = outer.iterations_used;
        result.inner_iterations_used = inner.iterations_used;
        result.outer_final_max_update = outer.final_update;
        result.inner_final_max_update = inner.final_update;
        result.outer_fit_final_residual = outer.final_residual;
        result.inner_fit_final_residual = inner.final_residual;
        result.outer_fit_max_abs_residual = outer.max_abs_residual;
        result.inner_fit_max_abs_residual = inner.max_abs_residual;
        result.outer_fit_bending_energy = outer.bending_energy;
        result.inner_fit_bending_energy = inner.bending_energy;
        result.outer_solver_iterations_used = outer.iterations_used;
        result.inner_solver_iterations_used = inner.iterations_used;
        result.outer_solver_final_update = outer.final_update;
        result.inner_solver_final_update = inner.final_update;
    }

    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (stage6_result.reconstructed_mask[idx] == 0 || !std::isfinite(result.z_outer_smooth[idx]) ||
            !std::isfinite(result.z_inner_smooth[idx])) {
            continue;
        }
        if (config.stage7_enforce_non_crossing &&
            result.z_outer_smooth[idx] < (result.z_inner_smooth[idx] + config.stage7_min_separation)) {
            const double mid = 0.5 * (result.z_outer_smooth[idx] + result.z_inner_smooth[idx]);
            result.z_inner_smooth[idx] = mid - (0.5 * config.stage7_min_separation);
            result.z_outer_smooth[idx] = mid + (0.5 * config.stage7_min_separation);
            result.smooth_non_crossing_adjustment_mask[idx] = 1;
            ++result.smooth_non_crossing_adjusted_node_count;
        }
        result.smooth_valid_mask[idx] = 1;
        ++result.smooth_valid_node_count;
        if (result.reliable_core_mask[idx] != 0) {
            result.metric_domain_mask[idx] = 1;
            ++result.metric_domain_node_count;
        }
    }
    if (result.metric_domain_node_count == 0) {
        throw std::runtime_error("Stage 7 produced zero metric-domain nodes after reliable-core restriction");
    }

    bool first_sep = true;
    double sum_sep = 0.0;
    std::size_t sep_count = 0;
    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.metric_domain_mask[idx] == 0) {
            continue;
        }
        const double separation = result.z_outer_smooth[idx] - result.z_inner_smooth[idx];
        if (!std::isfinite(separation)) {
            continue;
        }
        if (first_sep) {
            result.min_smooth_separation = separation;
            result.max_smooth_separation = separation;
            first_sep = false;
        } else {
            result.min_smooth_separation = std::min(result.min_smooth_separation, separation);
            result.max_smooth_separation = std::max(result.max_smooth_separation, separation);
        }
        sum_sep += separation;
        ++sep_count;
    }
    if (sep_count == 0) {
        throw std::runtime_error("Stage 7 could not compute finite smooth separation on metric domain");
    }
    result.mean_smooth_separation = sum_sep / static_cast<double>(sep_count);

    if (config.stage7_export_s7s6_deltas) {
        double sum_abs_outer_delta = 0.0;
        double sum_abs_inner_delta = 0.0;
        double sum_abs_thickness_delta = 0.0;
        std::size_t finite_outer_delta_count = 0;
        std::size_t finite_inner_delta_count = 0;
        std::size_t finite_thickness_delta_count = 0;
        for (std::size_t idx = 0; idx < result.node_count; ++idx) {
            const double z_outer_stage6 = stage6_result.z_outer_reconstructed[idx];
            const double z_outer_stage7 = result.z_outer_smooth[idx];
            if (std::isfinite(z_outer_stage6) && std::isfinite(z_outer_stage7)) {
                const double abs_delta = std::fabs(z_outer_stage7 - z_outer_stage6);
                result.max_abs_outer_delta = std::max(result.max_abs_outer_delta, abs_delta);
                sum_abs_outer_delta += abs_delta;
                ++finite_outer_delta_count;
            }

            const double z_inner_stage6 = stage6_result.z_inner_reconstructed[idx];
            const double z_inner_stage7 = result.z_inner_smooth[idx];
            if (std::isfinite(z_inner_stage6) && std::isfinite(z_inner_stage7)) {
                const double abs_delta = std::fabs(z_inner_stage7 - z_inner_stage6);
                result.max_abs_inner_delta = std::max(result.max_abs_inner_delta, abs_delta);
                sum_abs_inner_delta += abs_delta;
                ++finite_inner_delta_count;
            }

            if (std::isfinite(z_outer_stage6) && std::isfinite(z_inner_stage6) && std::isfinite(z_outer_stage7) &&
                std::isfinite(z_inner_stage7)) {
                const double t_stage6 = z_outer_stage6 - z_inner_stage6;
                const double t_stage7 = z_outer_stage7 - z_inner_stage7;
                const double abs_delta = std::fabs(t_stage7 - t_stage6);
                result.max_abs_thickness_delta = std::max(result.max_abs_thickness_delta, abs_delta);
                sum_abs_thickness_delta += abs_delta;
                ++finite_thickness_delta_count;
            }
        }
        result.mean_abs_outer_delta =
            finite_outer_delta_count > 0 ? (sum_abs_outer_delta / static_cast<double>(finite_outer_delta_count)) : 0.0;
        result.mean_abs_inner_delta =
            finite_inner_delta_count > 0 ? (sum_abs_inner_delta / static_cast<double>(finite_inner_delta_count)) : 0.0;
        result.mean_abs_thickness_delta =
            finite_thickness_delta_count > 0 ? (sum_abs_thickness_delta / static_cast<double>(finite_thickness_delta_count))
                                             : 0.0;
    }

    if (config.stage7_export_s7s6_deltas) {
        result.outer_delta_csv_path = geometryArtifactPath(config, "_s7s6_outer_delta.csv");
        result.inner_delta_csv_path = geometryArtifactPath(config, "_s7s6_inner_delta.csv");
        result.thickness_delta_csv_path = geometryArtifactPath(config, "_s7s6_thickness_delta.csv");
        if (!writeStage7DeltaCsv(result.grid,
                                 stage5_result.inside_disk_mask,
                                 stage5_result,
                                 stage6_result,
                                 result,
                                 Stage7DeltaCsvKind::outer,
                                 result.outer_delta_csv_path,
                                 result.max_abs_outer_delta,
                                 result.mean_abs_outer_delta)) {
            throw std::runtime_error("Failed to write Stage 7-vs-Stage 6 outer delta CSV");
        }
        if (!writeStage7DeltaCsv(result.grid,
                                 stage5_result.inside_disk_mask,
                                 stage5_result,
                                 stage6_result,
                                 result,
                                 Stage7DeltaCsvKind::inner,
                                 result.inner_delta_csv_path,
                                 result.max_abs_inner_delta,
                                 result.mean_abs_inner_delta)) {
            throw std::runtime_error("Failed to write Stage 7-vs-Stage 6 inner delta CSV");
        }
        if (!writeStage7DeltaCsv(result.grid,
                                 stage5_result.inside_disk_mask,
                                 stage5_result,
                                 stage6_result,
                                 result,
                                 Stage7DeltaCsvKind::thickness,
                                 result.thickness_delta_csv_path,
                                 result.max_abs_thickness_delta,
                                 result.mean_abs_thickness_delta)) {
            throw std::runtime_error("Failed to write Stage 7-vs-Stage 6 thickness delta CSV");
        }
    }

    if (config.debug) {
        result.outer_smooth_csv_path = geometryArtifactPath(config, "_outer_smooth.csv");
        result.inner_smooth_csv_path = geometryArtifactPath(config, "_inner_smooth.csv");
        result.smooth_valid_mask_csv_path = geometryArtifactPath(config, "_smooth_valid_mask.csv");
        result.metric_domain_mask_csv_path = geometryArtifactPath(config, "_smooth_metric_domain_mask.csv");
        result.smooth_non_crossing_adjustment_mask_csv_path =
            geometryArtifactPath(config, "_smooth_non_crossing_adjustment_mask.csv");
        result.summary_csv_path = geometryArtifactPath(config, "_stage7_summary.csv");
        if (!writeStage7FieldCsv(result, result.outer_smooth_csv_path, "z_outer_smooth", result.z_outer_smooth)) {
            throw std::runtime_error("Failed to write Stage 7 outer smooth CSV");
        }
        if (!writeStage7FieldCsv(result, result.inner_smooth_csv_path, "z_inner_smooth", result.z_inner_smooth)) {
            throw std::runtime_error("Failed to write Stage 7 inner smooth CSV");
        }
        if (!writeStage7MaskCsv(result, result.smooth_valid_mask_csv_path, "smooth_valid", result.smooth_valid_mask)) {
            throw std::runtime_error("Failed to write Stage 7 smooth-valid mask CSV");
        }
        if (!writeStage7MaskCsv(
                result, result.metric_domain_mask_csv_path, "metric_domain", result.metric_domain_mask)) {
            throw std::runtime_error("Failed to write Stage 7 metric-domain mask CSV");
        }
        if (!writeStage7MaskCsv(result,
                                result.smooth_non_crossing_adjustment_mask_csv_path,
                                "smooth_non_crossing_adjusted",
                                result.smooth_non_crossing_adjustment_mask)) {
            throw std::runtime_error("Failed to write Stage 7 non-crossing-adjustment mask CSV");
        }
        if (!writeStage7SummaryCsv(result, config)) {
            throw std::runtime_error("Failed to write Stage 7 summary CSV");
        }
    }

    if (config.stage7_export_meshes) {
        const bool use_stl = config.stage6_mesh_export_format == FoldPatchAnalysisConfig::MeshExportFormat::stl;
        const std::string extension = use_stl ? ".stl" : ".obj";
        const std::string mesh_prefix = config.stage7_method == FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit
                                            ? "_thin_plate"
                                            : "_smooth";
        if (config.stage6_split_in_out_meshes) {
            result.outer_mesh_path = geometryArtifactPath(config, "stage7" + mesh_prefix + "_outer_surface" + extension);
            result.inner_mesh_path = geometryArtifactPath(config, "stage7" + mesh_prefix + "_inner_surface" + extension);
            if (use_stl) {
                const Stage6StlExportResult outer_mesh =
                    writeStage6StlMesh(result.grid, result.metric_domain_mask, result.z_outer_smooth, result.outer_mesh_path);
                const Stage6StlExportResult inner_mesh =
                    writeStage6StlMesh(result.grid, result.metric_domain_mask, result.z_inner_smooth, result.inner_mesh_path);
                result.outer_mesh_vertex_count = outer_mesh.vertex_count;
                result.outer_mesh_face_count = outer_mesh.face_count;
                result.inner_mesh_vertex_count = inner_mesh.vertex_count;
                result.inner_mesh_face_count = inner_mesh.face_count;
            } else {
                const Stage6ObjExportResult outer_mesh =
                    writeStage6ObjMesh(result.grid, result.metric_domain_mask, result.z_outer_smooth, result.outer_mesh_path);
                const Stage6ObjExportResult inner_mesh =
                    writeStage6ObjMesh(result.grid, result.metric_domain_mask, result.z_inner_smooth, result.inner_mesh_path);
                result.outer_mesh_vertex_count = outer_mesh.vertex_count;
                result.outer_mesh_face_count = outer_mesh.face_count;
                result.inner_mesh_vertex_count = inner_mesh.vertex_count;
                result.inner_mesh_face_count = inner_mesh.face_count;
            }
        } else {
            result.outer_mesh_path = geometryArtifactPath(config, "stage7" + mesh_prefix + "_surface" + extension);
            result.inner_mesh_path = result.outer_mesh_path;
            if (use_stl) {
                const Stage6DualMeshExportResult combined = writeStage6StlMeshesCombined(
                    result.grid, result.metric_domain_mask, result.z_outer_smooth, result.z_inner_smooth, result.outer_mesh_path);
                result.outer_mesh_vertex_count = combined.outer_vertex_count;
                result.outer_mesh_face_count = combined.outer_face_count;
                result.inner_mesh_vertex_count = combined.inner_vertex_count;
                result.inner_mesh_face_count = combined.inner_face_count;
            } else {
                const Stage6DualMeshExportResult combined = writeStage6ObjMeshesCombined(
                    result.grid, result.metric_domain_mask, result.z_outer_smooth, result.z_inner_smooth, result.outer_mesh_path);
                result.outer_mesh_vertex_count = combined.outer_vertex_count;
                result.outer_mesh_face_count = combined.outer_face_count;
                result.inner_mesh_vertex_count = combined.inner_vertex_count;
                result.inner_mesh_face_count = combined.inner_face_count;
            }
        }
        const std::size_t outer_scalar_nodes =
            countStage6ReconstructedScalarNodes(result.grid, result.metric_domain_mask, result.z_outer_smooth);
        const std::size_t inner_scalar_nodes =
            countStage6ReconstructedScalarNodes(result.grid, result.metric_domain_mask, result.z_inner_smooth);
        const std::size_t outer_used_by_faces =
            countStage6ReconstructedNodesUsedByFaces(result.grid, result.metric_domain_mask, result.z_outer_smooth);
        const std::size_t inner_used_by_faces =
            countStage6ReconstructedNodesUsedByFaces(result.grid, result.metric_domain_mask, result.z_inner_smooth);
        result.outer_mesh_unused_scalar_nodes =
            outer_scalar_nodes > outer_used_by_faces ? outer_scalar_nodes - outer_used_by_faces : 0;
        result.inner_mesh_unused_scalar_nodes =
            inner_scalar_nodes > inner_used_by_faces ? inner_scalar_nodes - inner_used_by_faces : 0;
    }

    result.messages.push_back("Geometry Stage 7 method: " + result.stage7_method_label);
    result.messages.push_back("Geometry Stage 7 smoothing weight: " + std::to_string(config.stage7_smoothing_weight));
    result.messages.push_back("Geometry Stage 7 smoothing max iterations: " +
                              std::to_string(config.stage7_max_iterations));
    result.messages.push_back("Geometry Stage 7 thin-plate solver max iterations: " +
                              std::to_string(config.stage7_solver_max_iterations));
    result.messages.push_back("Geometry Stage 7 convergence tolerance: " +
                              std::to_string(config.stage7_convergence_tolerance));
    result.messages.push_back("Geometry Stage 7 preserve seed values: " +
                              std::to_string(config.stage7_preserve_seed_values ? 1 : 0));
    result.messages.push_back("Geometry Stage 7 lambda: " + std::to_string(config.stage7_lambda));
    result.messages.push_back("Geometry Stage 7 data weight seed: " + std::to_string(config.stage7_data_weight_seed));
    result.messages.push_back("Geometry Stage 7 data weight interp: " + std::to_string(config.stage7_data_weight_interp));
    result.messages.push_back("Geometry Stage 7 boundary mode: " + result.stage7_boundary_condition_mode_label);
    result.messages.push_back("Geometry Stage 7 fit active node count: " + std::to_string(result.stage7_fit_active_node_count));
    result.messages.push_back("Geometry Stage 7 fit seed-like node count: " +
                              std::to_string(result.stage7_fit_seed_like_node_count));
    result.messages.push_back("Geometry Stage 7 fit interp-like node count: " +
                              std::to_string(result.stage7_fit_interp_like_node_count));
    result.messages.push_back("Geometry Stage 7 enforce non-crossing: " +
                              std::to_string(config.stage7_enforce_non_crossing ? 1 : 0));
    result.messages.push_back("Geometry Stage 7 minimum separation: " + std::to_string(config.stage7_min_separation));
    const std::string stage7_iteration_counter_label =
        config.stage7_method == FoldPatchAnalysisConfig::Stage7Method::smooth ? "smoothing iterations used"
                                                                               : "solver iterations used";
    result.messages.push_back("Geometry Stage 7 outer " + stage7_iteration_counter_label + ": " +
                              std::to_string(result.outer_iterations_used));
    result.messages.push_back("Geometry Stage 7 inner " + stage7_iteration_counter_label + ": " +
                              std::to_string(result.inner_iterations_used));
    result.messages.push_back("Geometry Stage 7 outer final update: " + formatScientific(result.outer_final_max_update));
    result.messages.push_back("Geometry Stage 7 inner final update: " + formatScientific(result.inner_final_max_update));
    result.messages.push_back("Geometry Stage 7 outer fit residual: " + formatScientific(result.outer_fit_final_residual));
    result.messages.push_back("Geometry Stage 7 inner fit residual: " + formatScientific(result.inner_fit_final_residual));
    result.messages.push_back("Geometry Stage 7 smooth valid node count: " + std::to_string(result.smooth_valid_node_count));
    result.messages.push_back("Geometry Stage 7 metric domain node count: " + std::to_string(result.metric_domain_node_count));
    result.messages.push_back("Geometry Stage 7 smooth non-crossing adjusted node count: " +
                              std::to_string(result.smooth_non_crossing_adjusted_node_count));
    result.messages.push_back("Geometry Stage 7 min smooth separation: " + std::to_string(result.min_smooth_separation));
    result.messages.push_back("Geometry Stage 7 max smooth separation: " + std::to_string(result.max_smooth_separation));
    result.messages.push_back("Geometry Stage 7 mean smooth separation: " + std::to_string(result.mean_smooth_separation));
    if (config.stage7_export_s7s6_deltas) {
        result.messages.push_back("Geometry Stage 7 max abs outer delta: " + formatScientific(result.max_abs_outer_delta));
        result.messages.push_back("Geometry Stage 7 max abs inner delta: " + formatScientific(result.max_abs_inner_delta));
        result.messages.push_back("Geometry Stage 7 max abs thickness delta: " +
                                  formatScientific(result.max_abs_thickness_delta));
        result.messages.push_back("Geometry Stage 7-vs-Stage 6 outer delta CSV: " + result.outer_delta_csv_path);
        result.messages.push_back("Geometry Stage 7-vs-Stage 6 inner delta CSV: " + result.inner_delta_csv_path);
        result.messages.push_back("Geometry Stage 7-vs-Stage 6 thickness delta CSV: " + result.thickness_delta_csv_path);
    }
    if (config.debug) {
        result.messages.push_back("Geometry Stage 7 outer smooth CSV: " + result.outer_smooth_csv_path);
        result.messages.push_back("Geometry Stage 7 inner smooth CSV: " + result.inner_smooth_csv_path);
        result.messages.push_back("Geometry Stage 7 smooth-valid mask CSV: " + result.smooth_valid_mask_csv_path);
        result.messages.push_back("Geometry Stage 7 metric-domain mask CSV: " + result.metric_domain_mask_csv_path);
        result.messages.push_back("Geometry Stage 7 non-crossing-adjustment mask CSV: " +
                                  result.smooth_non_crossing_adjustment_mask_csv_path);
        result.messages.push_back("Geometry Stage 7 summary CSV: " + result.summary_csv_path);
    }
    if (config.stage7_export_meshes) {
        const std::string mesh_label =
            config.stage6_mesh_export_format == FoldPatchAnalysisConfig::MeshExportFormat::stl ? "STL" : "OBJ";
        if (config.stage6_split_in_out_meshes) {
            result.messages.push_back("Geometry Stage 7 outer " + mesh_label + ": " + result.outer_mesh_path +
                                      " (mesh vertices emitted=" + std::to_string(result.outer_mesh_vertex_count) +
                                      ", mesh faces emitted=" + std::to_string(result.outer_mesh_face_count) +
                                      ", scalar nodes not used by any face=" +
                                      std::to_string(result.outer_mesh_unused_scalar_nodes) + ")");
            result.messages.push_back("Geometry Stage 7 inner " + mesh_label + ": " + result.inner_mesh_path +
                                      " (mesh vertices emitted=" + std::to_string(result.inner_mesh_vertex_count) +
                                      ", mesh faces emitted=" + std::to_string(result.inner_mesh_face_count) +
                                      ", scalar nodes not used by any face=" +
                                      std::to_string(result.inner_mesh_unused_scalar_nodes) + ")");
        } else {
            result.messages.push_back("Geometry Stage 7 " + mesh_label + ": " + result.outer_mesh_path +
                                      " (outer mesh vertices emitted=" + std::to_string(result.outer_mesh_vertex_count) +
                                      ", outer mesh faces emitted=" + std::to_string(result.outer_mesh_face_count) +
                                      ", outer scalar nodes not used by any face=" +
                                      std::to_string(result.outer_mesh_unused_scalar_nodes) +
                                      ", inner mesh vertices emitted=" + std::to_string(result.inner_mesh_vertex_count) +
                                      ", inner mesh faces emitted=" + std::to_string(result.inner_mesh_face_count) +
                                      ", inner scalar nodes not used by any face=" +
                                      std::to_string(result.inner_mesh_unused_scalar_nodes) + ")");
        }
    }
    result.messages.push_back("Geometry analysis: completed Stage 7 surface smoothing / regularization.");

    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

namespace {

enum class Stage8FailureReason : uint8_t {
    none = 0,
    insufficient_points = 1,
    boundary_neighbor_geometry = 2,
    rank_deficient = 3,
    poor_conditioning = 4,
    high_residual = 5
};

struct Stage8SupportNode {
    std::size_t idx = 0;
    double u = 0.0;
    double v = 0.0;
    double dist2 = 0.0;
};

struct Stage8QuadraticFitResult {
    bool solved = false;
    bool rank_deficient = false;
    double condition_indicator = std::numeric_limits<double>::infinity();
    double rms_residual = std::numeric_limits<double>::quiet_NaN();
    double max_abs_residual = std::numeric_limits<double>::quiet_NaN();
    // z(u,v) = a + b*u + c*v + d*u^2 + e*u*v + f*v^2
    // Derivatives at center (u=v=0):
    // dz/dx = b, dz/dy = c, d2z/dx2 = 2d, d2z/dy2 = 2f, d2z/dxdy = e.
    double a = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    double c = std::numeric_limits<double>::quiet_NaN();
    double d = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double f = std::numeric_limits<double>::quiet_NaN();
};

double stage8DefaultFitRadius(const Stage4GridDescriptor& grid) { return 2.5 * grid.spacing; }
double stage8DefaultMinDirectionalSpan(const Stage4GridDescriptor& grid) { return 1.5 * grid.spacing; }

Stage8QuadraticFitResult fitLocalQuadraticLeastSquares(const std::vector<Stage8SupportNode>& support,
                                                       const std::vector<double>& z_values,
                                                       double tolerance) {
    Stage8QuadraticFitResult result;
    if (support.size() < 6) {
        result.rank_deficient = true;
        return result;
    }

    double ata[6][6] = {};
    double atz[6] = {};
    for (const Stage8SupportNode& node : support) {
        const double row[6] = {1.0, node.u, node.v, node.u * node.u, node.u * node.v, node.v * node.v};
        const double z = z_values[node.idx];
        for (int r = 0; r < 6; ++r) {
            atz[r] += row[r] * z;
            for (int c = 0; c < 6; ++c) {
                ata[r][c] += row[r] * row[c];
            }
        }
    }

    double scale = 0.0;
    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 6; ++c) {
            scale = std::max(scale, std::fabs(ata[r][c]));
        }
    }
    if (scale <= tolerance) {
        result.rank_deficient = true;
        return result;
    }

    double aug[6][7] = {};
    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 6; ++c) {
            aug[r][c] = ata[r][c];
        }
        aug[r][6] = atz[r];
    }

    double min_pivot = std::numeric_limits<double>::infinity();
    double max_pivot = 0.0;
    for (int col = 0; col < 6; ++col) {
        int pivot_row = col;
        double pivot_abs = std::fabs(aug[col][col]);
        for (int r = col + 1; r < 6; ++r) {
            const double candidate = std::fabs(aug[r][col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = r;
            }
        }
        const double pivot_tol = tolerance * std::max(1.0, scale);
        if (pivot_abs <= pivot_tol) {
            result.rank_deficient = true;
            return result;
        }
        if (pivot_row != col) {
            for (int c = col; c < 7; ++c) {
                std::swap(aug[col][c], aug[pivot_row][c]);
            }
        }
        min_pivot = std::min(min_pivot, pivot_abs);
        max_pivot = std::max(max_pivot, pivot_abs);

        for (int r = col + 1; r < 6; ++r) {
            const double factor = aug[r][col] / aug[col][col];
            if (factor == 0.0) {
                continue;
            }
            for (int c = col; c < 7; ++c) {
                aug[r][c] -= factor * aug[col][c];
            }
        }
    }

    if (!(min_pivot > 0.0) || !std::isfinite(min_pivot) || !std::isfinite(max_pivot)) {
        result.rank_deficient = true;
        return result;
    }
    result.condition_indicator = max_pivot / min_pivot;

    double x[6] = {};
    for (int r = 5; r >= 0; --r) {
        double rhs = aug[r][6];
        for (int c = r + 1; c < 6; ++c) {
            rhs -= aug[r][c] * x[c];
        }
        x[r] = rhs / aug[r][r];
    }

    result.a = x[0];
    result.b = x[1];
    result.c = x[2];
    result.d = x[3];
    result.e = x[4];
    result.f = x[5];

    double sum_sq_res = 0.0;
    double max_abs_res = 0.0;
    for (const Stage8SupportNode& node : support) {
        const double pred = result.a + (result.b * node.u) + (result.c * node.v) + (result.d * node.u * node.u) +
                            (result.e * node.u * node.v) + (result.f * node.v * node.v);
        const double res = z_values[node.idx] - pred;
        sum_sq_res += res * res;
        max_abs_res = std::max(max_abs_res, std::fabs(res));
    }
    result.rms_residual = std::sqrt(sum_sq_res / static_cast<double>(support.size()));
    result.max_abs_residual = max_abs_res;
    result.solved = true;
    return result;
}

bool writeStage8DerivativeCsv(const GeometryStage8DerivativeEstimationResult& result,
                              const std::string& path,
                              const std::vector<double>& dz_dx,
                              const std::vector<double>& dz_dy,
                              const std::vector<double>& d2z_dx2,
                              const std::vector<double>& d2z_dy2,
                              const std::vector<double>& d2z_dxdy,
                              const std::vector<double>& fit_rms,
                              const std::vector<double>& fit_max_abs,
                              const std::vector<double>& fit_condition) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,reconstructed,reliable_core,metric_domain,derivative_fit_attempted,derivative_valid,"
           "neighbor_count,neighbor_max_radius,fit_rms_residual,fit_max_abs_residual,fit_condition_indicator,"
           "dz_dx,dz_dy,d2z_dx2,d2z_dy2,d2z_dxdy\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.reconstructed_mask[idx]) << ',' << static_cast<int>(result.reliable_core_mask[idx])
                << ',' << static_cast<int>(result.metric_domain_mask[idx]) << ','
                << static_cast<int>(result.derivative_fit_attempted_mask[idx]) << ','
                << static_cast<int>(result.derivative_valid_mask[idx]) << ',' << result.derivative_neighbor_count[idx] << ','
                << result.derivative_neighbor_max_radius[idx] << ',' << fit_rms[idx] << ',' << fit_max_abs[idx] << ','
                << fit_condition[idx] << ',' << dz_dx[idx] << ',' << dz_dy[idx] << ',' << d2z_dx2[idx] << ','
                << d2z_dy2[idx] << ',' << d2z_dxdy[idx] << '\n';
        }
    }
    return out.good();
}

bool writeStage8DerivativeValidMaskCsv(const GeometryStage8DerivativeEstimationResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,metric_domain,derivative_fit_attempted,derivative_valid\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.metric_domain_mask[idx]) << ','
                << static_cast<int>(result.derivative_fit_attempted_mask[idx]) << ','
                << static_cast<int>(result.derivative_valid_mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage8DerivativeFailureReasonCsv(const GeometryStage8DerivativeEstimationResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,metric_domain,fit_attempted,invalid_insufficient_points,invalid_rank_deficient,"
           "invalid_poor_conditioning,invalid_high_residual,invalid_boundary_neighbor_geometry\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.metric_domain_mask[idx]) << ','
                << static_cast<int>(result.derivative_fit_attempted_mask[idx]) << ','
                << static_cast<int>(result.derivative_invalid_insufficient_points_mask[idx]) << ','
                << static_cast<int>(result.derivative_invalid_rank_deficient_mask[idx]) << ','
                << static_cast<int>(result.derivative_invalid_poor_conditioning_mask[idx]) << ','
                << static_cast<int>(result.derivative_invalid_high_residual_mask[idx]) << ','
                << static_cast<int>(result.derivative_invalid_boundary_neighbor_geometry_mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage8SummaryCsv(const GeometryStage8DerivativeEstimationResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "node_count,metric_domain_node_count,derivative_fit_attempted_node_count,derivative_valid_node_count,"
           "derivative_invalid_node_count,derivative_invalid_insufficient_points_count,"
           "derivative_invalid_rank_deficient_count,derivative_invalid_poor_conditioning_count,"
           "derivative_invalid_high_residual_count,derivative_invalid_boundary_neighbor_geometry_count,"
           "mean_outer_rms_residual_valid,mean_inner_rms_residual_valid,mean_outer_condition_indicator_valid,"
           "mean_inner_condition_indicator_valid,mean_neighbor_count_valid\n";

    double sum_outer_rms = 0.0;
    double sum_inner_rms = 0.0;
    double sum_outer_cond = 0.0;
    double sum_inner_cond = 0.0;
    double sum_neighbor_count = 0.0;
    std::size_t count = 0;
    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.derivative_valid_mask[idx] == 0) {
            continue;
        }
        sum_outer_rms += result.outer_fit_rms_residual[idx];
        sum_inner_rms += result.inner_fit_rms_residual[idx];
        sum_outer_cond += result.outer_fit_condition_indicator[idx];
        sum_inner_cond += result.inner_fit_condition_indicator[idx];
        sum_neighbor_count += static_cast<double>(result.derivative_neighbor_count[idx]);
        ++count;
    }
    const double denom = count > 0 ? static_cast<double>(count) : 1.0;
    out << result.node_count << ',' << result.metric_domain_node_count << ',' << result.derivative_fit_attempted_node_count
        << ',' << result.derivative_valid_node_count << ',' << result.derivative_invalid_node_count << ','
        << result.derivative_invalid_insufficient_points_count << ',' << result.derivative_invalid_rank_deficient_count
        << ',' << result.derivative_invalid_poor_conditioning_count << ',' << result.derivative_invalid_high_residual_count
        << ',' << result.derivative_invalid_boundary_neighbor_geometry_count << ',' << (sum_outer_rms / denom) << ','
        << (sum_inner_rms / denom) << ',' << (sum_outer_cond / denom) << ',' << (sum_inner_cond / denom) << ','
        << (sum_neighbor_count / denom) << '\n';
    return out.good();
}

struct Stage9ScalarStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double median = std::numeric_limits<double>::quiet_NaN();
    double stddev = std::numeric_limits<double>::quiet_NaN();
    double min = std::numeric_limits<double>::quiet_NaN();
    double max = std::numeric_limits<double>::quiet_NaN();
};

double medianOfSortedValues(const std::vector<double>& sorted_values) {
    if (sorted_values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const std::size_t n = sorted_values.size();
    if ((n % 2U) == 1U) {
        return sorted_values[n / 2U];
    }
    return 0.5 * (sorted_values[(n / 2U) - 1U] + sorted_values[n / 2U]);
}

double medianOfValues(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(values.begin(), values.end());
    return medianOfSortedValues(values);
}

double robustSigmaFromValues(const std::vector<double>& values, double median) {
    if (values.empty() || !std::isfinite(median)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<double> abs_deviations;
    abs_deviations.reserve(values.size());
    for (const double value : values) {
        if (!std::isfinite(value)) {
            continue;
        }
        abs_deviations.push_back(std::abs(value - median));
    }
    if (abs_deviations.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double mad = medianOfValues(abs_deviations);
    return 1.4826 * mad;
}

Stage9ScalarStats computeStage9ScalarStats(const std::vector<double>& values, const std::vector<uint8_t>& valid_mask) {
    Stage9ScalarStats stats;
    std::vector<double> selected;
    selected.reserve(values.size());
    for (std::size_t idx = 0; idx < values.size(); ++idx) {
        if (valid_mask[idx] == 0 || !std::isfinite(values[idx])) {
            continue;
        }
        selected.push_back(values[idx]);
    }
    if (selected.empty()) {
        return stats;
    }

    double sum = 0.0;
    stats.min = selected[0];
    stats.max = selected[0];
    for (const double value : selected) {
        sum += value;
        stats.min = std::min(stats.min, value);
        stats.max = std::max(stats.max, value);
    }
    stats.mean = sum / static_cast<double>(selected.size());

    std::sort(selected.begin(), selected.end());
    const std::size_t n = selected.size();
    if ((n % 2U) == 1U) {
        stats.median = selected[n / 2U];
    } else {
        stats.median = 0.5 * (selected[(n / 2U) - 1U] + selected[n / 2U]);
    }

    double variance_sum = 0.0;
    for (const double value : selected) {
        const double delta = value - stats.mean;
        variance_sum += delta * delta;
    }
    stats.stddev = std::sqrt(variance_sum / static_cast<double>(n));
    return stats;
}

bool writeStage9CurvatureCsv(const GeometryStage9CurvatureComputationResult& result,
                             const std::string& path,
                             const std::vector<double>& raw_H,
                             const std::vector<double>& oriented_H,
                             const std::vector<double>& raw_K,
                             const std::vector<double>& graph_normal_dot_radial,
                             const std::vector<uint8_t>& orientation_flip_applied,
                             const std::vector<uint8_t>& global_tail_flag,
                             const std::vector<uint8_t>& local_spike_flag,
                             const std::vector<uint8_t>& qc_warn_flag,
                             const std::vector<uint8_t>& confidence_class) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,reconstructed,reliable_core,metric_domain,derivative_valid,curvature_valid,H_raw,H_oriented,K_raw,"
           "graph_normal_dot_radial,orientation_flip_applied,K_global_tail_flag,K_local_spike_flag,K_qc_warn_flag,"
           "K_confidence_class\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.reconstructed_mask[idx]) << ',' << static_cast<int>(result.reliable_core_mask[idx])
                << ',' << static_cast<int>(result.metric_domain_mask[idx]) << ','
                << static_cast<int>(result.derivative_valid_mask[idx]) << ','
                << static_cast<int>(result.curvature_valid_mask[idx]) << ',' << raw_H[idx] << ',' << oriented_H[idx] << ','
                << raw_K[idx] << ',' << graph_normal_dot_radial[idx] << ','
                << static_cast<int>(orientation_flip_applied[idx]) << ','
                << static_cast<int>(global_tail_flag[idx]) << ',' << static_cast<int>(local_spike_flag[idx]) << ','
                << static_cast<int>(qc_warn_flag[idx]) << ',' << static_cast<int>(confidence_class[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage9CurvatureMaskCsv(const GeometryStage9CurvatureComputationResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,metric_domain,derivative_valid,curvature_valid,invalid_nonfinite_input,invalid_nonfinite_output\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.metric_domain_mask[idx]) << ','
                << static_cast<int>(result.derivative_valid_mask[idx]) << ','
                << static_cast<int>(result.curvature_valid_mask[idx]) << ','
                << static_cast<int>(result.curvature_invalid_nonfinite_input_mask[idx]) << ','
                << static_cast<int>(result.curvature_invalid_nonfinite_output_mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage9SummaryCsv(const GeometryStage9CurvatureComputationResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "node_count,metric_domain_node_count,derivative_valid_node_count,curvature_valid_node_count,"
           "curvature_valid_fraction_of_metric_domain,curvature_valid_fraction_of_derivative_valid,"
           "curvature_invalid_nonfinite_input_count,curvature_invalid_nonfinite_output_count,"
           "outer_K_global_tail_count,outer_K_local_spike_count,outer_K_qc_warn_count,"
           "inner_K_global_tail_count,inner_K_local_spike_count,inner_K_qc_warn_count,"
           "outer_K_global_tail_fraction_of_curvature_valid,outer_K_local_spike_fraction_of_curvature_valid,"
           "outer_K_qc_warn_fraction_of_curvature_valid,inner_K_global_tail_fraction_of_curvature_valid,"
           "inner_K_local_spike_fraction_of_curvature_valid,inner_K_qc_warn_fraction_of_curvature_valid,"
           "outer_mean_H,outer_median_H,outer_stddev_H,outer_min_H,outer_max_H,"
           "outer_mean_oriented_H,outer_median_oriented_H,outer_stddev_oriented_H,outer_min_oriented_H,"
           "outer_max_oriented_H,"
           "outer_mean_K,outer_median_K,outer_stddev_K,outer_min_K,outer_max_K,"
           "inner_mean_H,inner_median_H,inner_stddev_H,inner_min_H,inner_max_H,"
           "inner_mean_oriented_H,inner_median_oriented_H,inner_stddev_oriented_H,inner_min_oriented_H,"
           "inner_max_oriented_H,"
           "inner_mean_K,inner_median_K,inner_stddev_K,inner_min_K,inner_max_K\n";
    out << result.node_count << ',' << result.metric_domain_node_count << ',' << result.derivative_valid_node_count << ','
        << result.curvature_valid_node_count << ',' << result.curvature_valid_fraction_of_metric_domain << ','
        << result.curvature_valid_fraction_of_derivative_valid << ',' << result.curvature_invalid_nonfinite_input_count
        << ',' << result.curvature_invalid_nonfinite_output_count << ',' << result.outer_K_global_tail_count << ','
        << result.outer_K_local_spike_count << ',' << result.outer_K_qc_warn_count << ','
        << result.inner_K_global_tail_count << ',' << result.inner_K_local_spike_count << ','
        << result.inner_K_qc_warn_count << ',' << result.outer_K_global_tail_fraction_of_curvature_valid << ','
        << result.outer_K_local_spike_fraction_of_curvature_valid << ','
        << result.outer_K_qc_warn_fraction_of_curvature_valid << ','
        << result.inner_K_global_tail_fraction_of_curvature_valid << ','
        << result.inner_K_local_spike_fraction_of_curvature_valid << ','
        << result.inner_K_qc_warn_fraction_of_curvature_valid << ',' << result.outer_mean_H << ','
        << result.outer_median_H << ',' << result.outer_stddev_H << ',' << result.outer_min_H << ',' << result.outer_max_H
        << ',' << result.outer_mean_oriented_H << ',' << result.outer_median_oriented_H << ','
        << result.outer_stddev_oriented_H << ',' << result.outer_min_oriented_H << ',' << result.outer_max_oriented_H
        << ',' << result.outer_mean_K << ',' << result.outer_median_K << ',' << result.outer_stddev_K << ','
        << result.outer_min_K << ',' << result.outer_max_K << ',' << result.inner_mean_H << ',' << result.inner_median_H
        << ',' << result.inner_stddev_H << ',' << result.inner_min_H << ',' << result.inner_max_H << ','
        << result.inner_mean_oriented_H << ',' << result.inner_median_oriented_H << ','
        << result.inner_stddev_oriented_H << ',' << result.inner_min_oriented_H << ',' << result.inner_max_oriented_H << ','
        << result.inner_mean_K << ',' << result.inner_median_K << ',' << result.inner_stddev_K << ','
        << result.inner_min_K << ',' << result.inner_max_K << '\n';
    return out.good();
}



struct Stage10ScalarStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double median = std::numeric_limits<double>::quiet_NaN();
    double stddev = std::numeric_limits<double>::quiet_NaN();
    double min = std::numeric_limits<double>::quiet_NaN();
    double max = std::numeric_limits<double>::quiet_NaN();
};

bool sampleStage7ScalarFieldBilinear(const Stage4GridDescriptor& grid,
                                     const std::vector<double>& field,
                                     const std::vector<uint8_t>& valid_mask,
                                     double x,
                                     double y,
                                     double* out_value) {
    if (out_value == nullptr || grid.nx < 2 || grid.ny < 2 || field.size() != grid.nx * grid.ny ||
        valid_mask.size() != field.size()) {
        return false;
    }
    if (!std::isfinite(x) || !std::isfinite(y) || x < grid.x_min || x > grid.x_max || y < grid.y_min || y > grid.y_max) {
        return false;
    }

    const auto x_it_hi = std::upper_bound(grid.x_values.begin(), grid.x_values.end(), x);
    const auto y_it_hi = std::upper_bound(grid.y_values.begin(), grid.y_values.end(), y);
    if (x_it_hi == grid.x_values.begin() || y_it_hi == grid.y_values.begin() || x_it_hi == grid.x_values.end() ||
        y_it_hi == grid.y_values.end()) {
        return false;
    }
    const std::size_t i1 = static_cast<std::size_t>(std::distance(grid.x_values.begin(), x_it_hi));
    const std::size_t j1 = static_cast<std::size_t>(std::distance(grid.y_values.begin(), y_it_hi));
    const std::size_t i0 = i1 - 1;
    const std::size_t j0 = j1 - 1;
    const double x0 = grid.x_values[i0];
    const double x1 = grid.x_values[i1];
    const double y0 = grid.y_values[j0];
    const double y1 = grid.y_values[j1];
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    if (dx <= 0.0 || dy <= 0.0) {
        return false;
    }

    const std::size_t idx00 = stage4NodeIndex(i0, j0, grid.nx);
    const std::size_t idx10 = stage4NodeIndex(i1, j0, grid.nx);
    const std::size_t idx01 = stage4NodeIndex(i0, j1, grid.nx);
    const std::size_t idx11 = stage4NodeIndex(i1, j1, grid.nx);
    const double z00 = field[idx00];
    const double z10 = field[idx10];
    const double z01 = field[idx01];
    const double z11 = field[idx11];
    if (valid_mask[idx00] == 0 || valid_mask[idx10] == 0 || valid_mask[idx01] == 0 || valid_mask[idx11] == 0 ||
        !std::isfinite(z00) || !std::isfinite(z10) || !std::isfinite(z01) || !std::isfinite(z11)) {
        return false;
    }

    const double tx = (x - x0) / dx;
    const double ty = (y - y0) / dy;
    *out_value = (1.0 - tx) * (1.0 - ty) * z00 + tx * (1.0 - ty) * z10 + (1.0 - tx) * ty * z01 + tx * ty * z11;
    return true;
}

enum class Stage10RadialIntersectionStatus : uint8_t {
    success = 0,
    outside_inner_sampling_domain = 1,
    no_bracket = 2,
    root_failure = 3
};

struct Stage10RadialIntersectionResult {
    Stage10RadialIntersectionStatus status = Stage10RadialIntersectionStatus::root_failure;
    double s_inner = std::numeric_limits<double>::quiet_NaN();
};

Stage10RadialIntersectionResult findRadialInnerIntersectionOnStage7Surface(const Stage4GridDescriptor& grid,
                                                                            const std::vector<double>& z_inner,
                                                                            const std::vector<uint8_t>& metric_mask,
                                                                            double ux,
                                                                            double uy,
                                                                            double uz,
                                                                            double s_out,
                                                                            double tolerance) {
    Stage10RadialIntersectionResult result;
    if (!std::isfinite(ux) || !std::isfinite(uy) || !std::isfinite(uz) || !std::isfinite(s_out) || s_out <= tolerance) {
        return result;
    }
    auto eval_f = [&](double s, double* out_f) -> bool {
        const double x = s * ux;
        const double y = s * uy;
        double z_inner_interp = std::numeric_limits<double>::quiet_NaN();
        if (!sampleStage7ScalarFieldBilinear(grid, z_inner, metric_mask, x, y, &z_inner_interp)) {
            return false;
        }
        *out_f = s * uz - z_inner_interp;
        return std::isfinite(*out_f);
    };

    double f_hi = 0.0;
    if (!eval_f(s_out, &f_hi)) {
        result.status = Stage10RadialIntersectionStatus::outside_inner_sampling_domain;
        return result;
    }
    if (std::fabs(f_hi) <= 1e-9) {
        result.status = Stage10RadialIntersectionStatus::success;
        result.s_inner = s_out;
        return result;
    }

    constexpr std::size_t kScanSteps = 64;
    const double s_min = std::max(0.0, 0.01 * s_out);
    bool saw_domain_failure = false;
    bool bracket_found = false;
    double s_a = s_out;
    double s_b = s_out;
    double f_a = f_hi;
    double f_b = f_hi;

    double prev_s = s_out;
    double prev_f = f_hi;
    for (std::size_t step = 1; step <= kScanSteps; ++step) {
        const double alpha = static_cast<double>(step) / static_cast<double>(kScanSteps);
        const double s = s_out - alpha * (s_out - s_min);
        double f_s = 0.0;
        if (!eval_f(s, &f_s)) {
            saw_domain_failure = true;
            continue;
        }
        if (std::fabs(f_s) <= 1e-9) {
            result.status = Stage10RadialIntersectionStatus::success;
            result.s_inner = s;
            return result;
        }
        if ((prev_f > 0.0 && f_s < 0.0) || (prev_f < 0.0 && f_s > 0.0)) {
            s_a = prev_s;
            s_b = s;
            f_a = prev_f;
            f_b = f_s;
            bracket_found = true;
            break;
        }
        prev_s = s;
        prev_f = f_s;
    }

    if (!bracket_found) {
        result.status = saw_domain_failure ? Stage10RadialIntersectionStatus::outside_inner_sampling_domain
                                           : Stage10RadialIntersectionStatus::no_bracket;
        return result;
    }

    constexpr std::size_t kMaxBisectionIterations = 64;
    for (std::size_t iter = 0; iter < kMaxBisectionIterations; ++iter) {
        const double s_mid = 0.5 * (s_a + s_b);
        double f_mid = 0.0;
        if (!eval_f(s_mid, &f_mid)) {
            result.status = Stage10RadialIntersectionStatus::outside_inner_sampling_domain;
            return result;
        }
        if (std::fabs(f_mid) <= 1e-9 || std::fabs(s_b - s_a) <= std::max(1e-8, tolerance)) {
            result.status = Stage10RadialIntersectionStatus::success;
            result.s_inner = s_mid;
            return result;
        }
        if ((f_a > 0.0 && f_mid < 0.0) || (f_a < 0.0 && f_mid > 0.0)) {
            s_b = s_mid;
            f_b = f_mid;
        } else {
            s_a = s_mid;
            f_a = f_mid;
        }
    }
    (void)f_b;
    result.status = Stage10RadialIntersectionStatus::root_failure;
    return result;
}

Stage10ScalarStats computeStage10ScalarStats(const std::vector<double>& values, const std::vector<uint8_t>& valid_mask) {
    Stage10ScalarStats stats;
    std::vector<double> selected;
    selected.reserve(values.size());
    for (std::size_t idx = 0; idx < values.size(); ++idx) {
        if (valid_mask[idx] == 0 || !std::isfinite(values[idx])) {
            continue;
        }
        selected.push_back(values[idx]);
    }
    if (selected.empty()) {
        return stats;
    }

    double sum = 0.0;
    stats.min = selected[0];
    stats.max = selected[0];
    for (double value : selected) {
        sum += value;
        stats.min = std::min(stats.min, value);
        stats.max = std::max(stats.max, value);
    }
    stats.mean = sum / static_cast<double>(selected.size());

    std::sort(selected.begin(), selected.end());
    stats.median = medianOfSortedValues(selected);

    if (selected.size() >= 2) {
        double variance_sum = 0.0;
        for (double value : selected) {
            const double delta = value - stats.mean;
            variance_sum += delta * delta;
        }
        stats.stddev = std::sqrt(variance_sum / static_cast<double>(selected.size() - 1));
    }
    return stats;
}

bool writeStage10ThicknessCsv(const GeometryStage10ThicknessResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,z_outer_smooth,in_metric_domain,curvature_valid,thickness_attempted,thickness_valid,thickness_vertical,thickness_radial,s_outer,s_inner\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << result.thickness_outer_z_stage7[idx] << ','
                << static_cast<int>(result.metric_domain_mask[idx]) << ','
                << static_cast<int>(result.curvature_valid_mask[idx]) << ','
                << static_cast<int>(result.thickness_attempted_mask[idx]) << ','
                << static_cast<int>(result.thickness_valid_mask[idx]) << ',' << result.thickness_vertical[idx] << ','
                << result.thickness_radial[idx] << ',' << result.thickness_radial_s_outer[idx] << ','
                << result.thickness_radial_s_inner[idx] << '\n';
        }
    }
    return out.good();
}

bool writeStage10ValidMaskCsv(const GeometryStage10ThicknessResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,thickness_valid\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.thickness_valid_mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage10InvalidReasonCsv(const GeometryStage10ThicknessResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,x,y,outside_domain,outside_inner_sampling_domain,no_bracket,root_failure,nonfinite_surface,negative_or_zero,below_min_threshold,above_max_threshold,qc_warn\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            out << i << ',' << j << ',' << result.grid.x_values[i] << ',' << result.grid.y_values[j] << ','
                << static_cast<int>(result.thickness_invalid_outside_domain_mask[idx]) << ','
                << static_cast<int>(result.thickness_radial_invalid_outside_inner_domain_mask[idx]) << ','
                << static_cast<int>(result.thickness_radial_invalid_no_bracket_mask[idx]) << ','
                << static_cast<int>(result.thickness_radial_invalid_root_failure_mask[idx]) << ','
                << static_cast<int>(result.thickness_invalid_nonfinite_surface_mask[idx]) << ','
                << static_cast<int>(result.thickness_invalid_negative_or_zero_mask[idx]) << ','
                << static_cast<int>(result.thickness_invalid_below_min_threshold_mask[idx]) << ','
                << static_cast<int>(result.thickness_invalid_above_max_threshold_mask[idx]) << ','
                << static_cast<int>(result.thickness_qc_warn_mask[idx]) << '\n';
        }
    }
    return out.good();
}

bool writeStage10SummaryCsv(const GeometryStage10ThicknessResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "stage10_method,thickness_method_label,local_thickness_definition,thickness_input_surface_definition,"
           "thickness_contract_note,thickness_domain_definition,thickness_domain_limitation_note,"
           "radial_invalid_outside_inner_domain_definition,node_count,thickness_excluded_outside_domain_count,"
           "thickness_attempted_node_count,thickness_valid_node_count,thickness_attempted_invalid_node_count,"
           "thickness_invalid_nonfinite_surface_count,thickness_invalid_negative_or_zero_count,"
           "thickness_invalid_below_min_threshold_count,thickness_invalid_above_max_threshold_count,"
           "thickness_qc_warn_count,thickness_valid_and_curvature_valid_node_count,"
           "thickness_valid_fraction_of_metric_domain,thickness_valid_intersection_fraction_of_metric_domain,"
           "thickness_radial_valid_fraction_of_metric_domain,"
           "thickness_radial_invalid_outside_inner_domain_fraction_of_metric_domain,"
           "thickness_radial_attempted_node_count,thickness_radial_valid_node_count,"
           "thickness_radial_invalid_no_bracket_count,thickness_radial_invalid_root_failure_count,"
           "thickness_radial_invalid_outside_inner_domain_count,mean_thickness_active,median_thickness_active,"
           "stddev_thickness_active,min_thickness_active,max_thickness_active\n";
    out << result.stage10_method << ',' << result.thickness_method_label << ',' << result.local_thickness_definition << ','
        << result.thickness_input_surface_definition << ',' << result.thickness_contract_note << ','
        << result.thickness_domain_definition << ',' << result.thickness_domain_limitation_note << ','
        << result.radial_invalid_outside_inner_domain_definition << ',' << result.node_count << ','
        << result.thickness_excluded_outside_domain_count << ',' << result.thickness_attempted_node_count << ','
        << result.thickness_valid_node_count << ',' << result.thickness_attempted_invalid_node_count << ','
        << result.thickness_invalid_nonfinite_surface_count << ','
        << result.thickness_invalid_negative_or_zero_count << ','
        << result.thickness_invalid_below_min_threshold_count << ','
        << result.thickness_invalid_above_max_threshold_count << ',' << result.thickness_qc_warn_count << ','
        << result.thickness_valid_and_curvature_valid_node_count << ','
        << result.thickness_valid_fraction_of_metric_domain << ','
        << result.thickness_valid_intersection_fraction_of_metric_domain << ','
        << result.thickness_radial_valid_fraction_of_metric_domain << ','
        << result.thickness_radial_invalid_outside_inner_domain_fraction_of_metric_domain << ','
        << result.thickness_radial_attempted_node_count << ',' << result.thickness_radial_valid_node_count << ','
        << result.thickness_radial_invalid_no_bracket_count << ','
        << result.thickness_radial_invalid_root_failure_count << ','
        << result.thickness_radial_invalid_outside_inner_domain_count << ',' << result.mean_thickness_active << ','
        << result.median_thickness_active << ',' << result.stddev_thickness_active << ','
        << result.min_thickness_active << ',' << result.max_thickness_active << '\n';
    return out.good();
}

bool writeStage10RadialPOutInCsv(const GeometryStage10ThicknessResult& result, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }
    out << "i,j,P_out x,P_out y,P_out z,P_in x,P_in y,P_in z,t_radial\n";
    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            const double x = result.grid.x_values[i];
            const double y = result.grid.y_values[j];
            const double z_out = result.thickness_outer_z_stage7[idx];
            const double s_out = result.thickness_radial_s_outer[idx];
            const double s_in = result.thickness_radial_s_inner[idx];
            double p_in_x = std::numeric_limits<double>::quiet_NaN();
            double p_in_y = std::numeric_limits<double>::quiet_NaN();
            double p_in_z = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(s_out) && s_out > 0.0 && std::isfinite(s_in)) {
                const double ux = x / s_out;
                const double uy = y / s_out;
                const double uz = z_out / s_out;
                p_in_x = s_in * ux;
                p_in_y = s_in * uy;
                p_in_z = s_in * uz;
            }
            out << i << ',' << j << ',' << x << ',' << y << ',' << z_out << ',' << p_in_x << ',' << p_in_y << ','
                << p_in_z << ',' << result.thickness_radial[idx] << '\n';
        }
    }
    return out.good();
}
} // namespace

GeometryStage8DerivativeEstimationResult runGeometryAnalysisStage8DerivativeEstimation(
    const GeometryStage7SmoothedSurfaceResult& stage7_result,
    const GeometryStage5SurfacePrepResult& stage5_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance) {
    GeometryStage8DerivativeEstimationResult result;
    if (!stage7_result.success) {
        throw std::runtime_error("Stage 8 cannot run before successful Stage 7 surface smoothing");
    }
    if (stage7_result.node_count == 0 || stage7_result.grid.nx == 0 || stage7_result.grid.ny == 0) {
        throw std::runtime_error("Stage 8 requires a non-empty Stage 7 grid");
    }
    if (config.stage8_min_points < 6) {
        throw std::runtime_error("Stage 8 requires stage8_min_points >= 6 for quadratic fitting");
    }
    if (config.stage8_max_rms_residual < 0.0 || config.stage8_max_abs_residual < 0.0 ||
        config.stage8_max_condition_indicator <= 0.0) {
        throw std::runtime_error("Stage 8 requires non-negative residual thresholds and positive condition threshold");
    }

    const std::size_t node_count = stage7_result.node_count;
    result.messages.push_back("Geometry Stage 8");
    result.messages.push_back("Geometry analysis: starting Stage 8 local quadratic derivative estimation.");
    result.grid = stage7_result.grid;
    result.node_count = node_count;
    result.reconstructed_mask = stage7_result.reconstructed_mask;
    result.reliable_core_mask = stage5_result.reliable_core_mask;
    result.metric_domain_mask = stage7_result.metric_domain_mask;

    const double nan = std::numeric_limits<double>::quiet_NaN();
    result.derivative_fit_attempted_mask.assign(node_count, 0);
    result.derivative_valid_mask.assign(node_count, 0);
    result.derivative_invalid_insufficient_points_mask.assign(node_count, 0);
    result.derivative_invalid_rank_deficient_mask.assign(node_count, 0);
    result.derivative_invalid_poor_conditioning_mask.assign(node_count, 0);
    result.derivative_invalid_high_residual_mask.assign(node_count, 0);
    result.derivative_invalid_boundary_neighbor_geometry_mask.assign(node_count, 0);

    result.outer_dz_dx.assign(node_count, nan);
    result.outer_dz_dy.assign(node_count, nan);
    result.outer_d2z_dx2.assign(node_count, nan);
    result.outer_d2z_dy2.assign(node_count, nan);
    result.outer_d2z_dxdy.assign(node_count, nan);
    result.inner_dz_dx.assign(node_count, nan);
    result.inner_dz_dy.assign(node_count, nan);
    result.inner_d2z_dx2.assign(node_count, nan);
    result.inner_d2z_dy2.assign(node_count, nan);
    result.inner_d2z_dxdy.assign(node_count, nan);

    result.outer_fit_rms_residual.assign(node_count, nan);
    result.inner_fit_rms_residual.assign(node_count, nan);
    result.outer_fit_max_abs_residual.assign(node_count, nan);
    result.inner_fit_max_abs_residual.assign(node_count, nan);
    result.outer_fit_condition_indicator.assign(node_count, nan);
    result.inner_fit_condition_indicator.assign(node_count, nan);
    result.derivative_neighbor_count.assign(node_count, 0);
    result.derivative_neighbor_max_radius.assign(node_count, nan);

    const double fit_radius = config.stage8_fit_radius > 0.0 ? config.stage8_fit_radius : stage8DefaultFitRadius(result.grid);
    const double min_directional_span = config.stage8_min_directional_span > 0.0
                                            ? config.stage8_min_directional_span
                                            : stage8DefaultMinDirectionalSpan(result.grid);
    const double fit_radius2 = fit_radius * fit_radius;

    for (const uint8_t m : result.metric_domain_mask) {
        if (m != 0) {
            ++result.metric_domain_node_count;
        }
    }

    for (std::size_t j0 = 0; j0 < result.grid.ny; ++j0) {
        for (std::size_t i0 = 0; i0 < result.grid.nx; ++i0) {
            const std::size_t center_idx = stage4NodeIndex(i0, j0, result.grid.nx);
            const bool candidate = result.metric_domain_mask[center_idx] != 0 && stage7_result.smooth_valid_mask[center_idx] != 0 &&
                                   std::isfinite(stage7_result.z_outer_smooth[center_idx]) &&
                                   std::isfinite(stage7_result.z_inner_smooth[center_idx]);
            if (!candidate) {
                continue;
            }

            result.derivative_fit_attempted_mask[center_idx] = 1;
            ++result.derivative_fit_attempted_node_count;

            const double x0 = result.grid.x_values[i0];
            const double y0 = result.grid.y_values[j0];
            std::vector<Stage8SupportNode> support;
            support.reserve(64);
            for (std::size_t j = 0; j < result.grid.ny; ++j) {
                const double y = result.grid.y_values[j];
                const double v = y - y0;
                for (std::size_t i = 0; i < result.grid.nx; ++i) {
                    const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
                    if (result.metric_domain_mask[idx] == 0 || stage7_result.smooth_valid_mask[idx] == 0 ||
                        !std::isfinite(stage7_result.z_outer_smooth[idx]) || !std::isfinite(stage7_result.z_inner_smooth[idx])) {
                        continue;
                    }
                    const double x = result.grid.x_values[i];
                    const double u = x - x0;
                    const double dist2 = (u * u) + (v * v);
                    if (dist2 <= fit_radius2 + tolerance) {
                        support.push_back({idx, u, v, dist2});
                    }
                }
            }

            std::sort(support.begin(), support.end(), [](const Stage8SupportNode& a, const Stage8SupportNode& b) {
                if (a.dist2 != b.dist2) {
                    return a.dist2 < b.dist2;
                }
                return a.idx < b.idx;
            });
            if (config.stage8_max_points > 0 && support.size() > config.stage8_max_points) {
                support.resize(config.stage8_max_points);
            }

            result.derivative_neighbor_count[center_idx] = static_cast<int>(support.size());
            result.derivative_neighbor_max_radius[center_idx] =
                support.empty() ? nan : std::sqrt(support.back().dist2);

            Stage8FailureReason reason = Stage8FailureReason::none;
            if (support.size() < config.stage8_min_points) {
                reason = Stage8FailureReason::insufficient_points;
            }

            bool has_pos_u = false;
            bool has_neg_u = false;
            bool has_pos_v = false;
            bool has_neg_v = false;
            double min_u = std::numeric_limits<double>::infinity();
            double max_u = -std::numeric_limits<double>::infinity();
            double min_v = std::numeric_limits<double>::infinity();
            double max_v = -std::numeric_limits<double>::infinity();
            for (const Stage8SupportNode& node : support) {
                min_u = std::min(min_u, node.u);
                max_u = std::max(max_u, node.u);
                min_v = std::min(min_v, node.v);
                max_v = std::max(max_v, node.v);
                has_pos_u = has_pos_u || (node.u > tolerance);
                has_neg_u = has_neg_u || (node.u < -tolerance);
                has_pos_v = has_pos_v || (node.v > tolerance);
                has_neg_v = has_neg_v || (node.v < -tolerance);
            }
            const bool centered_ok = (!config.stage8_require_centered_support) || (has_pos_u && has_neg_u && has_pos_v && has_neg_v);
            const bool span_ok = support.empty() ? false : ((max_u - min_u) >= min_directional_span &&
                                                            (max_v - min_v) >= min_directional_span);
            if (reason == Stage8FailureReason::none && (!centered_ok || !span_ok)) {
                reason = Stage8FailureReason::boundary_neighbor_geometry;
            }

            Stage8QuadraticFitResult outer_fit;
            Stage8QuadraticFitResult inner_fit;
            if (reason == Stage8FailureReason::none) {
                outer_fit = fitLocalQuadraticLeastSquares(support, stage7_result.z_outer_smooth, tolerance);
                inner_fit = fitLocalQuadraticLeastSquares(support, stage7_result.z_inner_smooth, tolerance);
                result.outer_fit_rms_residual[center_idx] = outer_fit.rms_residual;
                result.inner_fit_rms_residual[center_idx] = inner_fit.rms_residual;
                result.outer_fit_max_abs_residual[center_idx] = outer_fit.max_abs_residual;
                result.inner_fit_max_abs_residual[center_idx] = inner_fit.max_abs_residual;
                result.outer_fit_condition_indicator[center_idx] = outer_fit.condition_indicator;
                result.inner_fit_condition_indicator[center_idx] = inner_fit.condition_indicator;

                if (outer_fit.rank_deficient || inner_fit.rank_deficient || !outer_fit.solved || !inner_fit.solved) {
                    reason = Stage8FailureReason::rank_deficient;
                } else if (!std::isfinite(outer_fit.condition_indicator) || !std::isfinite(inner_fit.condition_indicator) ||
                           outer_fit.condition_indicator > config.stage8_max_condition_indicator ||
                           inner_fit.condition_indicator > config.stage8_max_condition_indicator) {
                    reason = Stage8FailureReason::poor_conditioning;
                } else if (!std::isfinite(outer_fit.rms_residual) || !std::isfinite(inner_fit.rms_residual) ||
                           !std::isfinite(outer_fit.max_abs_residual) || !std::isfinite(inner_fit.max_abs_residual) ||
                           outer_fit.rms_residual > config.stage8_max_rms_residual ||
                           inner_fit.rms_residual > config.stage8_max_rms_residual ||
                           outer_fit.max_abs_residual > config.stage8_max_abs_residual ||
                           inner_fit.max_abs_residual > config.stage8_max_abs_residual) {
                    reason = Stage8FailureReason::high_residual;
                }
            }

            if (reason == Stage8FailureReason::none) {
                result.derivative_valid_mask[center_idx] = 1;
                ++result.derivative_valid_node_count;
                result.outer_dz_dx[center_idx] = outer_fit.b;
                result.outer_dz_dy[center_idx] = outer_fit.c;
                result.outer_d2z_dx2[center_idx] = 2.0 * outer_fit.d;
                result.outer_d2z_dy2[center_idx] = 2.0 * outer_fit.f;
                result.outer_d2z_dxdy[center_idx] = outer_fit.e;
                result.inner_dz_dx[center_idx] = inner_fit.b;
                result.inner_dz_dy[center_idx] = inner_fit.c;
                result.inner_d2z_dx2[center_idx] = 2.0 * inner_fit.d;
                result.inner_d2z_dy2[center_idx] = 2.0 * inner_fit.f;
                result.inner_d2z_dxdy[center_idx] = inner_fit.e;
            } else {
                ++result.derivative_invalid_node_count;
                switch (reason) {
                case Stage8FailureReason::insufficient_points:
                    result.derivative_invalid_insufficient_points_mask[center_idx] = 1;
                    ++result.derivative_invalid_insufficient_points_count;
                    break;
                case Stage8FailureReason::boundary_neighbor_geometry:
                    result.derivative_invalid_boundary_neighbor_geometry_mask[center_idx] = 1;
                    ++result.derivative_invalid_boundary_neighbor_geometry_count;
                    break;
                case Stage8FailureReason::rank_deficient:
                    result.derivative_invalid_rank_deficient_mask[center_idx] = 1;
                    ++result.derivative_invalid_rank_deficient_count;
                    break;
                case Stage8FailureReason::poor_conditioning:
                    result.derivative_invalid_poor_conditioning_mask[center_idx] = 1;
                    ++result.derivative_invalid_poor_conditioning_count;
                    break;
                case Stage8FailureReason::high_residual:
                    result.derivative_invalid_high_residual_mask[center_idx] = 1;
                    ++result.derivative_invalid_high_residual_count;
                    break;
                case Stage8FailureReason::none:
                    break;
                }
            }
        }
    }

    if (config.stage8_export_csv) {
        result.outer_derivatives_csv_path = geometryArtifactPath(config, "_outer_derivatives.csv");
        result.inner_derivatives_csv_path = geometryArtifactPath(config, "_inner_derivatives.csv");
        result.derivative_valid_mask_csv_path = geometryArtifactPath(config, "_derivative_valid_mask.csv");
        result.derivative_failure_reason_csv_path = geometryArtifactPath(config, "_derivative_failure_reasons.csv");
        result.derivative_summary_csv_path = geometryArtifactPath(config, "_stage8_summary.csv");
        if (!writeStage8DerivativeCsv(result,
                                      result.outer_derivatives_csv_path,
                                      result.outer_dz_dx,
                                      result.outer_dz_dy,
                                      result.outer_d2z_dx2,
                                      result.outer_d2z_dy2,
                                      result.outer_d2z_dxdy,
                                      result.outer_fit_rms_residual,
                                      result.outer_fit_max_abs_residual,
                                      result.outer_fit_condition_indicator)) {
            throw std::runtime_error("Failed to write Stage 8 outer derivatives CSV");
        }
        if (!writeStage8DerivativeCsv(result,
                                      result.inner_derivatives_csv_path,
                                      result.inner_dz_dx,
                                      result.inner_dz_dy,
                                      result.inner_d2z_dx2,
                                      result.inner_d2z_dy2,
                                      result.inner_d2z_dxdy,
                                      result.inner_fit_rms_residual,
                                      result.inner_fit_max_abs_residual,
                                      result.inner_fit_condition_indicator)) {
            throw std::runtime_error("Failed to write Stage 8 inner derivatives CSV");
        }
        if (!writeStage8DerivativeValidMaskCsv(result, result.derivative_valid_mask_csv_path)) {
            throw std::runtime_error("Failed to write Stage 8 derivative-valid mask CSV");
        }
        if (!writeStage8DerivativeFailureReasonCsv(result, result.derivative_failure_reason_csv_path)) {
            throw std::runtime_error("Failed to write Stage 8 derivative failure-reason CSV");
        }
        if (!writeStage8SummaryCsv(result, result.derivative_summary_csv_path)) {
            throw std::runtime_error("Failed to write Stage 8 summary CSV");
        }
    }

    result.messages.push_back("Geometry Stage 8 fit radius: " + std::to_string(fit_radius));
    result.messages.push_back("Geometry Stage 8 minimum support points: " + std::to_string(config.stage8_min_points));
    result.messages.push_back("Geometry Stage 8 metric-domain node count: " + std::to_string(result.metric_domain_node_count));
    result.messages.push_back("Geometry Stage 8 fit-attempted node count: " +
                              std::to_string(result.derivative_fit_attempted_node_count));
    result.messages.push_back("Geometry Stage 8 derivative-valid node count: " +
                              std::to_string(result.derivative_valid_node_count));
    result.messages.push_back("Geometry Stage 8 invalid insufficient-points count: " +
                              std::to_string(result.derivative_invalid_insufficient_points_count));
    result.messages.push_back("Geometry Stage 8 invalid rank-deficient count: " +
                              std::to_string(result.derivative_invalid_rank_deficient_count));
    result.messages.push_back("Geometry Stage 8 invalid poor-conditioning count: " +
                              std::to_string(result.derivative_invalid_poor_conditioning_count));
    result.messages.push_back("Geometry Stage 8 invalid high-residual count: " +
                              std::to_string(result.derivative_invalid_high_residual_count));
    result.messages.push_back("Geometry Stage 8 invalid boundary-neighbor-geometry count: " +
                              std::to_string(result.derivative_invalid_boundary_neighbor_geometry_count));
    if (config.stage8_export_csv) {
        result.messages.push_back("Geometry Stage 8 outer derivatives CSV: " + result.outer_derivatives_csv_path);
        result.messages.push_back("Geometry Stage 8 inner derivatives CSV: " + result.inner_derivatives_csv_path);
        result.messages.push_back("Geometry Stage 8 derivative-valid mask CSV: " + result.derivative_valid_mask_csv_path);
        result.messages.push_back("Geometry Stage 8 derivative failure-reason CSV: " +
                                  result.derivative_failure_reason_csv_path);
        result.messages.push_back("Geometry Stage 8 summary CSV: " + result.derivative_summary_csv_path);
    }
    result.messages.push_back("Geometry analysis: completed Stage 8 local quadratic derivative estimation.");

    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryStage9CurvatureComputationResult runGeometryAnalysisStage9CurvatureComputation(
    const GeometryStage7SmoothedSurfaceResult& stage7_result,
    const GeometryStage8DerivativeEstimationResult& stage8_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance) {
    GeometryStage9CurvatureComputationResult result;
    if (!stage7_result.success) {
        throw std::runtime_error("Stage 9 cannot run before successful Stage 7 surface smoothing");
    }
    if (!stage8_result.success) {
        throw std::runtime_error("Stage 9 cannot run before successful Stage 8 derivative estimation");
    }
    if (stage8_result.node_count == 0 || stage8_result.grid.nx == 0 || stage8_result.grid.ny == 0) {
        throw std::runtime_error("Stage 9 requires a non-empty Stage 8 derivative grid");
    }
    if (config.stage9_qc_n_tail <= 0.0 || config.stage9_qc_n_spike <= 0.0) {
        throw std::runtime_error("Stage 9 QC thresholds must be positive");
    }
    if (config.stage9_qc_min_neighbors == 0) {
        throw std::runtime_error("Stage 9 QC minimum neighborhood size must be >= 1");
    }
    if (config.stage9_qc_abs_scale_floor <= 0.0) {
        throw std::runtime_error("Stage 9 QC absolute scale floor must be > 0");
    }

    result.messages.push_back("Geometry Stage 9");
    result.messages.push_back("Geometry analysis: starting Stage 9 curvature computation.");
    result.grid = stage8_result.grid;
    result.node_count = stage8_result.node_count;
    result.reconstructed_mask = stage8_result.reconstructed_mask;
    result.reliable_core_mask = stage8_result.reliable_core_mask;
    result.metric_domain_mask = stage8_result.metric_domain_mask;
    result.derivative_valid_mask = stage8_result.derivative_valid_mask;

    const double nan = std::numeric_limits<double>::quiet_NaN();
    result.curvature_valid_mask.assign(result.node_count, 0);
    result.curvature_invalid_nonfinite_input_mask.assign(result.node_count, 0);
    result.curvature_invalid_nonfinite_output_mask.assign(result.node_count, 0);
    result.outer_K_global_tail_flag.assign(result.node_count, 0);
    result.inner_K_global_tail_flag.assign(result.node_count, 0);
    result.outer_K_local_spike_flag.assign(result.node_count, 0);
    result.inner_K_local_spike_flag.assign(result.node_count, 0);
    result.outer_K_qc_warn_flag.assign(result.node_count, 0);
    result.inner_K_qc_warn_flag.assign(result.node_count, 0);
    result.outer_K_confidence_class.assign(result.node_count, 0);
    result.inner_K_confidence_class.assign(result.node_count, 0);
    result.outer_mean_curvature_H.assign(result.node_count, nan);
    result.outer_oriented_mean_curvature_H.assign(result.node_count, nan);
    result.outer_gaussian_curvature_K.assign(result.node_count, nan);
    result.outer_graph_normal_dot_radial.assign(result.node_count, nan);
    result.outer_outward_normal_alignment_flag.assign(result.node_count, 0);
    result.outer_orientation_flip_applied_flag.assign(result.node_count, 0);
    result.inner_mean_curvature_H.assign(result.node_count, nan);
    result.inner_oriented_mean_curvature_H.assign(result.node_count, nan);
    result.inner_gaussian_curvature_K.assign(result.node_count, nan);
    result.inner_graph_normal_dot_radial.assign(result.node_count, nan);
    result.inner_outward_normal_alignment_flag.assign(result.node_count, 0);
    result.inner_orientation_flip_applied_flag.assign(result.node_count, 0);

    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.metric_domain_mask[idx] != 0) {
            ++result.metric_domain_node_count;
        }
        if (result.derivative_valid_mask[idx] != 0) {
            ++result.derivative_valid_node_count;
        }
    }

    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.derivative_valid_mask[idx] == 0) {
            continue;
        }

        const bool finite_inputs = std::isfinite(stage8_result.outer_dz_dx[idx]) && std::isfinite(stage8_result.outer_dz_dy[idx]) &&
                                   std::isfinite(stage8_result.outer_d2z_dx2[idx]) &&
                                   std::isfinite(stage8_result.outer_d2z_dy2[idx]) &&
                                   std::isfinite(stage8_result.outer_d2z_dxdy[idx]) &&
                                   std::isfinite(stage8_result.inner_dz_dx[idx]) && std::isfinite(stage8_result.inner_dz_dy[idx]) &&
                                   std::isfinite(stage8_result.inner_d2z_dx2[idx]) &&
                                   std::isfinite(stage8_result.inner_d2z_dy2[idx]) &&
                                   std::isfinite(stage8_result.inner_d2z_dxdy[idx]);
        if (!finite_inputs) {
            result.curvature_invalid_nonfinite_input_mask[idx] = 1;
            ++result.curvature_invalid_nonfinite_input_count;
            continue;
        }

        const double outer_zx = stage8_result.outer_dz_dx[idx];
        const double outer_zy = stage8_result.outer_dz_dy[idx];
        const double outer_zxx = stage8_result.outer_d2z_dx2[idx];
        const double outer_zyy = stage8_result.outer_d2z_dy2[idx];
        const double outer_zxy = stage8_result.outer_d2z_dxdy[idx];

        const double inner_zx = stage8_result.inner_dz_dx[idx];
        const double inner_zy = stage8_result.inner_dz_dy[idx];
        const double inner_zxx = stage8_result.inner_d2z_dx2[idx];
        const double inner_zyy = stage8_result.inner_d2z_dy2[idx];
        const double inner_zxy = stage8_result.inner_d2z_dxdy[idx];

        // Monge-patch Gaussian curvature formula:
        // K = (z_xx z_yy - z_xy^2) / (1 + z_x^2 + z_y^2)^2
        const double outer_denom_base = 1.0 + (outer_zx * outer_zx) + (outer_zy * outer_zy);
        const double inner_denom_base = 1.0 + (inner_zx * inner_zx) + (inner_zy * inner_zy);
        const double outer_K =
            ((outer_zxx * outer_zyy) - (outer_zxy * outer_zxy)) / (outer_denom_base * outer_denom_base);
        const double inner_K =
            ((inner_zxx * inner_zyy) - (inner_zxy * inner_zxy)) / (inner_denom_base * inner_denom_base);

        // Monge-patch mean curvature formula:
        // H = ((1+z_x^2) z_yy - 2 z_x z_y z_xy + (1+z_y^2) z_xx) / (2(1+z_x^2+z_y^2)^(3/2))
        const double outer_H_num = ((1.0 + (outer_zx * outer_zx)) * outer_zyy) - (2.0 * outer_zx * outer_zy * outer_zxy) +
                                   ((1.0 + (outer_zy * outer_zy)) * outer_zxx);
        const double inner_H_num = ((1.0 + (inner_zx * inner_zx)) * inner_zyy) - (2.0 * inner_zx * inner_zy * inner_zxy) +
                                   ((1.0 + (inner_zy * inner_zy)) * inner_zxx);
        const double outer_H = outer_H_num / (2.0 * std::pow(outer_denom_base, 1.5));
        const double inner_H = inner_H_num / (2.0 * std::pow(inner_denom_base, 1.5));

        const std::size_t i = idx % result.grid.nx;
        const std::size_t j = idx / result.grid.nx;
        const double x = result.grid.x_values[i];
        const double y = result.grid.y_values[j];
        const double outer_z = stage7_result.z_outer_smooth[idx];
        const double inner_z = stage7_result.z_inner_smooth[idx];

        const double outer_graph_normal_x = -outer_zx;
        const double outer_graph_normal_y = -outer_zy;
        const double outer_graph_normal_z = 1.0;
        const double inner_graph_normal_x = -inner_zx;
        const double inner_graph_normal_y = -inner_zy;
        const double inner_graph_normal_z = 1.0;

        const double outer_graph_normal_norm = std::sqrt((outer_graph_normal_x * outer_graph_normal_x) +
                                                         (outer_graph_normal_y * outer_graph_normal_y) +
                                                         (outer_graph_normal_z * outer_graph_normal_z));
        const double inner_graph_normal_norm = std::sqrt((inner_graph_normal_x * inner_graph_normal_x) +
                                                         (inner_graph_normal_y * inner_graph_normal_y) +
                                                         (inner_graph_normal_z * inner_graph_normal_z));
        const double outer_radial_norm = std::sqrt((x * x) + (y * y) + (outer_z * outer_z));
        const double inner_radial_norm = std::sqrt((x * x) + (y * y) + (inner_z * inner_z));
        const double outer_graph_normal_dot_radial_raw =
            (outer_graph_normal_x * x) + (outer_graph_normal_y * y) + (outer_graph_normal_z * outer_z);
        const double inner_graph_normal_dot_radial_raw =
            (inner_graph_normal_x * x) + (inner_graph_normal_y * y) + (inner_graph_normal_z * inner_z);
        const double outer_graph_normal_dot_radial_denom = outer_graph_normal_norm * outer_radial_norm;
        const double inner_graph_normal_dot_radial_denom = inner_graph_normal_norm * inner_radial_norm;
        const double outer_graph_normal_dot_radial =
            (outer_graph_normal_dot_radial_denom > tolerance) ? (outer_graph_normal_dot_radial_raw / outer_graph_normal_dot_radial_denom) : 0.0;
        const double inner_graph_normal_dot_radial =
            (inner_graph_normal_dot_radial_denom > tolerance) ? (inner_graph_normal_dot_radial_raw / inner_graph_normal_dot_radial_denom) : 0.0;

        const bool outer_graph_normal_is_outward = outer_graph_normal_dot_radial > 0.0;
        const bool inner_graph_normal_is_outward = inner_graph_normal_dot_radial > 0.0;
        const double outer_oriented_H = outer_graph_normal_is_outward ? outer_H : -outer_H;
        const double inner_oriented_H = inner_graph_normal_is_outward ? inner_H : -inner_H;

        const bool finite_outputs =
            std::isfinite(outer_H) && std::isfinite(outer_K) && std::isfinite(inner_H) && std::isfinite(inner_K) &&
            std::isfinite(outer_graph_normal_dot_radial) && std::isfinite(inner_graph_normal_dot_radial) &&
            std::isfinite(outer_oriented_H) && std::isfinite(inner_oriented_H);
        if (!finite_outputs) {
            result.curvature_invalid_nonfinite_output_mask[idx] = 1;
            ++result.curvature_invalid_nonfinite_output_count;
            continue;
        }

        result.curvature_valid_mask[idx] = 1;
        ++result.curvature_valid_node_count;
        result.outer_mean_curvature_H[idx] = outer_H;
        result.outer_oriented_mean_curvature_H[idx] = outer_oriented_H;
        result.outer_gaussian_curvature_K[idx] = outer_K;
        result.outer_graph_normal_dot_radial[idx] = outer_graph_normal_dot_radial;
        result.outer_outward_normal_alignment_flag[idx] = static_cast<uint8_t>(outer_graph_normal_is_outward ? 1 : 0);
        result.outer_orientation_flip_applied_flag[idx] = static_cast<uint8_t>(outer_graph_normal_is_outward ? 0 : 1);
        result.inner_mean_curvature_H[idx] = inner_H;
        result.inner_oriented_mean_curvature_H[idx] = inner_oriented_H;
        result.inner_gaussian_curvature_K[idx] = inner_K;
        result.inner_graph_normal_dot_radial[idx] = inner_graph_normal_dot_radial;
        result.inner_outward_normal_alignment_flag[idx] = static_cast<uint8_t>(inner_graph_normal_is_outward ? 1 : 0);
        result.inner_orientation_flip_applied_flag[idx] = static_cast<uint8_t>(inner_graph_normal_is_outward ? 0 : 1);
        result.outer_K_confidence_class[idx] = 2;
        result.inner_K_confidence_class[idx] = 2;
    }

    std::vector<double> outer_valid_K;
    std::vector<double> inner_valid_K;
    outer_valid_K.reserve(result.curvature_valid_node_count);
    inner_valid_K.reserve(result.curvature_valid_node_count);
    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.curvature_valid_mask[idx] == 0) {
            continue;
        }
        outer_valid_K.push_back(result.outer_gaussian_curvature_K[idx]);
        inner_valid_K.push_back(result.inner_gaussian_curvature_K[idx]);
    }

    const double outer_global_median = medianOfValues(outer_valid_K);
    const double inner_global_median = medianOfValues(inner_valid_K);
    const double outer_global_sigma = robustSigmaFromValues(outer_valid_K, outer_global_median);
    const double inner_global_sigma = robustSigmaFromValues(inner_valid_K, inner_global_median);

    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.curvature_valid_mask[idx] == 0) {
            continue;
        }
        if (std::isfinite(outer_global_sigma) && outer_global_sigma > 0.0 &&
            std::abs(result.outer_gaussian_curvature_K[idx] - outer_global_median) > config.stage9_qc_n_tail * outer_global_sigma) {
            result.outer_K_global_tail_flag[idx] = 1;
            ++result.outer_K_global_tail_count;
        }
        if (std::isfinite(inner_global_sigma) && inner_global_sigma > 0.0 &&
            std::abs(result.inner_gaussian_curvature_K[idx] - inner_global_median) > config.stage9_qc_n_tail * inner_global_sigma) {
            result.inner_K_global_tail_flag[idx] = 1;
            ++result.inner_K_global_tail_count;
        }
    }

    for (std::size_t j = 0; j < result.grid.ny; ++j) {
        for (std::size_t i = 0; i < result.grid.nx; ++i) {
            const std::size_t idx = stage4NodeIndex(i, j, result.grid.nx);
            if (result.curvature_valid_mask[idx] == 0) {
                continue;
            }

            std::vector<double> outer_neighbor_K;
            std::vector<double> inner_neighbor_K;
            outer_neighbor_K.reserve(8);
            inner_neighbor_K.reserve(8);
            for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    if (di == 0 && dj == 0) {
                        continue;
                    }
                    const long ni = static_cast<long>(i) + di;
                    const long nj = static_cast<long>(j) + dj;
                    if (ni < 0 || nj < 0 || ni >= static_cast<long>(result.grid.nx) || nj >= static_cast<long>(result.grid.ny)) {
                        continue;
                    }
                    const std::size_t nidx =
                        stage4NodeIndex(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj), result.grid.nx);
                    if (result.curvature_valid_mask[nidx] == 0) {
                        continue;
                    }
                    outer_neighbor_K.push_back(result.outer_gaussian_curvature_K[nidx]);
                    inner_neighbor_K.push_back(result.inner_gaussian_curvature_K[nidx]);
                }
            }

            if (outer_neighbor_K.size() >= config.stage9_qc_min_neighbors) {
                const double local_median = medianOfValues(outer_neighbor_K);
                const double local_sigma = robustSigmaFromValues(outer_neighbor_K, local_median);
                const double local_scale = std::max(config.stage9_qc_abs_scale_floor, local_sigma);
                if (std::isfinite(local_scale) &&
                    std::abs(result.outer_gaussian_curvature_K[idx] - local_median) > config.stage9_qc_n_spike * local_scale) {
                    result.outer_K_local_spike_flag[idx] = 1;
                    ++result.outer_K_local_spike_count;
                }
            }
            if (inner_neighbor_K.size() >= config.stage9_qc_min_neighbors) {
                const double local_median = medianOfValues(inner_neighbor_K);
                const double local_sigma = robustSigmaFromValues(inner_neighbor_K, local_median);
                const double local_scale = std::max(config.stage9_qc_abs_scale_floor, local_sigma);
                if (std::isfinite(local_scale) &&
                    std::abs(result.inner_gaussian_curvature_K[idx] - local_median) > config.stage9_qc_n_spike * local_scale) {
                    result.inner_K_local_spike_flag[idx] = 1;
                    ++result.inner_K_local_spike_count;
                }
            }
        }
    }

    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.curvature_valid_mask[idx] == 0) {
            continue;
        }
        const bool outer_warn = result.outer_K_global_tail_flag[idx] != 0 || result.outer_K_local_spike_flag[idx] != 0;
        const bool inner_warn = result.inner_K_global_tail_flag[idx] != 0 || result.inner_K_local_spike_flag[idx] != 0;
        if (outer_warn) {
            result.outer_K_qc_warn_flag[idx] = 1;
            result.outer_K_confidence_class[idx] = 1;
            ++result.outer_K_qc_warn_count;
        }
        if (inner_warn) {
            result.inner_K_qc_warn_flag[idx] = 1;
            result.inner_K_confidence_class[idx] = 1;
            ++result.inner_K_qc_warn_count;
        }
    }

    if (result.curvature_valid_node_count > 0) {
        const double denom = static_cast<double>(result.curvature_valid_node_count);
        result.outer_K_global_tail_fraction_of_curvature_valid = static_cast<double>(result.outer_K_global_tail_count) / denom;
        result.inner_K_global_tail_fraction_of_curvature_valid = static_cast<double>(result.inner_K_global_tail_count) / denom;
        result.outer_K_local_spike_fraction_of_curvature_valid = static_cast<double>(result.outer_K_local_spike_count) / denom;
        result.inner_K_local_spike_fraction_of_curvature_valid = static_cast<double>(result.inner_K_local_spike_count) / denom;
        result.outer_K_qc_warn_fraction_of_curvature_valid = static_cast<double>(result.outer_K_qc_warn_count) / denom;
        result.inner_K_qc_warn_fraction_of_curvature_valid = static_cast<double>(result.inner_K_qc_warn_count) / denom;
    }

    const Stage9ScalarStats outer_H_stats = computeStage9ScalarStats(result.outer_mean_curvature_H, result.curvature_valid_mask);
    result.outer_mean_H = outer_H_stats.mean;
    result.outer_median_H = outer_H_stats.median;
    result.outer_stddev_H = outer_H_stats.stddev;
    result.outer_min_H = outer_H_stats.min;
    result.outer_max_H = outer_H_stats.max;
    const Stage9ScalarStats outer_oriented_H_stats =
        computeStage9ScalarStats(result.outer_oriented_mean_curvature_H, result.curvature_valid_mask);
    result.outer_mean_oriented_H = outer_oriented_H_stats.mean;
    result.outer_median_oriented_H = outer_oriented_H_stats.median;
    result.outer_stddev_oriented_H = outer_oriented_H_stats.stddev;
    result.outer_min_oriented_H = outer_oriented_H_stats.min;
    result.outer_max_oriented_H = outer_oriented_H_stats.max;

    const Stage9ScalarStats outer_K_stats = computeStage9ScalarStats(result.outer_gaussian_curvature_K, result.curvature_valid_mask);
    result.outer_mean_K = outer_K_stats.mean;
    result.outer_median_K = outer_K_stats.median;
    result.outer_stddev_K = outer_K_stats.stddev;
    result.outer_min_K = outer_K_stats.min;
    result.outer_max_K = outer_K_stats.max;

    const Stage9ScalarStats inner_H_stats = computeStage9ScalarStats(result.inner_mean_curvature_H, result.curvature_valid_mask);
    result.inner_mean_H = inner_H_stats.mean;
    result.inner_median_H = inner_H_stats.median;
    result.inner_stddev_H = inner_H_stats.stddev;
    result.inner_min_H = inner_H_stats.min;
    result.inner_max_H = inner_H_stats.max;
    const Stage9ScalarStats inner_oriented_H_stats =
        computeStage9ScalarStats(result.inner_oriented_mean_curvature_H, result.curvature_valid_mask);
    result.inner_mean_oriented_H = inner_oriented_H_stats.mean;
    result.inner_median_oriented_H = inner_oriented_H_stats.median;
    result.inner_stddev_oriented_H = inner_oriented_H_stats.stddev;
    result.inner_min_oriented_H = inner_oriented_H_stats.min;
    result.inner_max_oriented_H = inner_oriented_H_stats.max;

    const Stage9ScalarStats inner_K_stats = computeStage9ScalarStats(result.inner_gaussian_curvature_K, result.curvature_valid_mask);
    result.inner_mean_K = inner_K_stats.mean;
    result.inner_median_K = inner_K_stats.median;
    result.inner_stddev_K = inner_K_stats.stddev;
    result.inner_min_K = inner_K_stats.min;
    result.inner_max_K = inner_K_stats.max;

    if (result.metric_domain_node_count > 0) {
        result.curvature_valid_fraction_of_metric_domain =
            static_cast<double>(result.curvature_valid_node_count) / static_cast<double>(result.metric_domain_node_count);
    }
    if (result.derivative_valid_node_count > 0) {
        result.curvature_valid_fraction_of_derivative_valid =
            static_cast<double>(result.curvature_valid_node_count) / static_cast<double>(result.derivative_valid_node_count);
    }

    if (config.stage9_export_csv) {
        result.outer_curvature_csv_path = geometryArtifactPath(config, "_outer_curvature.csv");
        result.inner_curvature_csv_path = geometryArtifactPath(config, "_inner_curvature.csv");
        result.curvature_valid_mask_csv_path = geometryArtifactPath(config, "_curvature_valid_mask.csv");
        result.curvature_summary_csv_path = geometryArtifactPath(config, "_curvature_summary.csv");
        if (!writeStage9CurvatureCsv(result,
                                     result.outer_curvature_csv_path,
                                     result.outer_mean_curvature_H,
                                     result.outer_oriented_mean_curvature_H,
                                     result.outer_gaussian_curvature_K,
                                     result.outer_graph_normal_dot_radial,
                                     result.outer_orientation_flip_applied_flag,
                                     result.outer_K_global_tail_flag,
                                     result.outer_K_local_spike_flag,
                                     result.outer_K_qc_warn_flag,
                                     result.outer_K_confidence_class)) {
            throw std::runtime_error("Failed to write Stage 9 outer curvature CSV");
        }
        if (!writeStage9CurvatureCsv(result,
                                     result.inner_curvature_csv_path,
                                     result.inner_mean_curvature_H,
                                     result.inner_oriented_mean_curvature_H,
                                     result.inner_gaussian_curvature_K,
                                     result.inner_graph_normal_dot_radial,
                                     result.inner_orientation_flip_applied_flag,
                                     result.inner_K_global_tail_flag,
                                     result.inner_K_local_spike_flag,
                                     result.inner_K_qc_warn_flag,
                                     result.inner_K_confidence_class)) {
            throw std::runtime_error("Failed to write Stage 9 inner curvature CSV");
        }
        if (!writeStage9CurvatureMaskCsv(result, result.curvature_valid_mask_csv_path)) {
            throw std::runtime_error("Failed to write Stage 9 curvature-valid mask CSV");
        }
        if (!writeStage9SummaryCsv(result, result.curvature_summary_csv_path)) {
            throw std::runtime_error("Failed to write Stage 9 curvature summary CSV");
        }
    }

    (void)tolerance;
    result.messages.push_back("Geometry Stage 9 metric-domain node count: " + std::to_string(result.metric_domain_node_count));
    result.messages.push_back("Geometry Stage 9 derivative-valid node count: " + std::to_string(result.derivative_valid_node_count));
    result.messages.push_back("Geometry Stage 9 curvature-valid node count: " + std::to_string(result.curvature_valid_node_count));
    result.messages.push_back("Geometry Stage 9 invalid nonfinite-input count: " +
                              std::to_string(result.curvature_invalid_nonfinite_input_count));
    result.messages.push_back("Geometry Stage 9 invalid nonfinite-output count: " +
                              std::to_string(result.curvature_invalid_nonfinite_output_count));
    result.messages.push_back("Geometry Stage 9 outer K global-tail count: " + std::to_string(result.outer_K_global_tail_count));
    result.messages.push_back("Geometry Stage 9 outer K local-spike count: " + std::to_string(result.outer_K_local_spike_count));
    result.messages.push_back("Geometry Stage 9 outer K warn count: " + std::to_string(result.outer_K_qc_warn_count));
    result.messages.push_back("Geometry Stage 9 inner K global-tail count: " + std::to_string(result.inner_K_global_tail_count));
    result.messages.push_back("Geometry Stage 9 inner K local-spike count: " + std::to_string(result.inner_K_local_spike_count));
    result.messages.push_back("Geometry Stage 9 inner K warn count: " + std::to_string(result.inner_K_qc_warn_count));
    result.messages.push_back("Geometry Stage 9 curvature-valid fraction of metric domain: " +
                              std::to_string(result.curvature_valid_fraction_of_metric_domain));
    result.messages.push_back("Geometry Stage 9 curvature-valid fraction of derivative valid: " +
                              std::to_string(result.curvature_valid_fraction_of_derivative_valid));
    result.messages.push_back("Geometry Stage 9 outer H stats (mean/median/std/min/max): " + std::to_string(result.outer_mean_H) +
                              ", " + std::to_string(result.outer_median_H) + ", " + std::to_string(result.outer_stddev_H) + ", " +
                              std::to_string(result.outer_min_H) + ", " + std::to_string(result.outer_max_H));
    result.messages.push_back("Geometry Stage 9 outer oriented H stats (mean/median/std/min/max): " +
                              std::to_string(result.outer_mean_oriented_H) + ", " +
                              std::to_string(result.outer_median_oriented_H) + ", " +
                              std::to_string(result.outer_stddev_oriented_H) + ", " +
                              std::to_string(result.outer_min_oriented_H) + ", " +
                              std::to_string(result.outer_max_oriented_H));
    result.messages.push_back("Geometry Stage 9 outer K stats (mean/median/std/min/max): " + std::to_string(result.outer_mean_K) +
                              ", " + std::to_string(result.outer_median_K) + ", " + std::to_string(result.outer_stddev_K) + ", " +
                              std::to_string(result.outer_min_K) + ", " + std::to_string(result.outer_max_K));
    result.messages.push_back("Geometry Stage 9 inner H stats (mean/median/std/min/max): " + std::to_string(result.inner_mean_H) +
                              ", " + std::to_string(result.inner_median_H) + ", " + std::to_string(result.inner_stddev_H) + ", " +
                              std::to_string(result.inner_min_H) + ", " + std::to_string(result.inner_max_H));
    result.messages.push_back("Geometry Stage 9 inner oriented H stats (mean/median/std/min/max): " +
                              std::to_string(result.inner_mean_oriented_H) + ", " +
                              std::to_string(result.inner_median_oriented_H) + ", " +
                              std::to_string(result.inner_stddev_oriented_H) + ", " +
                              std::to_string(result.inner_min_oriented_H) + ", " +
                              std::to_string(result.inner_max_oriented_H));
    result.messages.push_back("Geometry Stage 9 inner K stats (mean/median/std/min/max): " + std::to_string(result.inner_mean_K) +
                              ", " + std::to_string(result.inner_median_K) + ", " + std::to_string(result.inner_stddev_K) + ", " +
                              std::to_string(result.inner_min_K) + ", " + std::to_string(result.inner_max_K));
    if (config.stage9_export_csv) {
        result.messages.push_back("Geometry Stage 9 outer curvature CSV: " + result.outer_curvature_csv_path);
        result.messages.push_back("Geometry Stage 9 inner curvature CSV: " + result.inner_curvature_csv_path);
        result.messages.push_back("Geometry Stage 9 curvature-valid mask CSV: " + result.curvature_valid_mask_csv_path);
        result.messages.push_back("Geometry Stage 9 summary CSV: " + result.curvature_summary_csv_path);
    }
    result.messages.push_back("Geometry analysis: completed Stage 9 curvature computation.");
    result.success = true;
    logMessages(result.messages, logger);
    return result;
}

GeometryStage10ThicknessResult runGeometryAnalysisStage10ThicknessComputation(
    const GeometryStage7SmoothedSurfaceResult& stage7_result,
    const GeometryStage8DerivativeEstimationResult& stage8_result,
    const GeometryStage9CurvatureComputationResult& stage9_result,
    const FoldPatchAnalysisConfig& config,
    Logger* logger,
    double tolerance) {
    GeometryStage10ThicknessResult result;
    if (!stage7_result.success) {
        throw std::runtime_error("Stage 10 cannot run before successful Stage 7 surface smoothing");
    }
    if (!stage9_result.success) {
        throw std::runtime_error("Stage 10 cannot run before successful Stage 9 curvature computation");
    }
    if (stage7_result.node_count == 0 || stage7_result.grid.nx == 0 || stage7_result.grid.ny == 0) {
        throw std::runtime_error("Stage 10 requires a non-empty Stage 7 grid");
    }
    result.messages.push_back("Geometry Stage 10");
    result.messages.push_back("Geometry analysis: starting Stage 10 thickness computation.");
    result.grid = stage7_result.grid;
    result.node_count = stage7_result.node_count;
    result.reconstructed_mask = stage7_result.reconstructed_mask;
    result.reliable_core_mask = stage7_result.reliable_core_mask;
    result.metric_domain_mask = stage7_result.metric_domain_mask;
    result.derivative_valid_mask = stage8_result.derivative_valid_mask;
    if (result.derivative_valid_mask.size() != result.node_count) {
        result.derivative_valid_mask.assign(result.node_count, 0);
    }
    result.curvature_valid_mask = stage9_result.curvature_valid_mask;
    if (result.curvature_valid_mask.size() != result.node_count) {
        result.curvature_valid_mask.assign(result.node_count, 0);
    }

    result.thickness_domain_definition = config.stage10_use_curvature_valid_domain_only
                                             ? "stage7_metric_domain_intersect_stage9_curvature_valid"
                                             : "stage7_metric_domain";

    const double nan = std::numeric_limits<double>::quiet_NaN();
    result.thickness_vertical.assign(result.node_count, nan);
    result.thickness_outer_z_stage7.assign(result.node_count, nan);
    result.thickness_radial.assign(result.node_count, nan);
    result.thickness_radial_s_outer.assign(result.node_count, nan);
    result.thickness_radial_s_inner.assign(result.node_count, nan);
    result.thickness_attempted_mask.assign(result.node_count, 0);
    result.thickness_valid_mask.assign(result.node_count, 0);
    result.thickness_invalid_outside_domain_mask.assign(result.node_count, 0);
    result.thickness_invalid_nonfinite_surface_mask.assign(result.node_count, 0);
    result.thickness_invalid_negative_or_zero_mask.assign(result.node_count, 0);
    result.thickness_invalid_below_min_threshold_mask.assign(result.node_count, 0);
    result.thickness_invalid_above_max_threshold_mask.assign(result.node_count, 0);
    result.thickness_qc_warn_mask.assign(result.node_count, 0);
    result.thickness_radial_valid_mask.assign(result.node_count, 0);
    result.thickness_radial_invalid_no_bracket_mask.assign(result.node_count, 0);
    result.thickness_radial_invalid_root_failure_mask.assign(result.node_count, 0);
    result.thickness_radial_invalid_outside_inner_domain_mask.assign(result.node_count, 0);

    if (!stage8_result.success) {
        result.messages.push_back("Geometry Stage 10 note: Stage 8 derivative estimation unavailable; continuing with Stage 7 surfaces.");
    }

    std::size_t metric_domain_node_count = 0;
    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        if (result.metric_domain_mask[idx] != 0) {
            ++metric_domain_node_count;
        }
    }

    for (std::size_t idx = 0; idx < result.node_count; ++idx) {
        const bool in_domain = result.metric_domain_mask[idx] != 0 &&
                               (!config.stage10_use_curvature_valid_domain_only || result.curvature_valid_mask[idx] != 0);
        if (!in_domain) {
            result.thickness_invalid_outside_domain_mask[idx] = 1;
            ++result.thickness_excluded_outside_domain_count;
            continue;
        }

        result.thickness_attempted_mask[idx] = 1;
        ++result.thickness_attempted_node_count;

        const double z_outer = stage7_result.z_outer_smooth[idx];
        const double z_inner = stage7_result.z_inner_smooth[idx];
        if (!std::isfinite(z_outer) || !std::isfinite(z_inner)) {
            result.thickness_invalid_nonfinite_surface_mask[idx] = 1;
            ++result.thickness_attempted_invalid_node_count;
            ++result.thickness_invalid_nonfinite_surface_count;
            continue;
        }

        const double t_vertical = z_outer - z_inner;
        result.thickness_vertical[idx] = t_vertical;
        result.thickness_outer_z_stage7[idx] = z_outer;

        double t_selected = t_vertical;
        if (config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial) {
            const std::size_t i = idx % result.grid.nx;
            const std::size_t j = idx / result.grid.nx;
            const double x = result.grid.x_values[i];
            const double y = result.grid.y_values[j];
            const double s_out = std::sqrt(x * x + y * y + z_outer * z_outer);
            result.thickness_radial_s_outer[idx] = s_out;
            if (s_out <= tolerance) {
                result.thickness_radial_invalid_root_failure_mask[idx] = 1;
                ++result.thickness_attempted_invalid_node_count;
                ++result.thickness_radial_invalid_root_failure_count;
                continue;
            }
            const double ux = x / s_out;
            const double uy = y / s_out;
            const double uz = z_outer / s_out;
            const Stage10RadialIntersectionResult intersection = findRadialInnerIntersectionOnStage7Surface(
                result.grid, stage7_result.z_inner_smooth, result.metric_domain_mask, ux, uy, uz, s_out, tolerance);
            ++result.thickness_radial_attempted_node_count;
            if (intersection.status != Stage10RadialIntersectionStatus::success || !std::isfinite(intersection.s_inner)) {
                ++result.thickness_attempted_invalid_node_count;
                if (intersection.status == Stage10RadialIntersectionStatus::outside_inner_sampling_domain) {
                    result.thickness_radial_invalid_outside_inner_domain_mask[idx] = 1;
                    ++result.thickness_radial_invalid_outside_inner_domain_count;
                } else if (intersection.status == Stage10RadialIntersectionStatus::no_bracket) {
                    result.thickness_radial_invalid_no_bracket_mask[idx] = 1;
                    ++result.thickness_radial_invalid_no_bracket_count;
                } else {
                    result.thickness_radial_invalid_root_failure_mask[idx] = 1;
                    ++result.thickness_radial_invalid_root_failure_count;
                }
                continue;
            }
            result.thickness_radial_s_inner[idx] = intersection.s_inner;
            const double t_radial = s_out - intersection.s_inner;
            result.thickness_radial[idx] = t_radial;
            t_selected = t_radial;
        }

        if (t_selected <= tolerance) {
            result.thickness_invalid_negative_or_zero_mask[idx] = 1;
            ++result.thickness_attempted_invalid_node_count;
            ++result.thickness_invalid_negative_or_zero_count;
            continue;
        }
        if (config.stage10_min_thickness > 0.0 && t_selected < config.stage10_min_thickness) {
            result.thickness_invalid_below_min_threshold_mask[idx] = 1;
            ++result.thickness_attempted_invalid_node_count;
            ++result.thickness_invalid_below_min_threshold_count;
            continue;
        }
        if (config.stage10_max_thickness > 0.0 && t_selected > config.stage10_max_thickness) {
            result.thickness_invalid_above_max_threshold_mask[idx] = 1;
            ++result.thickness_attempted_invalid_node_count;
            ++result.thickness_invalid_above_max_threshold_count;
            continue;
        }

        result.thickness_valid_mask[idx] = 1;
        if (config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial) {
            result.thickness_radial_valid_mask[idx] = 1;
            ++result.thickness_radial_valid_node_count;
        }
        ++result.thickness_valid_node_count;
        if (result.curvature_valid_mask[idx] != 0) {
            ++result.thickness_valid_and_curvature_valid_node_count;
        }

        if (!config.stage10_use_curvature_valid_domain_only && result.curvature_valid_mask[idx] == 0) {
            result.thickness_qc_warn_mask[idx] = 1;
            ++result.thickness_qc_warn_count;
        }
    }

    const Stage10ScalarStats vertical_stats = computeStage10ScalarStats(result.thickness_vertical, result.thickness_valid_mask);
    result.mean_thickness_vertical = vertical_stats.mean;
    result.median_thickness_vertical = vertical_stats.median;
    result.stddev_thickness_vertical = vertical_stats.stddev;
    result.min_thickness_vertical = vertical_stats.min;
    result.max_thickness_vertical = vertical_stats.max;
    const Stage10ScalarStats radial_stats = computeStage10ScalarStats(result.thickness_radial, result.thickness_valid_mask);
    result.mean_thickness_radial = radial_stats.mean;
    result.median_thickness_radial = radial_stats.median;
    result.stddev_thickness_radial = radial_stats.stddev;
    result.min_thickness_radial = radial_stats.min;
    result.max_thickness_radial = radial_stats.max;

    result.stage10_method = stage10ThicknessMethodLabel(config.stage10_thickness_method);
    if (config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial) {
        result.thickness_method_label = "stage10_radial_origin_centered_ray_intersection";
        result.local_thickness_definition =
            "distance_between_stage7_outer_point_and_stage7_inner_intersection_along_origin_centered_ray";
        result.thickness_contract_note =
            "stage10_radial_thickness_is_computed_by_tracing_the_origin_to_outer_surface_ray_and_intersecting_it_with_the_stage7_inner_surface";
        result.thickness_domain_limitation_note =
            "stage10_radial_validity_requires_not_only_outer_node_membership_in_stage7_metric_domain_but_also_inner_surface_sampleability_along_the_origin_centered_radial_ray";
        result.radial_invalid_outside_inner_domain_definition =
            "outer_node_is_in_stage7_metric_domain_but_required_inner_surface_evaluation_points_along_the_radial_ray_leave_the_current_valid_inner_sampling_or_interpolation_domain";
        result.mean_thickness_active = result.mean_thickness_radial;
        result.median_thickness_active = result.median_thickness_radial;
        result.stddev_thickness_active = result.stddev_thickness_radial;
        result.min_thickness_active = result.min_thickness_radial;
        result.max_thickness_active = result.max_thickness_radial;
        result.messages.push_back(
            "Stage10 radial P_out/P_in export is only defined for nodes where the inner surface remains sampleable along the origin-centered radial ray.");
    } else {
        result.mean_thickness_active = result.mean_thickness_vertical;
        result.median_thickness_active = result.median_thickness_vertical;
        result.stddev_thickness_active = result.stddev_thickness_vertical;
        result.min_thickness_active = result.min_thickness_vertical;
        result.max_thickness_active = result.max_thickness_vertical;
    }

    if (metric_domain_node_count > 0) {
        result.thickness_valid_fraction_of_metric_domain =
            static_cast<double>(result.thickness_valid_node_count) / static_cast<double>(metric_domain_node_count);
        result.thickness_valid_intersection_fraction_of_metric_domain =
            static_cast<double>(result.thickness_valid_and_curvature_valid_node_count) /
            static_cast<double>(metric_domain_node_count);
        if (config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial) {
            result.thickness_radial_valid_fraction_of_metric_domain =
                static_cast<double>(result.thickness_radial_valid_node_count) / static_cast<double>(metric_domain_node_count);
            result.thickness_radial_invalid_outside_inner_domain_fraction_of_metric_domain =
                static_cast<double>(result.thickness_radial_invalid_outside_inner_domain_count) /
                static_cast<double>(metric_domain_node_count);
        }
    }

    if (config.stage10_export_csv) {
        const std::string method_suffix =
            config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial ? "radial" : "vertical";
        result.thickness_vertical_csv_path = geometryArtifactPath(config, "_thickness_") + method_suffix + ".csv";
        result.thickness_valid_mask_csv_path = geometryArtifactPath(config, "_thickness_") + method_suffix + "_valid_mask.csv";
        result.thickness_invalid_reason_csv_path = geometryArtifactPath(config, "_thickness_") + method_suffix + "_invalid_reason.csv";
        result.thickness_summary_csv_path = geometryArtifactPath(config, "_thickness_") + method_suffix + "_summary.csv";
        if (!writeStage10ThicknessCsv(result, result.thickness_vertical_csv_path)) {
            throw std::runtime_error("Failed to write Stage 10 thickness CSV");
        }
        if (!writeStage10ValidMaskCsv(result, result.thickness_valid_mask_csv_path)) {
            throw std::runtime_error("Failed to write Stage 10 thickness valid-mask CSV");
        }
        if (!writeStage10InvalidReasonCsv(result, result.thickness_invalid_reason_csv_path)) {
            throw std::runtime_error("Failed to write Stage 10 thickness invalid-reason CSV");
        }
        if (!writeStage10SummaryCsv(result, result.thickness_summary_csv_path)) {
            throw std::runtime_error("Failed to write Stage 10 thickness summary CSV");
        }
        if (config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial) {
            result.thickness_radial_p_out_in_csv_path = geometryArtifactPath(config, "_thickness_radial_P_out_P_in_t.csv");
            if (!writeStage10RadialPOutInCsv(result, result.thickness_radial_p_out_in_csv_path)) {
                throw std::runtime_error("Failed to write Stage 10 radial P_out/P_in CSV");
            }
        }
    }

    result.messages.push_back("Geometry Stage 10 method: " + result.stage10_method);
    result.messages.push_back("Geometry Stage 10 domain mode: " + result.thickness_domain_definition);
    result.messages.push_back("Geometry Stage 10 attempted node count: " + std::to_string(result.thickness_attempted_node_count));
    result.messages.push_back("Geometry Stage 10 valid node count: " + std::to_string(result.thickness_valid_node_count));
    result.messages.push_back("Geometry Stage 10 excluded outside-domain node count: " +
                              std::to_string(result.thickness_excluded_outside_domain_count));
    result.messages.push_back("Geometry Stage 10 attempted-invalid node count: " +
                              std::to_string(result.thickness_attempted_invalid_node_count));
    result.messages.push_back("Geometry Stage 10 thickness-valid and curvature-valid intersection node count: " +
                              std::to_string(result.thickness_valid_and_curvature_valid_node_count));

    if (result.thickness_valid_node_count > 0) {
        if (config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial) {
            result.messages.push_back("Geometry Stage 10 radial attempted node count: " +
                                      std::to_string(result.thickness_radial_attempted_node_count));
            result.messages.push_back("Geometry Stage 10 radial valid node count: " +
                                      std::to_string(result.thickness_radial_valid_node_count));
            result.messages.push_back("Geometry Stage 10 radial invalid no-bracket count: " +
                                      std::to_string(result.thickness_radial_invalid_no_bracket_count));
            result.messages.push_back("Geometry Stage 10 radial invalid root-failure count: " +
                                      std::to_string(result.thickness_radial_invalid_root_failure_count));
            result.messages.push_back("Geometry Stage 10 radial invalid outside-inner-domain count: " +
                                      std::to_string(result.thickness_radial_invalid_outside_inner_domain_count));
            result.messages.push_back("Geometry Stage 10 radial active thickness stats (mean/min/max): " +
                                      std::to_string(result.mean_thickness_active) + ", " +
                                      std::to_string(result.min_thickness_active) + ", " +
                                      std::to_string(result.max_thickness_active));
        } else {
            result.messages.push_back("Geometry Stage 10 active thickness stats (mean/min/max): " +
                                      std::to_string(result.mean_thickness_active) + ", " +
                                      std::to_string(result.min_thickness_active) + ", " +
                                      std::to_string(result.max_thickness_active));
        }
        result.success = true;
    } else {
        result.messages.push_back("Geometry Stage 10 failure: no valid thickness nodes were produced.");
        result.success = false;
    }

    if (config.stage10_export_csv) {
        result.messages.push_back("Geometry Stage 10 thickness CSV: " + result.thickness_vertical_csv_path);
        result.messages.push_back("Geometry Stage 10 valid-mask CSV: " + result.thickness_valid_mask_csv_path);
        result.messages.push_back("Geometry Stage 10 invalid-reason CSV: " + result.thickness_invalid_reason_csv_path);
        result.messages.push_back("Geometry Stage 10 summary CSV: " + result.thickness_summary_csv_path);
        if (config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial) {
            result.messages.push_back("Geometry Stage 10 radial P_out/P_in CSV: " + result.thickness_radial_p_out_in_csv_path);
        }
    }
    result.messages.push_back("Geometry analysis: completed Stage 10 thickness computation.");
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
    result.stage3_patch = runGeometryAnalysisStage3PatchNormalization(result.stage2_patch, config, logger);
    result.stage4_raw =
        runGeometryAnalysisStage4RawSheetDetection(capsid, config, parser_config, result.stage3_patch, logger);
    result.stage5_prep = runGeometryAnalysisStage5SurfacePreparation(result.stage4_raw, config, logger);
    result.stage6_surfaces = runGeometryAnalysisStage6SurfaceReconstruction(result.stage5_prep, config, logger);
    if (config.stage7_enabled) {
        result.stage7_smooth =
            runGeometryAnalysisStage7SurfaceSmoothing(result.stage6_surfaces, result.stage5_prep, config, logger);
    } else {
        result.stage7_smooth.success = true;
        result.stage7_smooth.messages.push_back("Geometry Stage 7 smoothing disabled by configuration.");
    }
    const bool stage8_allowed_by_stage7_method =
        config.stage7_method == FoldPatchAnalysisConfig::Stage7Method::thin_plate_grid_fit;
    if (config.stage8_enabled && stage8_allowed_by_stage7_method) {
        result.stage8_derivatives =
            runGeometryAnalysisStage8DerivativeEstimation(result.stage7_smooth, result.stage5_prep, config, logger);
    } else {
        result.stage8_derivatives.success = true;
        if (!config.stage8_enabled) {
            result.stage8_derivatives.messages.push_back("Geometry Stage 8 derivative estimation disabled by configuration.");
        } else {
            result.stage8_derivatives.messages.push_back(
                "Geometry Stage 8 derivative estimation skipped: requires Stage 7 thin_plate_grid_fit method.");
        }
        result.messages.insert(result.messages.end(),
                               result.stage8_derivatives.messages.begin(),
                               result.stage8_derivatives.messages.end());
    }
    if (config.stage9_enabled && config.stage8_enabled && stage8_allowed_by_stage7_method) {
        result.stage9_curvature =
            runGeometryAnalysisStage9CurvatureComputation(result.stage7_smooth, result.stage8_derivatives, config, logger);
    } else {
        result.stage9_curvature.success = true;
        if (!config.stage9_enabled) {
            result.stage9_curvature.messages.push_back("Geometry Stage 9 curvature computation disabled by configuration.");
        } else if (!config.stage8_enabled) {
            result.stage9_curvature.messages.push_back(
                "Geometry Stage 9 curvature computation skipped: requires Stage 8 derivative estimation.");
        } else {
            result.stage9_curvature.messages.push_back(
                "Geometry Stage 9 curvature computation skipped: requires Stage 7 thin_plate_grid_fit method.");
        }
        result.messages.insert(result.messages.end(),
                               result.stage9_curvature.messages.begin(),
                               result.stage9_curvature.messages.end());
    }
    if (config.stage10_enabled) {
        result.stage10_thickness = runGeometryAnalysisStage10ThicknessComputation(
            result.stage7_smooth, result.stage8_derivatives, result.stage9_curvature, config, logger);
    } else {
        result.stage10_thickness.messages.push_back("Geometry Stage 10 thickness computation disabled by configuration.");
        result.messages.insert(result.messages.end(),
                               result.stage10_thickness.messages.begin(),
                               result.stage10_thickness.messages.end());
    }
    result.success = result.preparation.success && result.stage2_patch.success && result.stage3_patch.success &&
                     result.stage4_raw.success && result.stage5_prep.success && result.stage6_surfaces.success &&
                     result.stage7_smooth.success && result.stage8_derivatives.success && result.stage9_curvature.success &&
                     (!config.stage10_enabled || result.stage10_thickness.success);
    const std::string run_summary_json_path = geometryArtifactPath(config, "stage_run_summary.json");
    if (!writeGeometryRunSummaryJson(config, parser_config, run_summary_json_path)) {
        throw std::runtime_error("Failed to write geometry run summary JSON: " + run_summary_json_path);
    }
    result.run_summary_json_path = run_summary_json_path;
    logMessages(result.messages, logger);

    return result;
}

std::string buildGeometrySummaryReport(const GeometryAnalysisResult& result,
                                       const FoldPatchAnalysisConfig& config,
                                       const std::string& input_name) {
    const auto fmtCount = [](std::size_t value) -> std::string { return std::to_string(value); };
    const auto fmtFixed = [](double value, int precision = 3) -> std::string {
        if (!std::isfinite(value)) {
            return "n/a";
        }
        std::ostringstream out;
        out << std::fixed << std::setprecision(precision) << value;
        return out.str();
    };
    const auto fmtFraction = [](std::size_t numerator, std::size_t denominator) -> std::string {
        if (denominator == 0) {
            return "n/a";
        }
        std::ostringstream out;
        out << std::fixed << std::setprecision(3) << (static_cast<double>(numerator) / static_cast<double>(denominator));
        return out.str();
    };
    const auto appendKeyValue = [](std::ostringstream& out, const std::string& key, const std::string& value) {
        out << std::left << std::setw(18) << (key + ":") << value << '\n';
    };

    const std::size_t patch_atoms = result.stage3_patch.analytical_patch.atom_count > 0
                                        ? result.stage3_patch.analytical_patch.atom_count
                                        : result.stage2_patch.selected_atoms_count;
    const std::size_t inside_disk_count = result.stage4_raw.inside_disk_count;
    const std::size_t metric_domain_count = result.stage7_smooth.metric_domain_node_count;
    const std::size_t derivative_valid_count = result.stage8_derivatives.derivative_valid_node_count;
    const std::size_t curvature_valid_count = result.stage9_curvature.curvature_valid_node_count;
    const std::size_t thickness_valid_count = result.stage10_thickness.thickness_valid_node_count;

    std::vector<std::string> notes;
    notes.reserve(5);
    const bool radial_method = config.stage10_thickness_method == FoldPatchAnalysisConfig::Stage10ThicknessMethod::radial;
    notes.push_back(radial_method
                        ? "Thickness is reported using the radial method (origin-centered ray from outer to inner surface), which is generally more physically meaningful for capsid shell geometry over larger patches."
                        : "Thickness is reported using the vertical graph-gap method (outer z minus inner z); this is most reliable for local patches near the aligned fold axis.");

    const double thickness_frac = result.stage10_thickness.thickness_valid_fraction_of_metric_domain;
    const double curvature_frac = result.stage9_curvature.curvature_valid_fraction_of_metric_domain;
    const double max_k_warn_frac =
        std::max(result.stage9_curvature.outer_K_qc_warn_fraction_of_curvature_valid,
                 result.stage9_curvature.inner_K_qc_warn_fraction_of_curvature_valid);

    if (radial_method && result.stage10_thickness.thickness_radial_invalid_outside_inner_domain_count > 0 &&
        notes.size() < 3) {
        notes.push_back(
            "Some radial-thickness nodes were excluded because the inner surface could not be sampled along the full radial ray near the patch edge.");
    }
    if (radial_method && result.stage10_thickness.thickness_radial_invalid_outside_inner_domain_count > 0 &&
        result.stage10_thickness.thickness_radial_invalid_no_bracket_count == 0 &&
        result.stage10_thickness.thickness_radial_invalid_root_failure_count == 0 && notes.size() < 5) {
        notes.push_back(
            "Radial-thickness failures are currently dominated by inner-surface sampling-domain loss along the radial ray, not by bracketing or root-convergence failure.");
    }
    if (result.stage4_raw.unique_both_contact_atom_count > 0 && notes.size() < 5) {
        notes.push_back(
            "Warning. Raw sheet-contact overlap indicators were detected in Stage 4; inspect raw-contact diagnostics if needed.");
    }
    if (notes.size() < 3) {
        if (std::isfinite(thickness_frac) && thickness_frac >= 0.98) {
            notes.push_back("Thickness coverage is effectively complete over the metric domain.");
        } else if (std::isfinite(thickness_frac) && thickness_frac >= 0.90) {
            notes.push_back("Thickness coverage is high, but a small fraction of the metric domain was excluded.");
        } else if (std::isfinite(thickness_frac)) {
            notes.push_back(
                "Thickness coverage is reduced; interpret full-domain thickness summaries cautiously and prefer the trusted core.");
        }
    }
    if (notes.size() < 3) {
        if (std::isfinite(curvature_frac) && curvature_frac >= 0.95) {
            notes.push_back("Curvature coverage is high over the metric domain.");
        } else if (std::isfinite(curvature_frac) && curvature_frac >= 0.85) {
            notes.push_back("Curvature coverage is good but not complete; edge-region fits are beginning to drop out.");
        } else if (std::isfinite(curvature_frac)) {
            notes.push_back(
                "Curvature coverage is limited; curvature summaries increasingly reflect the trusted interior rather than the full patch.");
        }
    }
    if (notes.size() < 3 && std::isfinite(max_k_warn_frac)) {
        if (max_k_warn_frac < 0.05) {
            notes.push_back("Gaussian curvature is in a low-warning regime; K summaries are comparatively stable.");
        } else if (max_k_warn_frac < 0.15) {
            notes.push_back("Gaussian curvature shows a moderate warning burden; interpret K tails and local extremes with care.");
        } else {
            notes.push_back(
                "Gaussian curvature is in a high-warning regime; large-|K| values are likely the least robust curvature quantities in this run.");
        }
    }
    if (notes.size() < 3) {
        if (config.cylinder_radius > 60.0) {
            notes.push_back(
                "This is a large patch relative to the fold axis; whole-domain averages combine multiple local geometric regimes rather than a single local cap.");
        } else if (config.cylinder_radius > 30.0) {
            notes.push_back(
                "This patch is no longer purely local; curvature and thickness summarize a broader fold-centered region.");
        }
    }
    if (notes.size() < 5) {
        notes.push_back(
            "Mean curvature H is reported with the surface-orientation sign convention used internally; its magnitude is usually more robust than its sign when summarizing broad capsid patches.");
    }
    notes.push_back(
        "Gaussian curvature K distinguishes local shape class: positive K is dome-like, negative K is saddle-like, and values near zero are weakly curved or cylindrical-like.");
    while (notes.size() < 2) {
        static const char* kFallbacks[] = {
            "Detailed per-node values and masks are available in the exported CSV files.",
            "Summary statistics are reported over the valid analysis domain, not over excluded nodes.",
        };
        notes.push_back(kFallbacks[notes.size()]);
    }

    std::ostringstream out;
    out << "==================== Geometry Summary ====================\n";
    appendKeyValue(out, "Input", input_name);
    appendKeyValue(out,
                   "Fold",
                   result.preparation.resolved_fold_name + " (type=" + std::to_string(config.fold_type) +
                       ", index=" + std::to_string(config.fold_index) + ")");
    appendKeyValue(out, "Working frame", result.preparation.final_frame_description);
    appendKeyValue(out, "Cylinder radius", fmtFixed(config.cylinder_radius) + " \u00C5");
    appendKeyValue(out, "Grid spacing", fmtFixed(config.grid_spacing) + " \u00C5");
    appendKeyValue(out, "Thickness method", radial_method ? "radial" : "vertical");
    appendKeyValue(out, "Patch atoms", fmtCount(patch_atoms));
    appendKeyValue(out, "Grid", fmtCount(result.stage4_raw.grid.nx) + " x " + fmtCount(result.stage4_raw.grid.ny));
    appendKeyValue(out, "Metric domain", fmtCount(metric_domain_count) + " nodes");
    appendKeyValue(out, "Curvature-valid", fmtCount(curvature_valid_count) + " nodes");
    appendKeyValue(out,
                   "Run summary JSON",
                   result.run_summary_json_path.empty() ? "n/a" : result.run_summary_json_path);
    appendKeyValue(out, "Output root", config.output_root_dir);
    out << '\n';

    out << "Coverage / trust\n";
    out << "---------------------------------------------------------\n";
    out << std::left << std::setw(38) << "Inside-disk nodes" << fmtCount(inside_disk_count) << '\n';
    out << std::left << std::setw(38) << "Metric-domain fraction" << fmtFraction(metric_domain_count, inside_disk_count) << '\n';
    out << std::left << std::setw(38) << "Derivative-valid / metric-domain"
        << fmtFraction(derivative_valid_count, metric_domain_count) << '\n';
    out << std::left << std::setw(38) << "Curvature-valid / metric-domain"
        << (std::isfinite(curvature_frac) ? fmtFixed(curvature_frac) : fmtFraction(curvature_valid_count, metric_domain_count))
        << '\n';
    out << std::left << std::setw(38) << "Thickness-valid / metric-domain"
        << fmtFraction(thickness_valid_count, metric_domain_count) << '\n';
    out << std::left << std::setw(38) << "Outer K warn / curvature-valid"
        << fmtFraction(result.stage9_curvature.outer_K_qc_warn_count, curvature_valid_count) << '\n';
    out << std::left << std::setw(38) << "Inner K warn / curvature-valid"
        << fmtFraction(result.stage9_curvature.inner_K_qc_warn_count, curvature_valid_count) << '\n';
    if (radial_method) {
        out << std::left << std::setw(38) << "Radial outside-inner-domain frac"
            << fmtFraction(result.stage10_thickness.thickness_radial_invalid_outside_inner_domain_count, metric_domain_count)
            << '\n';
    }
    out << '\n';

    out << "Thickness\n";
    out << "---------------------------------------------------------\n";
    out << std::left << std::setw(11) << "Method" << std::setw(9) << "Valid" << std::setw(10) << "Mean" << std::setw(10)
        << "Median" << std::setw(10) << "StdDev" << std::setw(10) << "Min" << "Max\n";
    out << std::left << std::setw(11) << (radial_method ? "radial" : "vertical") << std::setw(9)
        << fmtCount(thickness_valid_count) << std::setw(10) << fmtFixed(result.stage10_thickness.mean_thickness_active)
        << std::setw(10) << fmtFixed(result.stage10_thickness.median_thickness_active) << std::setw(10)
        << fmtFixed(result.stage10_thickness.stddev_thickness_active) << std::setw(10)
        << fmtFixed(result.stage10_thickness.min_thickness_active) << fmtFixed(result.stage10_thickness.max_thickness_active)
        << '\n';
    out << '\n';

    out << "Curvature - outer surface\n";
    out << "---------------------------------------------------------\n";
    out << std::left << std::setw(8) << "Metric" << std::setw(10) << "Mean" << std::setw(10) << "Median" << std::setw(10)
        << "StdDev" << std::setw(10) << "Min" << "Max\n";
    out << std::left << std::setw(8) << "H" << std::setw(10) << fmtFixed(result.stage9_curvature.outer_mean_oriented_H)
        << std::setw(10) << fmtFixed(result.stage9_curvature.outer_median_oriented_H) << std::setw(10)
        << fmtFixed(result.stage9_curvature.outer_stddev_oriented_H) << std::setw(10)
        << fmtFixed(result.stage9_curvature.outer_min_oriented_H) << fmtFixed(result.stage9_curvature.outer_max_oriented_H)
        << '\n';
    out << std::left << std::setw(8) << "K" << std::setw(10) << fmtFixed(result.stage9_curvature.outer_mean_K)
        << std::setw(10) << fmtFixed(result.stage9_curvature.outer_median_K) << std::setw(10)
        << fmtFixed(result.stage9_curvature.outer_stddev_K) << std::setw(10) << fmtFixed(result.stage9_curvature.outer_min_K)
        << fmtFixed(result.stage9_curvature.outer_max_K) << '\n';
    out << "K QC warn fraction: " << fmtFraction(result.stage9_curvature.outer_K_qc_warn_count, curvature_valid_count)
        << "\n\n";

    out << "Curvature - inner surface\n";
    out << "---------------------------------------------------------\n";
    out << std::left << std::setw(8) << "Metric" << std::setw(10) << "Mean" << std::setw(10) << "Median" << std::setw(10)
        << "StdDev" << std::setw(10) << "Min" << "Max\n";
    out << std::left << std::setw(8) << "H" << std::setw(10) << fmtFixed(result.stage9_curvature.inner_mean_oriented_H)
        << std::setw(10) << fmtFixed(result.stage9_curvature.inner_median_oriented_H) << std::setw(10)
        << fmtFixed(result.stage9_curvature.inner_stddev_oriented_H) << std::setw(10)
        << fmtFixed(result.stage9_curvature.inner_min_oriented_H) << fmtFixed(result.stage9_curvature.inner_max_oriented_H)
        << '\n';
    out << std::left << std::setw(8) << "K" << std::setw(10) << fmtFixed(result.stage9_curvature.inner_mean_K)
        << std::setw(10) << fmtFixed(result.stage9_curvature.inner_median_K) << std::setw(10)
        << fmtFixed(result.stage9_curvature.inner_stddev_K) << std::setw(10) << fmtFixed(result.stage9_curvature.inner_min_K)
        << fmtFixed(result.stage9_curvature.inner_max_K) << '\n';
    out << "K QC warn fraction: " << fmtFraction(result.stage9_curvature.inner_K_qc_warn_count, curvature_valid_count)
        << "\n\n";

    out << "Notes\n";
    out << "---------------------------------------------------------\n";
    for (const std::string& note : notes) {
        out << "- " << note << '\n';
    }
    out << "==========================================================\n";
    return out.str();
}
