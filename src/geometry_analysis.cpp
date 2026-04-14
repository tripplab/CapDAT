#include "geometry_analysis.hpp"

#include "export_capsid.hpp"
#include "timer.hpp"

#include <cmath>
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
    out << "    \"output_prefix\": \"" << jsonEscape(config.output_prefix) << "\"\n";
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
        << result.max_smooth_separation << ',' << result.mean_smooth_separation << ',' << result.normal_orientation_outer
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
        result.stage3_normalized_atoms_csv_path = config.output_prefix + "_stage3_normalized_atoms.csv";
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
        result.outer_csv_path = config.output_prefix + "_outer_raw.csv";
        result.inner_csv_path = config.output_prefix + "_inner_raw.csv";
        result.valid_mask_csv_path = config.output_prefix + "_valid_mask_raw.csv";
        result.outer_only_mask_csv_path = config.output_prefix + "_outer_only_mask_raw.csv";
        result.inner_only_mask_csv_path = config.output_prefix + "_inner_only_mask_raw.csv";
        result.negative_thickness_mask_csv_path = config.output_prefix + "_negative_thickness_mask_raw.csv";
        result.summary_csv_path = config.output_prefix + "_stage4_summary.csv";
    }
    result.contact_atoms_pdb_path = config.output_prefix + "_raw_contact_atoms.pdb";

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
    if (result.unique_both_contact_atom_count == 0) {
        result.messages.push_back("Geometry Stage 4 unique both-contact patch atoms (index-based): " +
                                  std::to_string(result.unique_both_contact_atom_count));
    } else {
        result.messages.push_back(std::string(kWarningMessagePrefix) +
                                  "Geometry Stage 4 unique both-contact patch atoms (index-based): " +
                                  std::to_string(result.unique_both_contact_atom_count));
    }
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
    if (result.zero_thickness_node_count > 0 || result.negative_thickness_node_count > 0 ||
        result.unique_both_contact_atom_count > 0) {
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
        result.outer_seed_csv_path = config.output_prefix + "_outer_seed.csv";
        result.inner_seed_csv_path = config.output_prefix + "_inner_seed.csv";
        result.paired_seed_mask_csv_path = config.output_prefix + "_paired_seed_mask.csv";
        result.boundary_exclusion_mask_csv_path = config.output_prefix + "_boundary_exclusion_mask.csv";
        result.interp_allowed_mask_csv_path = config.output_prefix + "_interp_allowed_mask.csv";
        result.hard_invalid_mask_csv_path = config.output_prefix + "_hard_invalid_mask.csv";
        result.reliable_core_mask_csv_path = config.output_prefix + "_reliable_core_mask.csv";
        result.summary_csv_path = config.output_prefix + "_stage5_summary.csv";

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
        result.outer_reconstructed_csv_path = config.output_prefix + "_outer_reconstructed.csv";
        result.inner_reconstructed_csv_path = config.output_prefix + "_inner_reconstructed.csv";
        result.reconstructed_mask_csv_path = config.output_prefix + "_reconstructed_mask.csv";
        result.final_valid_analysis_mask_csv_path = config.output_prefix + "_final_valid_analysis_mask.csv";
        result.non_crossing_adjustment_mask_csv_path = config.output_prefix + "_non_crossing_adjustment_mask.csv";
        result.summary_csv_path = config.output_prefix + "_stage6_summary.csv";

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
            result.outer_obj_path = config.output_prefix + "_outer_surface" + extension;
            result.inner_obj_path = config.output_prefix + "_inner_surface" + extension;
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
            result.outer_obj_path = config.output_prefix + "_surface" + extension;
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

    if (config.debug) {
        result.outer_smooth_csv_path = config.output_prefix + "_outer_smooth.csv";
        result.inner_smooth_csv_path = config.output_prefix + "_inner_smooth.csv";
        result.smooth_valid_mask_csv_path = config.output_prefix + "_smooth_valid_mask.csv";
        result.metric_domain_mask_csv_path = config.output_prefix + "_smooth_metric_domain_mask.csv";
        result.smooth_non_crossing_adjustment_mask_csv_path =
            config.output_prefix + "_smooth_non_crossing_adjustment_mask.csv";
        result.summary_csv_path = config.output_prefix + "_stage7_summary.csv";
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
            result.outer_mesh_path = config.output_prefix + mesh_prefix + "_outer_surface" + extension;
            result.inner_mesh_path = config.output_prefix + mesh_prefix + "_inner_surface" + extension;
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
            result.outer_mesh_path = config.output_prefix + mesh_prefix + "_surface" + extension;
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
    result.messages.push_back("Geometry Stage 7 max iterations: " + std::to_string(config.stage7_max_iterations));
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
    result.messages.push_back("Geometry Stage 7 outer iterations used: " + std::to_string(result.outer_iterations_used));
    result.messages.push_back("Geometry Stage 7 inner iterations used: " + std::to_string(result.inner_iterations_used));
    result.messages.push_back("Geometry Stage 7 outer final update: " + std::to_string(result.outer_final_max_update));
    result.messages.push_back("Geometry Stage 7 inner final update: " + std::to_string(result.inner_final_max_update));
    result.messages.push_back("Geometry Stage 7 outer fit residual: " + formatScientific(result.outer_fit_final_residual));
    result.messages.push_back("Geometry Stage 7 inner fit residual: " + formatScientific(result.inner_fit_final_residual));
    result.messages.push_back("Geometry Stage 7 smooth valid node count: " + std::to_string(result.smooth_valid_node_count));
    result.messages.push_back("Geometry Stage 7 metric domain node count: " + std::to_string(result.metric_domain_node_count));
    result.messages.push_back("Geometry Stage 7 smooth non-crossing adjusted node count: " +
                              std::to_string(result.smooth_non_crossing_adjusted_node_count));
    result.messages.push_back("Geometry Stage 7 min smooth separation: " + std::to_string(result.min_smooth_separation));
    result.messages.push_back("Geometry Stage 7 max smooth separation: " + std::to_string(result.max_smooth_separation));
    result.messages.push_back("Geometry Stage 7 mean smooth separation: " + std::to_string(result.mean_smooth_separation));
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
    result.success = result.preparation.success && result.stage2_patch.success && result.stage3_patch.success &&
                     result.stage4_raw.success && result.stage5_prep.success && result.stage6_surfaces.success &&
                     result.stage7_smooth.success;
    const std::string run_summary_json_path = config.output_prefix + "_run_summary.json";
    if (!writeGeometryRunSummaryJson(config, parser_config, run_summary_json_path)) {
        throw std::runtime_error("Failed to write geometry run summary JSON: " + run_summary_json_path);
    }
    result.messages.push_back("Geometry run summary JSON: " + run_summary_json_path);
    logMessages(result.messages, logger);

    return result;
}
