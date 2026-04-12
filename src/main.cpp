#include "export_capsid.hpp"
#include "geometry_analysis.hpp"
#include "logger.hpp"
#include "pdb_parser.hpp"
#include "reorientation_workflow.hpp"
#include "structural_summary.hpp"
#include "summary_reporter.hpp"
#include "timer.hpp"

#include <exception>
#include <iostream>
#include <cstddef>
#include <string>

/**
 * @brief Print the CapDAT help message to standard output.
 */
void printHelp(const std::string& program_name) {
    std::cout
        << "CapDAT v0.1.0\n"
        << "Capsid Data Analysis Toolkit - v01 foundation release\n\n"
        << "Usage:\n"
        << "  " << program_name << " --input <file> [options]\n\n"
        << "Required options:\n"
        << "  -i, --input <file>      Input PDB file\n\n"
        << "Optional options:\n"
        << "  -l, --log <file>        Write log output to file\n"
        << "  -v, --verbose           Increase terminal verbosity\n"
        << "      --include-hetatm    Include HETATM records\n"
        << "      --export-final <f>  Export current accepted Capsid coordinates to file\n"
        << "      --reorient          Enable post-parse in-place reorientation workflow\n"
        << "      --align-fold <name> Reorientation source: canonical fold (2_0,2_1,3_0,3_1,5_0)\n"
        << "      --align-vector <v>  Reorientation source: custom direction x,y,z from origin\n"
        << "      --align-axis <a>    Target alignment axis x|y|z (default: z)\n"
        << "      --geometry-analysis Run geometry analysis Stage 1 preparation\n"
        << "      --geometry_fold_type <n>   Geometry fold type 2|3|5 (default: 2)\n"
        << "      --geometry_fold_index <n>  Geometry fold index for selected type (default: 0)\n"
        << "      --geometry_cylinder_radius <A>  Geometry cylinder radius in angstroms (default: 12.0)\n"
        << "      --geometry_grid_spacing <A>  Geometry Stage 4 XY grid spacing in angstroms (default: 2.0)\n"
        << "      --geometry_min_atoms_in_patch <n>  Minimum selected atoms required (default: 20)\n"
        << "      --geometry_out_prefix <path>  Prefix for geometry analysis outputs (default: geometry)\n"
        << "      --quiet             Reduce terminal output\n"
        << "  -h, --help              Show this help message\n"
        << "      --version           Show version information\n\n"
        << "Notes:\n"
        << "  - --write-clean-pdb has been replaced by --export-final.\n"
        << "  - Reorientation applies a pure rotation (no translation) in place to\n"
        << "    current Capsid coordinates, changing frame identity from original\n"
        << "    parsed frame to a derived frame tracked in Capsid orientation state.\n"
        << "  - --align-axis defaults to positive Z when omitted.\n\n"
        << "Examples:\n"
        << "  " << program_name << " --input capsid.pdb --export-final accepted.pdb\n"
        << "  " << program_name << " -i capsid.pdb --reorient --align-fold 5_0 --align-axis x --export-final aligned.pdb\n";
}

void printVersion() {
    std::cout << "CapDAT v" << CAPDAT_VERSION << '\n';
}

void printSummary(const Capsid& capsid,
                  const ParseStats& stats,
                  const StructuralSummary& structural_summary) {
    std::cout << "Input file:              " << capsid.sourcePath() << '\n';
    std::cout << "Total lines read:        " << stats.total_lines_read << '\n';
    std::cout << "Coordinate records seen: " << stats.total_coordinate_records_detected << '\n';
    std::cout << "Accepted atoms:          " << capsid.atomCount() << '\n';
    std::cout << "Accepted residues:       " << capsid.residueCount() << '\n';
    std::cout << "Internal subunits:       " << capsid.subunitCount() << '\n';
    std::cout << "Accepted HETATM:         " << capsid.hetatmCount() << '\n';
    std::cout << "Alternate locations:     " << capsid.altLocCount() << '\n';
    std::cout << "Malformed records:       " << stats.total_malformed_records << '\n';
    std::cout << "Skipped records:         " << capsid.skippedRecordCount() << '\n';
    printStructuralSummaryBlock(std::cout, structural_summary);
}

int main(int argc, char* argv[]) {
    std::string input_path;
    std::string log_path;
    std::string export_final_output_path;
    bool verbose = false;
    bool quiet = false;
    bool include_hetatm = false;

    // CLI flags for strict opt-in reorientation workflow.
    bool reorient_requested = false;
    bool align_fold_given = false;
    bool align_vector_given = false;
    std::string align_fold_name;
    std::string align_vector_text;
    char align_axis = 'z';

    bool geometry_analysis_requested = false;
    int geometry_fold_type = 2;
    int geometry_fold_index = 0;
    double geometry_cylinder_radius = 12.0;
    double geometry_grid_spacing = 2.0;
    std::size_t geometry_min_atoms_in_patch = 20;
    std::string geometry_output_prefix = "geometry";

    const std::string program_name = (argc > 0) ? argv[0] : "capsid_analyzer";

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printHelp(program_name);
            return 0;
        }
        if (arg == "--version") {
            printVersion();
            return 0;
        }
        if (arg == "-i" || arg == "--input") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for " << arg << '\n';
                return 1;
            }
            input_path = argv[++i];
            continue;
        }
        if (arg == "-l" || arg == "--log") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for " << arg << '\n';
                return 1;
            }
            log_path = argv[++i];
            continue;
        }
        if (arg == "-v" || arg == "--verbose") {
            verbose = true;
            continue;
        }
        if (arg == "--quiet") {
            quiet = true;
            continue;
        }
        if (arg == "--include-hetatm") {
            include_hetatm = true;
            continue;
        }
        if (arg == "--export-final") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for " << arg << '\n';
                return 1;
            }
            export_final_output_path = argv[++i];
            continue;
        }
        if (arg == "--write-clean-pdb") {
            std::cerr << "Error: --write-clean-pdb is deprecated. Use --export-final <file>.\n";
            return 1;
        }
        if (arg == "--reorient") {
            reorient_requested = true;
            continue;
        }
        if (arg == "--align-fold") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --align-fold\n";
                return 1;
            }
            align_fold_given = true;
            align_fold_name = argv[++i];
            continue;
        }
        if (arg == "--align-vector") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --align-vector\n";
                return 1;
            }
            align_vector_given = true;
            align_vector_text = argv[++i];
            continue;
        }
        if (arg == "--align-axis") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --align-axis\n";
                return 1;
            }
            const std::string axis_arg = argv[++i];
            if (axis_arg.size() != 1 || (axis_arg[0] != 'x' && axis_arg[0] != 'y' && axis_arg[0] != 'z')) {
                std::cerr << "Error: --align-axis must be one of x, y, z\n";
                return 1;
            }
            align_axis = axis_arg[0];
            continue;
        }
        if (arg == "--geometry-analysis") {
            geometry_analysis_requested = true;
            continue;
        }
        if (arg == "--geometry_fold_type") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --geometry_fold_type\n";
                return 1;
            }
            geometry_fold_type = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--geometry_fold_index") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --geometry_fold_index\n";
                return 1;
            }
            geometry_fold_index = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--geometry_cylinder_radius") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --geometry_cylinder_radius\n";
                return 1;
            }
            geometry_cylinder_radius = std::stod(argv[++i]);
            continue;
        }
        if (arg == "--geometry_grid_spacing") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --geometry_grid_spacing\n";
                return 1;
            }
            geometry_grid_spacing = std::stod(argv[++i]);
            continue;
        }
        if (arg == "--geometry_min_atoms_in_patch") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --geometry_min_atoms_in_patch\n";
                return 1;
            }
            geometry_min_atoms_in_patch = static_cast<std::size_t>(std::stoul(argv[++i]));
            continue;
        }
        if (arg == "--geometry_out_prefix") {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for --geometry_out_prefix\n";
                return 1;
            }
            geometry_output_prefix = argv[++i];
            continue;
        }

        std::cerr << "Error: unknown argument: " << arg << '\n';
        return 1;
    }

    if (input_path.empty()) {
        std::cerr << "Error: input file is required. Use --input <file>.\n";
        return 1;
    }

    if (align_fold_given && align_vector_given) {
        std::cerr << "Error: --align-fold and --align-vector are mutually exclusive.\n";
        return 1;
    }
    if (reorient_requested && !align_fold_given && !align_vector_given) {
        std::cerr << "Error: --reorient requires exactly one source: --align-fold or --align-vector.\n";
        return 1;
    }
    if ((align_fold_given || align_vector_given) && !reorient_requested) {
        std::cerr << "Error: --align-fold/--align-vector require --reorient.\n";
        return 1;
    }
    if (geometry_analysis_requested && reorient_requested) {
        std::cerr << "Error: --geometry-analysis cannot be combined with --reorient.\n";
        return 1;
    }

    Logger logger;
    if (quiet) {
        logger.setVerbosity(LogLevel::WARNING);
    } else if (verbose) {
        logger.setVerbosity(LogLevel::DEBUG);
    } else {
        logger.setVerbosity(LogLevel::INFO);
    }

    try {
        if (!log_path.empty()) {
            logger.setLogFile(log_path);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 2;
    }

    ParserConfig config;
    config.include_hetatm = include_hetatm;
    config.strict_mode = false;
    config.keep_altloc_all = true;
    config.ignore_blank_chain_id = false;
    config.verbose_warnings = verbose;
    config.protein_only = true;

    Timer timer;
    timer.start();

    try {
        logger.info("Starting CapDAT");

        PdbParser parser(config, &logger);
        Capsid capsid = parser.parseFile(input_path);

        logger.info("Starting extended structural summary geometry");
        StructuralSummary structural_summary = computeStructuralSummary(capsid);
        printSummary(capsid, parser.stats(), structural_summary);
        logger.info("Completed extended structural summary geometry");

        ReorientationRequest reorient_request;
        reorient_request.request_reorientation = reorient_requested;
        reorient_request.source_mode = align_vector_given ? ReorientationSourceMode::custom_vector
                                                          : ReorientationSourceMode::fold;
        reorient_request.fold_name = align_fold_name;
        reorient_request.custom_vector_text = align_vector_text;
        reorient_request.target_axis = align_axis;
        reorient_request.request_export = !export_final_output_path.empty();
        reorient_request.export_path = export_final_output_path;
        reorient_request.verbose = verbose;

        const ReorientationResult reorient_result =
            applyReorientationWorkflow(capsid, reorient_request, &logger);
        (void)reorient_result;

        FoldPatchAnalysisConfig geometry_config;
        geometry_config.enabled = geometry_analysis_requested;
        geometry_config.fold_type = geometry_fold_type;
        geometry_config.fold_index = geometry_fold_index;
        geometry_config.cylinder_radius = geometry_cylinder_radius;
        geometry_config.grid_spacing = geometry_grid_spacing;
        geometry_config.min_atoms_in_patch = geometry_min_atoms_in_patch;
        geometry_config.export_rotated_capsid = verbose;
        geometry_config.output_prefix = geometry_output_prefix;

        const GeometryAnalysisResult geometry_result =
            runFoldPatchGeometryAnalysis(capsid, geometry_config, config, &logger);
        if (!geometry_result.success) {
            throw std::runtime_error("Geometry analysis failed in Stage 1/2/3/4 pipeline");
        }

        if (!export_final_output_path.empty()) {
            ExportCapsidConfig writer_config;
            writer_config.output_path = export_final_output_path;
            writer_config.emit_header_comments = true;
            writer_config.emit_ter_records = true;
            writer_config.emit_end_record = true;
            writer_config.preserve_atom_serial_numbers = true;
            writer_config.coordinate_records_only = false;

            ExportCapsidWriter writer(&logger);
            const ExportCapsidStats write_stats = writer.write(capsid, writer_config, config);
            logger.info("Final exported atoms written: " + std::to_string(write_stats.atoms_written));
            logger.info("Final structure exported successfully: " + export_final_output_path);
        }

        timer.stop();
        std::cout << "Runtime (s):             " << timer.elapsedSeconds() << '\n';
        logger.info("Run completed successfully");
        return 0;
    } catch (const std::runtime_error& e) {
        timer.stop();
        logger.error(e.what());
        return 4;
    } catch (const std::exception& e) {
        timer.stop();
        logger.error(std::string("Unexpected internal failure: ") + e.what());
        return 4;
    } catch (...) {
        timer.stop();
        logger.error("Unexpected internal failure: unknown exception");
        return 4;
    }
}
