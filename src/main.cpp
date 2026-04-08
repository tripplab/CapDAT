#include "logger.hpp"
#include "pdb_parser.hpp"
#include "timer.hpp"

#include <exception>
#include <iostream>
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
        << "      --quiet             Reduce terminal output\n"
        << "  -h, --help              Show this help message\n"
        << "      --version           Show version information\n\n"
        << "Examples:\n"
        << "  " << program_name << " --input capsid.pdb\n"
        << "  " << program_name << " -i capsid.pdb --log run.log --verbose\n";
}

/**
 * @brief Print the CapDAT version string to standard output.
 */
void printVersion() {
    std::cout << "CapDAT v" << CAPDAT_VERSION << '\n';
}

/**
 * @brief Print a compact final summary for a successful run.
 */
void printSummary(const Capsid& capsid, const ParseStats& stats, double elapsed_seconds) {
    std::cout << "\n=== CapDAT Summary ===\n";
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
    std::cout << "Runtime (s):             " << elapsed_seconds << '\n';
}

/**
 * @brief Program entry point for CapDAT v01.
 *
 * Exit-code policy:
 * - 0: success
 * - 1: user or argument error
 * - 2: file access error
 * - 3: parsing failure
 * - 4: unexpected internal failure
 */
int main(int argc, char* argv[]) {
    std::string input_path;
    std::string log_path;
    bool verbose = false;
    bool quiet = false;
    bool include_hetatm = false;

    const std::string program_name = (argc > 0) ? argv[0] : "capsid_analyzer";

    // -------------------------------------------------------------------------
    // CLI parsing
    // -------------------------------------------------------------------------

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

        std::cerr << "Error: unknown argument: " << arg << '\n';
        return 1;
    }

    if (input_path.empty()) {
        std::cerr << "Error: input file is required. Use --input <file>.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Logger setup
    // -------------------------------------------------------------------------

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

    // -------------------------------------------------------------------------
    // Parser configuration
    // -------------------------------------------------------------------------

    ParserConfig config;
    config.include_hetatm = include_hetatm;
    config.strict_mode = false;
    config.keep_altloc_all = true;
    config.ignore_blank_chain_id = false;
    config.verbose_warnings = verbose;
    config.protein_only = true;

    // -------------------------------------------------------------------------
    // Execution
    // -------------------------------------------------------------------------

    Timer timer;
    timer.start();

    try {
        logger.info("Starting CapDAT");
        logger.info("Input file: " + input_path);

        if (!log_path.empty()) {
            logger.info("Log file: " + log_path);
        }

        logger.info(std::string("Include HETATM: ") + (include_hetatm ? "true" : "false"));
        logger.info(std::string("Protein only: true"));

        PdbParser parser(config, &logger);
        Capsid capsid = parser.parseFile(input_path);

        timer.stop();

        printSummary(capsid, parser.stats(), timer.elapsedSeconds());

        logger.info("Run completed successfully");
        return 0;
    } catch (const std::runtime_error& e) {
        timer.stop();
        logger.error(e.what());

        const std::string message = e.what();

        // Simple v01 categorization of failure classes.
        if (message.find("open input file") != std::string::npos ||
            message.find("Input file is empty") != std::string::npos) {
            return 2;
        }

        if (message.find("No valid coordinate records") != std::string::npos ||
            message.find("parse") != std::string::npos ||
            message.find("Parsing") != std::string::npos) {
            return 3;
        }

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

// NOTE ON MANUAL CLI PARSING:
//
// The technical design for v01 explicitly avoids external dependencies, so the
// CLI is parsed manually here instead of using a library. This keeps the build
// simple and portable for the foundation release.
//
// If later versions need richer argument handling, validation, or subcommands,
// introducing a dedicated CLI library could become worthwhile.

// NOTE ON EXIT-CODE MAPPING:
//
// The current runtime_error-to-exit-code mapping is intentionally simple and
// string-based. This is acceptable for an early implementation skeleton, but it
// is not the ideal long-term mechanism.
//
// A future refinement should likely introduce more explicit exception types or
// a small error-category system so exit-code behavior becomes more structured
// and less dependent on message text.

// NOTE ON SUMMARY OUTPUT:
//
// The summary combines hierarchy-derived counts from Capsid with parser-event
// counts from ParseStats. This split matches the design intent of v01: some
// values describe the final stored structure, while others describe events that
// occurred during parsing.
//
// That separation is useful and should probably be preserved in later versions
// even if reporting becomes more sophisticated.
