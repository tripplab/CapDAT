#include "pdb_parser.hpp"

#include <fstream>
#include <stdexcept>
#include <utility>

/**
 * @brief Construct a PdbParser with explicit configuration and optional logger.
 */
PdbParser::PdbParser(ParserConfig config, Logger* logger)
    : config_(config),
      logger_(logger) {
}

/**
 * @brief Parse a PDB file into a Capsid object.
 *
 * This v01 scaffolding:
 * - opens the file,
 * - scans line by line,
 * - recognizes supported coordinate records,
 * - parses them into temporary PdbRecord objects,
 * - validates and filters them,
 * - reconstructs the Capsid hierarchy in a single pass,
 * - finalizes structure-derived counts before returning.
 *
 * @throws std::runtime_error on fatal input or parsing failure.
 */
Capsid PdbParser::parseFile(const std::string& path) {
    stats_ = ParseStats{};
    current_internal_subunit_id_ = 0;
    has_active_chain_ = false;
    last_chain_id_ = ' ';
    last_residue_seq_ = 0;
    last_insertion_code_ = ' ';

    if (logger_ != nullptr) {
        logger_->info("Opening input file: " + path);
    }

    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Unable to open input file: " + path);
    }

    Capsid capsid(path);

    std::string line;
    bool saw_any_content = false;

    while (std::getline(input, line)) {
        saw_any_content = true;
        ++stats_.total_lines_read;

        if (!isCoordinateRecord(line)) {
            ++stats_.total_skipped_non_coordinate_records;
            continue;
        }

        ++stats_.total_coordinate_records_detected;

        PdbRecord record = parseAtomLikeRecord(line);

        if (!validateRecord(record)) {
            ++stats_.total_malformed_records;
            capsid.incrementSkippedRecordCount();

            if (logger_ != nullptr && config_.verbose_warnings) {
                logger_->warning("Skipped malformed coordinate record at line " +
                                 std::to_string(stats_.total_lines_read));
            }
            continue;
        }

        if (record.alt_loc != ' ') {
            ++stats_.total_altloc_records;
            capsid.incrementAltLocCount();
        }

        if (!shouldAcceptRecord(record)) {
            capsid.incrementSkippedRecordCount();
            continue;
        }

        appendRecordToCapsid(record, capsid);

        ++stats_.total_atoms_accepted;
        capsid.incrementAtomCount();

        if (record.is_hetatm) {
            ++stats_.total_hetatm_accepted;
            capsid.incrementHetatmCount();
        }
    }

    if (!saw_any_content) {
        throw std::runtime_error("Input file is empty: " + path);
    }

    if (stats_.total_atoms_accepted == 0) {
        throw std::runtime_error("No valid coordinate records were accepted from input file: " + path);
    }

    capsid.finalizeCounts();

    if (logger_ != nullptr) {
        logger_->info("Parsing completed successfully");
        logger_->info("Accepted atoms: " + std::to_string(capsid.atomCount()));
        logger_->info("Accepted residues: " + std::to_string(capsid.residueCount()));
        logger_->info("Internal subunits: " + std::to_string(capsid.subunitCount()));
    }

    return capsid;
}

/**
 * @brief Return read-only access to parse statistics from the last run.
 */
const ParseStats& PdbParser::stats() const {
    return stats_;
}

/**
 * @brief Determine whether a line is a supported coordinate record candidate.
 *
 * In v01:
 * - ATOM is always recognized,
 * - HETATM is recognized only if include_hetatm is enabled.
 */
bool PdbParser::isCoordinateRecord(const std::string& line) const {
    if (line.size() < 6) {
        return false;
    }

    const std::string record_name = trim(slice(line, 0, 6));

    if (record_name == "ATOM") {
        return true;
    }

    if (config_.include_hetatm && record_name == "HETATM") {
        return true;
    }

    return false;
}

/**
 * @brief Parse one fixed-column ATOM/HETATM-like record into a temporary PdbRecord.
 *
 * This implementation follows standard PDB fixed-column conventions in a
 * minimal, v01-friendly way.
 */
PdbRecord PdbParser::parseAtomLikeRecord(const std::string& line) const {
    PdbRecord record;
    record.raw_line = line;

    record.record_name = trim(slice(line, 0, 6));
    record.is_hetatm = (record.record_name == "HETATM");

    record.serial = parseInt(slice(line, 6, 5), 0);
    record.atom_name = trim(slice(line, 12, 4));

    {
        const std::string alt = slice(line, 16, 1);
        record.alt_loc = alt.empty() ? ' ' : alt[0];
    }

    record.residue_name = trim(slice(line, 17, 3));

    {
        const std::string chain = slice(line, 21, 1);
        record.chain_id = chain.empty() ? ' ' : chain[0];
    }

    record.residue_seq = parseInt(slice(line, 22, 4), 0);

    {
        const std::string insertion = slice(line, 26, 1);
        record.insertion_code = insertion.empty() ? ' ' : insertion[0];
    }

    record.x = parseDouble(slice(line, 30, 8), 0.0);
    record.y = parseDouble(slice(line, 38, 8), 0.0);
    record.z = parseDouble(slice(line, 46, 8), 0.0);

    record.occupancy = parseDouble(slice(line, 54, 6), 0.0);
    record.temp_factor = parseDouble(slice(line, 60, 6), 0.0);

    record.element = trim(slice(line, 76, 2));
    record.charge = trim(slice(line, 78, 2));

    record.is_valid = true;
    return record;
}

/**
 * @brief Perform basic v01 record-level validation.
 *
 * Minimal validity conditions:
 * - supported record type,
 * - enough line length to contain core coordinate fields,
 * - atom serial number parsed,
 * - residue sequence number parsed,
 * - x/y/z coordinate fields present.
 *
 * This scaffolding uses simple fallback-based checks and does not yet attempt
 * a deeply strict interpretation of all PDB edge cases.
 */
bool PdbParser::validateRecord(const PdbRecord& record) const {
    if (record.record_name != "ATOM" && record.record_name != "HETATM") {
        return false;
    }

    if (record.raw_line.size() < 54) {
        return false;
    }

    // In this initial scaffold, serial/residue_seq default to 0 on parse
    // failure, which is imperfect because 0 may also be a parsed value in some
    // malformed contexts. For v01 foundation scaffolding, this is acceptable,
    // but later refinement should track field-specific parse success directly.
    if (record.serial == 0) {
        return false;
    }

    if (record.residue_seq == 0 && trim(slice(record.raw_line, 22, 4)) != "0") {
        return false;
    }

    if (trim(slice(record.raw_line, 30, 8)).empty() ||
        trim(slice(record.raw_line, 38, 8)).empty() ||
        trim(slice(record.raw_line, 46, 8)).empty()) {
        return false;
    }

    return record.is_valid;
}

/**
 * @brief Decide whether a parsed record should be retained internally.
 *
 * Current v01 scaffold behavior:
 * - ATOM records are accepted,
 * - HETATM records are accepted only if include_hetatm is enabled,
 * - if protein_only is enabled, a simple residue-name heuristic is applied.
 *
 * The protein-only filter is deliberately conservative and should be refined in
 * later iterations if needed.
 */
bool PdbParser::shouldAcceptRecord(const PdbRecord& record) const {
    if (record.is_hetatm && !config_.include_hetatm) {
        return false;
    }

    if (!config_.protein_only) {
        return true;
    }

    static const char* kProteinResidues[] = {
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
        "SEC", "PYL", "ASX", "GLX", "UNK"
    };

    for (const char* residue_name : kProteinResidues) {
        if (record.residue_name == residue_name) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Append one accepted record into the Capsid hierarchy.
 *
 * This method performs single-pass reconstruction under the v01 assumption that
 * the input is in conventional PDB order and that subunits can be tracked
 * contiguously.
 *
 * Current scaffold policy for starting a new internal subunit:
 * - first accepted record starts the first subunit,
 * - a change in raw PDB chain label starts a new internal subunit.
 *
 * This is intentionally simple and will likely need refinement later to handle
 * more complex chain-label reuse patterns in large capsids.
 */
void PdbParser::appendRecordToCapsid(const PdbRecord& record, Capsid& capsid) {
    const bool start_new_chain =
        !has_active_chain_ || (record.chain_id != last_chain_id_);

    if (start_new_chain) {
        ++current_internal_subunit_id_;
        has_active_chain_ = true;
        last_chain_id_ = record.chain_id;
        last_residue_seq_ = record.residue_seq;
        last_insertion_code_ = record.insertion_code;

        Chain chain(current_internal_subunit_id_, record.chain_id);
        Residue residue(record.residue_name,
                        record.residue_seq,
                        record.insertion_code,
                        record.chain_id,
                        current_internal_subunit_id_);

        residue.addAtom(Atom(record.serial,
                             record.atom_name,
                             record.alt_loc,
                             record.residue_name,
                             record.chain_id,
                             record.residue_seq,
                             record.insertion_code,
                             record.x,
                             record.y,
                             record.z,
                             record.occupancy,
                             record.temp_factor,
                             record.element,
                             record.charge,
                             record.is_hetatm));

        chain.addResidue(std::move(residue));
        capsid.addChain(std::move(chain));

        ++stats_.total_internal_subunits_created;
        ++stats_.total_residues_created;
        capsid.incrementSubunitCount();
        capsid.incrementResidueCount();

        return;
    }

    const bool start_new_residue =
        (record.residue_seq != last_residue_seq_) ||
        (record.insertion_code != last_insertion_code_);

    if (start_new_residue) {
        Residue residue(record.residue_name,
                        record.residue_seq,
                        record.insertion_code,
                        record.chain_id,
                        current_internal_subunit_id_);

        residue.addAtom(Atom(record.serial,
                             record.atom_name,
                             record.alt_loc,
                             record.residue_name,
                             record.chain_id,
                             record.residue_seq,
                             record.insertion_code,
                             record.x,
                             record.y,
                             record.z,
                             record.occupancy,
                             record.temp_factor,
                             record.element,
                             record.charge,
                             record.is_hetatm));

        capsid.lastChain().addResidue(std::move(residue));

        last_residue_seq_ = record.residue_seq;
        last_insertion_code_ = record.insertion_code;

        ++stats_.total_residues_created;
        capsid.incrementResidueCount();

        return;
    }

    capsid.addAtomToLastChainResidue(Atom(record.serial,
                                          record.atom_name,
                                          record.alt_loc,
                                          record.residue_name,
                                          record.chain_id,
                                          record.residue_seq,
                                          record.insertion_code,
                                          record.x,
                                          record.y,
                                          record.z,
                                          record.occupancy,
                                          record.temp_factor,
                                          record.element,
                                          record.charge,
                                          record.is_hetatm));
}

/**
 * @brief Trim leading and trailing whitespace from a string.
 */
std::string PdbParser::trim(const std::string& s) const {
    const std::size_t first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }

    const std::size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

/**
 * @brief Extract a substring using fixed-column indexing.
 *
 * If the requested start is beyond the available line length, an empty string
 * is returned. If the requested width extends past the line end, the remaining
 * tail is returned.
 */
std::string PdbParser::slice(const std::string& line,
                             std::size_t start,
                             std::size_t length) const {
    if (start >= line.size()) {
        return "";
    }

    return line.substr(start, length);
}

/**
 * @brief Parse an integer field from text, returning fallback on failure.
 */
int PdbParser::parseInt(const std::string& text, int fallback) const {
    const std::string cleaned = trim(text);
    if (cleaned.empty()) {
        return fallback;
    }

    try {
        return std::stoi(cleaned);
    } catch (...) {
        return fallback;
    }
}

/**
 * @brief Parse a floating-point field from text, returning fallback on failure.
 */
double PdbParser::parseDouble(const std::string& text, double fallback) const {
    const std::string cleaned = trim(text);
    if (cleaned.empty()) {
        return fallback;
    }

    try {
        return std::stod(cleaned);
    } catch (...) {
        return fallback;
    }
}

// NOTE ON FOUNDATION-LEVEL VALIDATION:
//
// This first parser scaffold uses simple fallback-based numeric parsing and
// lightweight validity checks. That is acceptable for a v01 skeleton because
// the immediate goal is a clean, testable parsing flow that compiles and
// expresses the intended architecture.
//
// A later refinement should track field-specific parse success more explicitly
// so that malformed numeric fields can be distinguished more rigorously from
// legitimate values.

// NOTE ON PROTEIN-ONLY FILTERING:
//
// The current protein-only policy uses a residue-name heuristic. This is a
// practical starting point for v01 because the technical design explicitly says
// that non-protein molecules should be ignored internally, ideally under CLI
// control.
//
// This heuristic is intentionally conservative and may need refinement later,
// especially for modified residues, unusual naming conventions, or edge cases
// where protein-derived atoms appear in HETATM records.

// NOTE ON INTERNAL SUBUNIT RECONSTRUCTION:
//
// The current scaffold starts a new internal subunit when the raw PDB chain
// label changes. This is simple and compilable, but it does not yet fully solve
// the larger capsid problem where chain labels may be reused across multiple
// distinct contiguous subunits.
//
// The architecture is ready for that refinement, but the precise grouping rule
// should be strengthened in a later pass once we finalize the intended parser
// behavior for repeated chain-label blocks.

// NOTE ON PARSER MUTATION API:
//
// This revision removes the earlier const_cast-based access pattern. The parser
// now mutates the active structure only through explicit Capsid and Chain
// methods, which is much cleaner and easier to audit.
//
// For v01, this is the right compromise: it preserves the single-pass parser
// design while keeping const-correctness and ownership boundaries clearer.
