#ifndef CAPDAT_PDB_PARSER_HPP
#define CAPDAT_PDB_PARSER_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "capsid.hpp"
#include "logger.hpp"

/**
 * @brief Configuration options controlling PDB parsing behavior.
 *
 * Keeping parser behavior explicit in a dedicated configuration structure helps
 * avoid hidden policy decisions inside the parser implementation.
 *
 * Not all fields need to be exposed on the CLI in v01, but they are useful to
 * define now so the parsing layer remains easy to extend in later versions.
 */
struct ParserConfig {
    /// Include HETATM records in the accepted structure.
    bool include_hetatm = false;

    /// Escalate selected warnings into stricter validation behavior.
    bool strict_mode = false;

    /// Keep all alternate-location records as independent atoms.
    bool keep_altloc_all = true;

    /// Treat blank chain identifiers as problematic if enabled.
    bool ignore_blank_chain_id = false;

    /// Emit more detailed parsing warnings.
    bool verbose_warnings = false;

    /**
     * @brief Keep only atoms belonging to capsid proteins.
     *
     * In v01, the implementation objective explicitly states that non-protein
     * molecules such as DNA, RNA, water, ions, ligands, and other non-capsid
     * components should be ignored internally, ideally under CLI control.
     */
    bool protein_only = true;
};

/**
 * @brief Aggregate statistics collected during parsing.
 *
 * This structure is intended for:
 * - final reporting,
 * - debugging,
 * - validation checks,
 * - future tests.
 *
 * The fields are intentionally simple and mostly count-based in v01.
 */
struct ParseStats {
    std::size_t total_lines_read = 0;
    std::size_t total_coordinate_records_detected = 0;
    std::size_t total_atoms_accepted = 0;
    std::size_t total_hetatm_accepted = 0;
    std::size_t total_malformed_records = 0;
    std::size_t total_skipped_non_coordinate_records = 0;
    std::size_t total_warnings = 0;
    std::size_t total_altloc_records = 0;
    std::size_t total_internal_subunits_created = 0;
    std::size_t total_residues_created = 0;
};

/**
 * @brief Temporary representation of one parsed ATOM/HETATM-like line.
 *
 * PdbRecord is not part of the persistent domain model. Its purpose is to keep
 * fixed-column extraction and early validation separate from the construction
 * of Atom/Residue/Chain/Capsid objects.
 *
 * This helps keep parser logic readable and easier to test.
 */
struct PdbRecord {
    std::string record_name;
    std::string raw_line;

    int serial = 0;
    std::string atom_name;
    char alt_loc = ' ';
    std::string residue_name;
    char chain_id = ' ';
    int residue_seq = 0;
    char insertion_code = ' ';

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double occupancy = 0.0;
    double temp_factor = 0.0;

    std::string element;
    std::string charge;

    bool is_hetatm = false;
    bool is_valid = false;
};

/**
 * @brief Parse a PDB file and reconstruct the v01 Capsid hierarchy.
 *
 * PdbParser is the only component in v01 that should understand the fixed-
 * column layout of PDB coordinate records. Its responsibilities include:
 *
 * - opening and reading the input file,
 * - recognizing supported coordinate records,
 * - extracting fields from fixed columns,
 * - applying parser policy and validation,
 * - reconstructing internal subunits and residues,
 * - populating a Capsid object,
 * - collecting parsing statistics and warnings.
 */
class PdbParser {
public:
    /**
     * @brief Construct a parser with explicit configuration and logging support.
     *
     * @param config Parser behavior configuration.
     * @param logger Logger used for user-facing and diagnostic messages.
     */
    PdbParser(ParserConfig config, Logger* logger);

    /**
     * @brief Parse a PDB file into a Capsid object.
     *
     * @param path Input file path.
     * @return Fully constructed Capsid object.
     *
     * Fatal errors should be reported clearly and signaled via exceptions or
     * another explicit failure mechanism in the implementation.
     */
    [[nodiscard]] Capsid parseFile(const std::string& path);

    /// @return Read-only access to parse statistics from the last run.
    [[nodiscard]] const ParseStats& stats() const;

private:
    /**
     * @brief Determine whether a line is a supported coordinate record.
     *
     * In v01, this primarily means ATOM and optionally HETATM depending on the
     * parser configuration.
     *
     * @param line Raw input line.
     * @return True if the line is a supported coordinate record candidate.
     */
    [[nodiscard]] bool isCoordinateRecord(const std::string& line) const;

    /**
     * @brief Parse one fixed-column ATOM/HETATM-like record.
     *
     * @param line Raw input line.
     * @return Temporary parsed record with validity information.
     */
    [[nodiscard]] PdbRecord parseAtomLikeRecord(const std::string& line) const;

    /**
     * @brief Perform basic record-level validation.
     *
     * This method should check whether required fields were extracted well
     * enough for v01 ingestion.
     *
     * @param record Parsed temporary record.
     * @return True if the record is valid enough to ingest.
     */
    [[nodiscard]] bool validateRecord(const PdbRecord& record) const;

    /**
     * @brief Decide whether a parsed record should be kept internally.
     *
     * This is the natural place to implement policies such as:
     * - include/exclude HETATM,
     * - protein-only filtering,
     * - future stricter selection rules.
     *
     * @param record Parsed temporary record.
     * @return True if the record should be accepted into the Capsid structure.
     */
    [[nodiscard]] bool shouldAcceptRecord(const PdbRecord& record) const;

    /**
     * @brief Append one accepted record into the Capsid hierarchy.
     *
     * This method is responsible for driving the single-pass reconstruction of:
     * - internal subunits,
     * - residues,
     * - atoms.
     *
     * @param record Parsed and accepted record.
     * @param capsid Capsid object being populated.
     */
    void appendRecordToCapsid(const PdbRecord& record, Capsid& capsid);

    /**
     * @brief Trim leading and trailing whitespace from a string.
     *
     * @param s Input string.
     * @return Trimmed string.
     */
    [[nodiscard]] std::string trim(const std::string& s) const;

    /**
     * @brief Extract a substring using PDB-style fixed-column indexing.
     *
     * This helper keeps fixed-column slicing explicit and centralized.
     *
     * @param line Raw input line.
     * @param start Zero-based start index.
     * @param length Number of characters to extract.
     * @return Extracted substring, possibly shorter if the line is short.
     */
    [[nodiscard]] std::string slice(const std::string& line,
                                    std::size_t start,
                                    std::size_t length) const;

    /**
     * @brief Parse an integer field from text.
     *
     * @param text Raw field text.
     * @param fallback Value returned if parsing fails.
     * @return Parsed integer or fallback.
     */
    [[nodiscard]] int parseInt(const std::string& text, int fallback = 0) const;

    /**
     * @brief Parse a floating-point field from text.
     *
     * @param text Raw field text.
     * @param fallback Value returned if parsing fails.
     * @return Parsed floating-point value or fallback.
     */
    [[nodiscard]] double parseDouble(const std::string& text,
                                     double fallback = 0.0) const;

private:
    // -------------------------------------------------------------------------
    // Parser state
    // -------------------------------------------------------------------------

    ParserConfig config_{};
    ParseStats stats_{};
    Logger* logger_ = nullptr;

    // -------------------------------------------------------------------------
    // Single-pass reconstruction state
    // -------------------------------------------------------------------------
    //
    // These fields allow the parser to detect when it must start a new internal
    // subunit or residue while reading the file line by line in conventional
    // PDB order.

    std::size_t current_internal_subunit_id_ = 0;
    bool has_active_chain_ = false;

    char last_chain_id_ = ' ';
    int last_residue_seq_ = 0;
    char last_insertion_code_ = ' ';
};

// NOTE ON PdbRecord AS A TEMPORARY TYPE:
//
// PdbRecord exists to separate low-level fixed-column parsing from domain-model
// construction. This is intentional: parsing code tends to become messy, and
// mixing raw text extraction directly with hierarchy-building logic usually
// makes the code harder to test, debug, and evolve.
//
// If later development introduces support for additional formats such as mmCIF,
// this separation will become even more valuable.

// NOTE ON LOGGER OWNERSHIP:
//
// The parser stores a raw Logger* instead of owning a logger instance. This is
// deliberate in v01: the parser should use logging infrastructure provided by
// the application layer, not control its lifetime.
//
// A future revision could switch to a reference, smart pointer, or a more
// abstract logging interface if that becomes architecturally useful. For v01,
// a non-owning pointer keeps the coupling low and the wiring simple.

// NOTE ON SINGLE-PASS PARSING STATE:
//
// The parser keeps a small amount of mutable state to support single-pass
// hierarchy reconstruction. This is a practical design choice for v01 because
// it avoids reparsing and keeps memory overhead low.
//
// If later requirements include out-of-order reconstruction, more complex
// validation, or parallel parsing experiments, this state model may need to be
// revised. For the foundation release, it is a reasonable and efficient choice.

// A small forward-looking comment: the protein_only flag is worth defining now even if the first parser implementation uses a simple heuristic, because that filtering requirement is already explicit in the technical spec.

#endif // CAPDAT_PDB_PARSER_HPP








