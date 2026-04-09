#ifndef CAPDAT_EXPORT_CAPSID_HPP
#define CAPDAT_EXPORT_CAPSID_HPP

#include <cstddef>
#include <string>

#include "capsid.hpp"
#include "logger.hpp"
#include "pdb_parser.hpp"

/**
 * @brief Configuration for exporting the current accepted Capsid structure state.
 *
 * Export semantics convention:
 * - This writer emits whatever coordinates are currently stored in Capsid.
 * - It does not promise byte-for-byte reproduction of the source file.
 * - If a reorientation workflow rotated coordinates in place, those transformed
 *   coordinates are exactly what will be written.
 */
struct ExportCapsidConfig {
    std::string output_path;
    bool emit_header_comments = true;
    bool emit_ter_records = true;
    bool emit_end_record = true;
    bool preserve_atom_serial_numbers = true;
    bool coordinate_records_only = false;
};

/**
 * @brief Summary statistics produced by final accepted-structure export.
 */
struct ExportCapsidStats {
    std::size_t atoms_written = 0;
    std::size_t residues_written = 0;
    std::size_t subunits_written = 0;
    std::size_t hetatm_written = 0;
    std::size_t altloc_written = 0;
    std::size_t blank_chain_id_count = 0;
    std::size_t repeated_chain_id_groups = 0;
};

/**
 * @brief Writes the current in-memory accepted Capsid structure to PDB-like text.
 */
class ExportCapsidWriter {
public:
    explicit ExportCapsidWriter(Logger* logger);

    [[nodiscard]] ExportCapsidStats write(const Capsid& capsid,
                                          const ExportCapsidConfig& config,
                                          const ParserConfig& parser_config) const;

private:
    [[nodiscard]] std::string formatAtomRecord(const Atom& atom, int serial_number) const;
    [[nodiscard]] std::string trimToWidth(const std::string& value, std::size_t width) const;

private:
    Logger* logger_ = nullptr;
};

#endif // CAPDAT_EXPORT_CAPSID_HPP
