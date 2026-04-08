#ifndef CAPDAT_ACCEPTED_STRUCTURE_WRITER_HPP
#define CAPDAT_ACCEPTED_STRUCTURE_WRITER_HPP

#include <cstddef>
#include <string>

#include "capsid.hpp"
#include "logger.hpp"
#include "pdb_parser.hpp"

/**
 * @brief Configuration for exporting accepted atoms as a clean PDB-like file.
 */
struct AcceptedStructureWriterConfig {
    std::string output_path;
    bool emit_header_comments = true;
    bool emit_ter_records = true;
    bool emit_end_record = true;
    bool preserve_atom_serial_numbers = true;
    bool coordinate_records_only = false;
};

/**
 * @brief Summary statistics produced by accepted-structure export.
 */
struct AcceptedStructureWriteStats {
    std::size_t atoms_written = 0;
    std::size_t residues_written = 0;
    std::size_t subunits_written = 0;
    std::size_t hetatm_written = 0;
    std::size_t altloc_written = 0;
    std::size_t blank_chain_id_count = 0;
    std::size_t repeated_chain_id_groups = 0;
};

/**
 * @brief Writes accepted atoms from Capsid hierarchy into a clean PDB-like file.
 */
class AcceptedStructureWriter {
public:
    explicit AcceptedStructureWriter(Logger* logger);

    [[nodiscard]] AcceptedStructureWriteStats write(const Capsid& capsid,
                                                    const AcceptedStructureWriterConfig& config,
                                                    const ParserConfig& parser_config) const;

private:
    [[nodiscard]] std::string formatAtomRecord(const Atom& atom, int serial_number) const;
    [[nodiscard]] std::string trimToWidth(const std::string& value, std::size_t width) const;

private:
    Logger* logger_ = nullptr;
};

#endif // CAPDAT_ACCEPTED_STRUCTURE_WRITER_HPP
