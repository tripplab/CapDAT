#include "accepted_structure_writer.hpp"

#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>

AcceptedStructureWriter::AcceptedStructureWriter(Logger* logger)
    : logger_(logger) {
}

AcceptedStructureWriteStats AcceptedStructureWriter::write(const Capsid& capsid,
                                                           const AcceptedStructureWriterConfig& config,
                                                           const ParserConfig& parser_config) const {
    if (config.output_path.empty()) {
        throw std::runtime_error("Accepted structure export requires a non-empty output path.");
    }

    if (capsid.chains().empty()) {
        throw std::runtime_error("Accepted structure export failed: Capsid hierarchy is empty.");
    }

    if (logger_ != nullptr) {
        logger_->info("Accepted-coordinate export requested");
        logger_->info("Export path: " + config.output_path);
    }

    std::ofstream output(config.output_path, std::ios::out | std::ios::trunc);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to open accepted-coordinate output file: " + config.output_path);
    }

    if (logger_ != nullptr) {
        logger_->info("Writing accepted-coordinate export");
    }

    AcceptedStructureWriteStats stats;
    stats.subunits_written = capsid.chains().size();
    stats.residues_written = capsid.residueCount();

    std::map<char, std::size_t> chain_id_frequency;

    if (config.emit_header_comments && !config.coordinate_records_only) {
        const std::time_t now = std::time(nullptr);
        std::tm tm_snapshot{};
#if defined(_WIN32)
        localtime_s(&tm_snapshot, &now);
#else
        localtime_r(&now, &tm_snapshot);
#endif

        output << "REMARK   1 CapDAT accepted-coordinate export\n";
        output << "REMARK   1 Source: " << capsid.sourcePath() << "\n";
        output << "REMARK   1 Version: " << CAPDAT_VERSION << "\n";
        output << "REMARK   1 Protein-only filter: " << (parser_config.protein_only ? "true" : "false") << "\n";
        output << "REMARK   1 Include HETATM: " << (parser_config.include_hetatm ? "true" : "false") << "\n";
        output << "REMARK   1 Internal subunit identity may not be fully representable in PDB chain field\n";
        output << "REMARK   1 Export time: " << std::put_time(&tm_snapshot, "%Y-%m-%d %H:%M:%S") << "\n";
    }

    int next_serial = 1;

    for (const Chain& chain : capsid.chains()) {
        for (const Residue& residue : chain.residues()) {
            for (const Atom& atom : residue.atoms()) {
                const int serial_to_write = config.preserve_atom_serial_numbers ? atom.serial() : next_serial;
                output << formatAtomRecord(atom, serial_to_write) << '\n';

                ++stats.atoms_written;
                if (atom.isHetatm()) {
                    ++stats.hetatm_written;
                }
                if (atom.hasAltLoc()) {
                    ++stats.altloc_written;
                }
                if (atom.chainId() == ' ') {
                    ++stats.blank_chain_id_count;
                }

                ++chain_id_frequency[atom.chainId()];
                ++next_serial;
            }
        }

        if (config.emit_ter_records) {
            output << "TER\n";
        }
    }

    if (config.emit_end_record) {
        output << "END\n";
    }

    if (!output.good()) {
        throw std::runtime_error("I/O failure while writing accepted-coordinate output: " + config.output_path);
    }

    for (const auto& [chain_id, count] : chain_id_frequency) {
        if (chain_id != ' ' && count > 1) {
            ++stats.repeated_chain_id_groups;
        }
    }

    if (logger_ != nullptr) {
        logger_->info("Accepted-coordinate export completed");
        logger_->info("Exported atoms: " + std::to_string(stats.atoms_written));
        logger_->info("Exported residues: " + std::to_string(stats.residues_written));
        logger_->info("Exported subunits: " + std::to_string(stats.subunits_written));

        if (stats.blank_chain_id_count > 0) {
            logger_->warning("Accepted-coordinate export contains blank chain identifiers: " +
                             std::to_string(stats.blank_chain_id_count));
        }
        if (stats.repeated_chain_id_groups > 0) {
            logger_->warning("Accepted-coordinate export contains repeated chain identifiers. Distinct repeated IDs: " +
                             std::to_string(stats.repeated_chain_id_groups));
        }
        if (stats.altloc_written > 0) {
            logger_->warning("Accepted-coordinate export contains alternate-location atoms: " +
                             std::to_string(stats.altloc_written));
        }
        if (stats.hetatm_written > 0) {
            logger_->warning("Accepted-coordinate export contains accepted HETATM records: " +
                             std::to_string(stats.hetatm_written));
        }
    }

    return stats;
}

std::string AcceptedStructureWriter::formatAtomRecord(const Atom& atom, int serial_number) const {
    std::ostringstream oss;
    oss << std::left << std::setw(6) << (atom.isHetatm() ? "HETATM" : "ATOM");
    oss << std::right << std::setw(5) << serial_number << ' ';
    oss << std::left << std::setw(4) << trimToWidth(atom.name(), 4);
    oss << atom.altLoc();
    oss << std::left << std::setw(3) << trimToWidth(atom.residueName(), 3) << ' ';
    oss << atom.chainId();
    oss << std::right << std::setw(4) << atom.residueSeq();
    oss << atom.insertionCode() << "   ";

    oss << std::fixed << std::setprecision(3)
        << std::setw(8) << atom.x()
        << std::setw(8) << atom.y()
        << std::setw(8) << atom.z();

    oss << std::setprecision(2)
        << std::setw(6) << atom.occupancy()
        << std::setw(6) << atom.tempFactor();

    oss << "          ";
    oss << std::right << std::setw(2) << trimToWidth(atom.element(), 2);
    oss << std::setw(2) << trimToWidth(atom.charge(), 2);

    return oss.str();
}

std::string AcceptedStructureWriter::trimToWidth(const std::string& value, std::size_t width) const {
    if (value.size() <= width) {
        return value;
    }
    return value.substr(0, width);
}
