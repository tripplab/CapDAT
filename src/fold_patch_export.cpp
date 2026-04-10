#include "fold_patch_export.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "fold_patch_analysis.hpp"

namespace {

std::string trimToWidth(const std::string& value, std::size_t width) {
    if (value.size() <= width) {
        return value;
    }
    return value.substr(0, width);
}

int parseResidueSeq(const std::string& residue_id) {
    const std::size_t first_colon = residue_id.find(':');
    if (first_colon == std::string::npos) {
        return 0;
    }
    const std::size_t second_colon = residue_id.find(':', first_colon + 1);
    const std::string seq_text =
        second_colon == std::string::npos
            ? residue_id.substr(first_colon + 1)
            : residue_id.substr(first_colon + 1, second_colon - first_colon - 1);

    try {
        return std::stoi(seq_text);
    } catch (const std::exception&) {
        return 0;
    }
}

char chainIdFromSubunitId(const std::string& subunit_id) {
    try {
        const int raw = std::stoi(subunit_id);
        const int value = raw >= 0 ? raw : -raw;
        if (value < 26) {
            return static_cast<char>('A' + value);
        }
        return static_cast<char>('0' + (value % 10));
    } catch (const std::exception&) {
        return 'A';
    }
}

std::string formatPatchAtomRecord(const PatchAtom& atom) {
    std::ostringstream oss;

    oss << std::left << std::setw(6) << "ATOM";
    oss << std::right << std::setw(5) << atom.serial << ' ';
    oss << std::left << std::setw(4) << trimToWidth(atom.atom_name, 4);
    oss << ' ';
    oss << std::left << std::setw(3) << "UNK" << ' ';
    oss << chainIdFromSubunitId(atom.subunit_id);
    oss << std::right << std::setw(4) << parseResidueSeq(atom.residue_id);
    oss << ' ' << "   ";

    oss << std::fixed << std::setprecision(3)
        << std::setw(8) << atom.x
        << std::setw(8) << atom.y
        << std::setw(8) << atom.z;

    oss << std::setprecision(2)
        << std::setw(6) << 1.00
        << std::setw(6) << 0.00;
    oss << "          ";
    oss << std::right << std::setw(2) << trimToWidth(atom.element, 2);

    std::string line = oss.str();
    if (line.size() < 80) {
        line.append(80 - line.size(), ' ');
    } else if (line.size() > 80) {
        line.resize(80);
    }
    return line;
}

} // namespace

void exportFoldPatch(const FoldPatchAnalysisResult& result,
                     const FoldPatchAnalysisConfig& config) {
    if (!config.export_patch_pdb) {
        return;
    }

    const std::string output_path = config.out_prefix + "_patch.pdb";
    std::ofstream output(output_path, std::ios::out | std::ios::trunc);
    if (!output.is_open()) {
        std::cerr << "[fold_patch] WARNING: failed to open patch PDB output: " << output_path << '\n';
        return;
    }

    for (const PatchAtom& atom : result.patch_atoms) {
        output << formatPatchAtomRecord(atom) << '\n';
    }
    output << "END\n";

    if (!output.good()) {
        std::cerr << "[fold_patch] WARNING: failed while writing patch PDB output: " << output_path << '\n';
        return;
    }

    const_cast<FoldPatchAnalysisResult&>(result).exported_files.push_back(output_path);
}
