#include "fold_patch_extraction.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <unordered_map>

#include "fold_patch_analysis.hpp"

namespace {

std::string trim(const std::string& value) {
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    }).base();

    if (first >= last) {
        return {};
    }
    return std::string(first, last);
}

std::string inferElementFromAtomName(const std::string& atom_name) {
    for (char ch : atom_name) {
        if (std::isalpha(static_cast<unsigned char>(ch)) != 0) {
            return std::string(1, static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
        }
    }
    return {};
}

std::string normalizeElementSymbol(const std::string& value) {
    std::string normalized = trim(value);
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return normalized;
}

} // namespace

std::vector<PatchAtom> extractPatch(const Capsid& capsid,
                                    const LocalFrame& frame,
                                    const FoldPatchAnalysisConfig& config) {
    std::vector<PatchAtom> patch;

    for (const Chain& chain : capsid.chains()) {
        for (const Residue& residue : chain.residues()) {
            for (const Atom& atom : residue.atoms()) {
                if (atom.isHetatm()) {
                    continue;
                }

                const double dx = atom.x() - frame.origin.x;
                const double dy = atom.y() - frame.origin.y;
                const double dz = atom.z() - frame.origin.z;

                const double lx = frame.axes.m[0][0] * dx + frame.axes.m[1][0] * dy + frame.axes.m[2][0] * dz;
                const double ly = frame.axes.m[0][1] * dx + frame.axes.m[1][1] * dy + frame.axes.m[2][1] * dz;
                const double lz = frame.axes.m[0][2] * dx + frame.axes.m[1][2] * dy + frame.axes.m[2][2] * dz;

                if ((lx * lx + ly * ly) > (config.R * config.R)) {
                    continue;
                }
                if (lz <= 0.0 || lz > config.z_margin) {
                    continue;
                }

                PatchAtom pa;
                pa.x = lx;
                pa.y = ly;
                pa.z = lz;
                pa.serial = atom.serial();
                pa.atom_name = atom.name();
                pa.element = normalizeElementSymbol(atom.element());
                if (pa.element.empty()) {
                    pa.element = inferElementFromAtomName(pa.atom_name);
                }
                pa.vdw_radius = vdwRadius(pa.element);
                pa.residue_id = residue.residueKey();
                pa.subunit_id = std::to_string(residue.internalSubunitId());

                patch.push_back(std::move(pa));
            }
        }
    }

    return patch;
}

double vdwRadius(const std::string& element) {
    static const std::unordered_map<std::string, double> kTable = {
        {"H", 1.20},
        {"C", 1.70},
        {"N", 1.55},
        {"O", 1.52},
        {"S", 1.80},
        {"P", 1.80},
    };

    const std::string normalized = normalizeElementSymbol(element);
    const auto it = kTable.find(normalized);
    if (it != kTable.end()) {
        return it->second;
    }

    std::cerr << "[fold_patch] WARNING: unknown element '" << element
              << "', using fallback vdW = 1.70 Å\n";
    return 1.70;
}
