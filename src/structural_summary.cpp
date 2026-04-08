#include "structural_summary.hpp"

#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace {

inline void updateRange(double value, double& min_value, double& max_value) {
    if (value < min_value) {
        min_value = value;
    }
    if (value > max_value) {
        max_value = value;
    }
}

} // namespace

StructuralSummary computeStructuralSummary(const Capsid& capsid) {
    StructuralSummary summary;

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;

    summary.x_min = std::numeric_limits<double>::infinity();
    summary.y_min = std::numeric_limits<double>::infinity();
    summary.z_min = std::numeric_limits<double>::infinity();

    summary.x_max = -std::numeric_limits<double>::infinity();
    summary.y_max = -std::numeric_limits<double>::infinity();
    summary.z_max = -std::numeric_limits<double>::infinity();

    summary.internal_subunit_count = capsid.chains().size();

    std::size_t min_atoms_per_subunit = std::numeric_limits<std::size_t>::max();
    std::size_t max_atoms_per_subunit = 0;
    std::size_t sum_atoms_per_subunit = 0;

    std::size_t min_residues_per_subunit = std::numeric_limits<std::size_t>::max();
    std::size_t max_residues_per_subunit = 0;
    std::size_t sum_residues_per_subunit = 0;
    std::set<char> unique_original_labels;

    for (const Chain& chain : capsid.chains()) {
        const std::size_t chain_atom_count = chain.atomCount();
        const std::size_t chain_residue_count = chain.residueCount();
        unique_original_labels.insert(chain.pdbChainId());

        if (chain_atom_count < min_atoms_per_subunit) {
            min_atoms_per_subunit = chain_atom_count;
        }
        if (chain_atom_count > max_atoms_per_subunit) {
            max_atoms_per_subunit = chain_atom_count;
        }
        sum_atoms_per_subunit += chain_atom_count;

        if (chain_residue_count < min_residues_per_subunit) {
            min_residues_per_subunit = chain_residue_count;
        }
        if (chain_residue_count > max_residues_per_subunit) {
            max_residues_per_subunit = chain_residue_count;
        }
        sum_residues_per_subunit += chain_residue_count;

        for (const Residue& residue : chain.residues()) {
            for (const Atom& atom : residue.atoms()) {
                ++summary.accepted_atom_count;
                sum_x += atom.x();
                sum_y += atom.y();
                sum_z += atom.z();

                updateRange(atom.x(), summary.x_min, summary.x_max);
                updateRange(atom.y(), summary.y_min, summary.y_max);
                updateRange(atom.z(), summary.z_min, summary.z_max);
            }
        }
    }

    if (summary.accepted_atom_count == 0) {
        throw std::runtime_error("Cannot compute structural summary with zero accepted atoms");
    }

    const double atom_count_d = static_cast<double>(summary.accepted_atom_count);
    summary.center_x = sum_x / atom_count_d;
    summary.center_y = sum_y / atom_count_d;
    summary.center_z = sum_z / atom_count_d;

    summary.x_span = summary.x_max - summary.x_min;
    summary.y_span = summary.y_max - summary.y_min;
    summary.z_span = summary.z_max - summary.z_min;

    summary.bbox_center_x = (summary.x_min + summary.x_max) / 2.0;
    summary.bbox_center_y = (summary.y_min + summary.y_max) / 2.0;
    summary.bbox_center_z = (summary.z_min + summary.z_max) / 2.0;

    double sum_r = 0.0;
    double sum_r_squared = 0.0;

    summary.r_min = std::numeric_limits<double>::infinity();
    summary.r_max = 0.0;

    for (const Chain& chain : capsid.chains()) {
        for (const Residue& residue : chain.residues()) {
            for (const Atom& atom : residue.atoms()) {
                const double dx = atom.x() - summary.center_x;
                const double dy = atom.y() - summary.center_y;
                const double dz = atom.z() - summary.center_z;
                const double r = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));

                if (r < summary.r_min) {
                    summary.r_min = r;
                }
                if (r > summary.r_max) {
                    summary.r_max = r;
                }

                sum_r += r;
                sum_r_squared += r * r;
            }
        }
    }

    summary.r_mean = sum_r / atom_count_d;

    const double mean_r_squared = sum_r_squared / atom_count_d;
    const double variance = mean_r_squared - (summary.r_mean * summary.r_mean);
    summary.r_stddev = std::sqrt((variance > 0.0) ? variance : 0.0);

    summary.shell_thickness_estimate = summary.r_max - summary.r_min;

    if (summary.internal_subunit_count > 0) {
        const double subunit_count_d = static_cast<double>(summary.internal_subunit_count);

        summary.atoms_per_subunit.min = static_cast<double>(min_atoms_per_subunit);
        summary.atoms_per_subunit.max = static_cast<double>(max_atoms_per_subunit);
        summary.atoms_per_subunit.mean = static_cast<double>(sum_atoms_per_subunit) / subunit_count_d;

        summary.residues_per_subunit.min = static_cast<double>(min_residues_per_subunit);
        summary.residues_per_subunit.max = static_cast<double>(max_residues_per_subunit);
        summary.residues_per_subunit.mean = static_cast<double>(sum_residues_per_subunit) / subunit_count_d;
    }

    summary.unique_original_label_count = unique_original_labels.size();
    summary.sorted_unique_original_labels.assign(unique_original_labels.begin(), unique_original_labels.end());

    return summary;
}
