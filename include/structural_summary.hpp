#ifndef CAPDAT_STRUCTURAL_SUMMARY_HPP
#define CAPDAT_STRUCTURAL_SUMMARY_HPP

#include <cstddef>

#include "capsid.hpp"

struct RangeStats {
    double min = 0.0;
    double max = 0.0;
    double mean = 0.0;
};

struct StructuralSummary {
    std::size_t accepted_atom_count = 0;

    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;

    double x_min = 0.0;
    double x_max = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;
    double z_min = 0.0;
    double z_max = 0.0;

    double x_span = 0.0;
    double y_span = 0.0;
    double z_span = 0.0;

    double r_min = 0.0;
    double r_max = 0.0;
    double r_mean = 0.0;
    double r_stddev = 0.0;

    double shell_thickness_estimate = 0.0;

    double bbox_center_x = 0.0;
    double bbox_center_y = 0.0;
    double bbox_center_z = 0.0;

    std::size_t internal_subunit_count = 0;
    RangeStats atoms_per_subunit{};
    RangeStats residues_per_subunit{};
};

[[nodiscard]] StructuralSummary computeStructuralSummary(const Capsid& capsid);

#endif // CAPDAT_STRUCTURAL_SUMMARY_HPP
