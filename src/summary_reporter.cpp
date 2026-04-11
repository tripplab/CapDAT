#include "summary_reporter.hpp"

#include <iomanip>
#include <ostream>
#include <string>

void printStructuralSummaryBlock(std::ostream& out, const StructuralSummary& summary) {
    const auto old_flags = out.flags();
    const auto old_precision = out.precision();

    out << "\n";
    out << std::fixed << std::setprecision(3);

    out << "Geometric center (accepted atoms): "
        << "(" << summary.center_x << ", "
        << summary.center_y << ", "
        << summary.center_z << ")\n";

    out << "Coordinate bounds:\n";
    out << "  x: [" << summary.x_min << ", " << summary.x_max << "]\n";
    out << "  y: [" << summary.y_min << ", " << summary.y_max << "]\n";
    out << "  z: [" << summary.z_min << ", " << summary.z_max << "]\n";

    out << "Axis spans: "
        << "x=" << summary.x_span << ", "
        << "y=" << summary.y_span << ", "
        << "z=" << summary.z_span << "\n";

    out << "Radial extent from geometric center: "
        << "r_min=" << summary.r_min << ", "
        << "r_max=" << summary.r_max << ", "
        << "r_mean=" << summary.r_mean << ", "
        << "r_stddev=" << summary.r_stddev << "\n";

    out << "Shell-thickness estimate (r_max - r_min): "
        << summary.shell_thickness_estimate << "\n";

    out << "Bounding-box center (supplementary): "
        << "(" << summary.bbox_center_x << ", "
        << summary.bbox_center_y << ", "
        << summary.bbox_center_z << ")\n";

    out << "Per-subunit counts (compact):\n";
    out << "  atoms/subunit: min=" << summary.atoms_per_subunit.min
        << ", max=" << summary.atoms_per_subunit.max
        << ", mean=" << summary.atoms_per_subunit.mean << "\n";
    out << "  residues/subunit: min=" << summary.residues_per_subunit.min
        << ", max=" << summary.residues_per_subunit.max
        << ", mean=" << summary.residues_per_subunit.mean << "\n";
    out << "  unique original labels: " << summary.unique_original_label_count << "\n";

    constexpr std::size_t kMaxInlineLabelList = 12;
    if (!summary.sorted_unique_original_labels.empty() &&
        summary.sorted_unique_original_labels.size() <= kMaxInlineLabelList) {
        out << "  labels: ";
        for (std::size_t i = 0; i < summary.sorted_unique_original_labels.size(); ++i) {
            if (i > 0) {
                out << ", ";
            }

            const char label = summary.sorted_unique_original_labels[i];
            if (label == ' ') {
                out << "<blank>";
            } else {
                out << std::string(1, label);
            }
        }
        out << "\n";
    }

    out.flags(old_flags);
    out.precision(old_precision);
}
