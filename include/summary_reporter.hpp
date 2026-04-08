#ifndef CAPDAT_SUMMARY_REPORTER_HPP
#define CAPDAT_SUMMARY_REPORTER_HPP

#include <iosfwd>

#include "structural_summary.hpp"

void printStructuralSummaryBlock(std::ostream& out, const StructuralSummary& summary);

#endif // CAPDAT_SUMMARY_REPORTER_HPP
