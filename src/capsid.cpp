#include "capsid.hpp"

#include <stdexcept>
#include <utility>

/**
 * @brief Construct a Capsid associated with a source file path.
 */
Capsid::Capsid(std::string source_path)
    : source_path_(std::move(source_path)) {
}

/**
 * @brief Set the source path associated with this Capsid.
 */
void Capsid::setSourcePath(std::string source_path) {
    source_path_ = std::move(source_path);
}

/**
 * @brief Return the source file path associated with this Capsid.
 */
const std::string& Capsid::sourcePath() const {
    return source_path_;
}

/**
 * @brief Append a reconstructed subunit to the Capsid.
 */
void Capsid::addChain(Chain chain) {
    chains_.push_back(std::move(chain));
}

/**
 * @brief Return read-only access to all reconstructed subunits.
 */
const std::vector<Chain>& Capsid::chains() const {
    return chains_;
}

/**
 * @brief Return mutable access to the most recently appended Chain.
 *
 * @throws std::out_of_range if the Capsid contains no chains.
 */
Chain& Capsid::lastChain() {
    if (chains_.empty()) {
        throw std::out_of_range("Capsid::lastChain() called on empty Capsid");
    }
    return chains_.back();
}

/**
 * @brief Append an atom to the last residue of the last Chain.
 *
 * @throws std::out_of_range if the Capsid contains no chains.
 */
void Capsid::addAtomToLastChainResidue(Atom atom) {
    if (chains_.empty()) {
        throw std::out_of_range("Capsid::addAtomToLastChainResidue() called on empty Capsid");
    }
    chains_.back().addAtomToLastResidue(std::move(atom));
}

/**
 * @brief Return the total number of reconstructed internal subunits.
 */
std::size_t Capsid::subunitCount() const {
    return total_subunits_;
}

/**
 * @brief Return the total number of atoms in the Capsid.
 */
std::size_t Capsid::atomCount() const {
    return total_atoms_;
}

/**
 * @brief Return the total number of residues in the Capsid.
 */
std::size_t Capsid::residueCount() const {
    return total_residues_;
}

/**
 * @brief Return the total number of accepted HETATM records.
 */
std::size_t Capsid::hetatmCount() const {
    return total_hetatm_;
}

/**
 * @brief Return the total number of alternate-location atoms encountered.
 */
std::size_t Capsid::altLocCount() const {
    return total_altloc_;
}

/**
 * @brief Return the total number of skipped or malformed records.
 */
std::size_t Capsid::skippedRecordCount() const {
    return total_skipped_records_;
}

/**
 * @brief Increment total atom count by one.
 */
void Capsid::incrementAtomCount() {
    ++total_atoms_;
}

/**
 * @brief Increment total residue count by one.
 */
void Capsid::incrementResidueCount() {
    ++total_residues_;
}

/**
 * @brief Increment total subunit count by one.
 */
void Capsid::incrementSubunitCount() {
    ++total_subunits_;
}

/**
 * @brief Increment total HETATM count by one.
 */
void Capsid::incrementHetatmCount() {
    ++total_hetatm_;
}

/**
 * @brief Increment total alternate-location atom count by one.
 */
void Capsid::incrementAltLocCount() {
    ++total_altloc_;
}

/**
 * @brief Increment total skipped/malformed record count by one.
 */
void Capsid::incrementSkippedRecordCount() {
    ++total_skipped_records_;
}

/**
 * @brief Recompute top-level counts from the stored hierarchy.
 *
 * In v01, this method provides a clean and explicit consistency point for
 * counters derived from the actual stored structure. It recomputes:
 * - total subunits,
 * - total residues,
 * - total atoms.
 *
 * Counters such as accepted HETATM records, alternate-location counts, and
 * skipped/malformed records are parser-observed event counts, so they are not
 * recomputed from the hierarchy here.
 */
void Capsid::finalizeCounts() {
    total_subunits_ = chains_.size();
    total_residues_ = 0;
    total_atoms_ = 0;

    for (const Chain& chain : chains_) {
        total_residues_ += chain.residueCount();
        total_atoms_ += chain.atomCount();
    }
}

// NOTE ON PARSER-SAFE MUTATION:
//
// lastChain() and addAtomToLastChainResidue() were added so the parser can
// continue single-pass hierarchy construction without breaking const-correctness
// or resorting to const_cast.
//
// This is a cleaner v01 design because mutation now passes through explicit
// Capsid/Chain APIs rather than through back-door access to stored containers.

// NOTE ON BOUNDARY CHECKING:
//
// Both parser-facing mutation helpers fail explicitly on empty Capsid state.
// This is deliberate defensive programming for v01: internal misuse should be
// obvious and fail clearly instead of silently producing undefined behavior.

// NOTE ON COUNT FINALIZATION:
//
// finalizeCounts() intentionally recomputes the hierarchy-derived totals from
// the stored Capsid contents instead of trusting only incremental updates made
// during parsing. This gives v01 a simple built-in consistency checkpoint and
// reduces the risk of silent counter drift.
//
// Parser-observed counters such as skipped records, altloc counts, and accepted
// HETATM events remain incremental because they are not purely properties of
// the final stored hierarchy.
