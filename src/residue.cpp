#include "residue.hpp"

#include <sstream>
#include <utility>

/**
 * @brief Construct a Residue with its identifying metadata.
 *
 * The atom container starts empty and is populated incrementally by the parser.
 */
Residue::Residue(std::string name,
                 int seq_number,
                 char insertion_code,
                 char pdb_chain_id,
                 std::size_t internal_subunit_id)
    : name_(std::move(name)),
      seq_number_(seq_number),
      insertion_code_(insertion_code),
      pdb_chain_id_(pdb_chain_id),
      internal_subunit_id_(internal_subunit_id) {
}

/**
 * @brief Return the residue name.
 */
const std::string& Residue::name() const {
    return name_;
}

/**
 * @brief Return the residue sequence number.
 */
int Residue::seqNumber() const {
    return seq_number_;
}

/**
 * @brief Return the insertion code.
 *
 * A blank character means no insertion code was explicitly assigned.
 */
char Residue::insertionCode() const {
    return insertion_code_;
}

/**
 * @brief Return the original PDB chain identifier as metadata.
 */
char Residue::pdbChainId() const {
    return pdb_chain_id_;
}

/**
 * @brief Return the internal subunit identifier.
 */
std::size_t Residue::internalSubunitId() const {
    return internal_subunit_id_;
}

/**
 * @brief Append an atom to this residue.
 *
 * In v01, membership validation is expected to be handled by the parser before
 * calling this method.
 */
void Residue::addAtom(Atom atom) {
    atoms_.push_back(std::move(atom));
}

/**
 * @brief Return the number of atoms in this residue.
 */
std::size_t Residue::atomCount() const {
    return atoms_.size();
}

/**
 * @brief Return read-only access to the residue's atom collection.
 */
const std::vector<Atom>& Residue::atoms() const {
    return atoms_;
}

/**
 * @brief Build and return a compact residue identity key.
 *
 * The v01 grouping logic is based primarily on:
 * - internal subunit identifier,
 * - residue sequence number,
 * - insertion code.
 *
 * This method is mainly useful for debugging, reporting, and future validation.
 */
std::string Residue::residueKey() const {
    std::ostringstream oss;
    oss << internal_subunit_id_ << ':'
        << seq_number_ << ':'
        << insertion_code_;
    return oss.str();
}

// NOTE ON RESIDUE KEY FORMAT:
//
// The residueKey() string is intentionally simple and primarily intended for
// debugging, reporting, and light validation in v01. It should not be treated
// as a deeply standardized serialization format.
//
// If later modules require a more formal residue identifier, we may want a
// stronger typed key or a dedicated identity structure instead of a compact
// string representation.

// NOTE ON ATOM STORAGE:
//
// Residue stores atoms in insertion order using std::vector<Atom>. This is a
// natural fit for single-pass parsing and preserves the original local order of
// records as encountered in the input.
//
// If future performance analysis suggests a different memory layout, that can
// be revisited later. For the foundation release, this is the clearest and
// most maintainable choice.
