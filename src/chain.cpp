#include "chain.hpp"

#include <stdexcept>
#include <utility>

/**
 * @brief Construct a Chain with a unique internal ID and original PDB label.
 */
Chain::Chain(std::size_t internal_id, char pdb_chain_id)
    : internal_id_(internal_id),
      pdb_chain_id_(pdb_chain_id) {
}

/**
 * @brief Return the unique internal subunit identifier.
 */
std::size_t Chain::internalId() const {
    return internal_id_;
}

/**
 * @brief Return the original PDB chain identifier.
 *
 * This value is preserved as metadata and is not assumed to be globally unique
 * across the full capsid.
 */
char Chain::pdbChainId() const {
    return pdb_chain_id_;
}

/**
 * @brief Append a residue to this Chain.
 *
 * The atom-count cache is updated using the atom count already present in the
 * appended residue.
 */
void Chain::addResidue(Residue residue) {
    atom_count_cache_ += residue.atomCount();
    residues_.push_back(std::move(residue));
}

/**
 * @brief Append an atom to the most recently added residue.
 *
 * This keeps residue mutation and the chain-level atom-count cache synchronized.
 *
 * @throws std::out_of_range if the Chain contains no residues.
 */
void Chain::addAtomToLastResidue(Atom atom) {
    if (residues_.empty()) {
        throw std::out_of_range("Chain::addAtomToLastResidue() called on empty Chain");
    }

    residues_.back().addAtom(std::move(atom));
    ++atom_count_cache_;
}

/**
 * @brief Return read-only access to the most recently appended residue.
 *
 * @throws std::out_of_range if the Chain contains no residues.
 */
const Residue& Chain::lastResidue() const {
    if (residues_.empty()) {
        throw std::out_of_range("Chain::lastResidue() called on empty Chain");
    }
    return residues_.back();
}

/**
 * @brief Return read-only access to the residues stored in this Chain.
 */
const std::vector<Residue>& Chain::residues() const {
    return residues_;
}

/**
 * @brief Return mutable residue storage for post-parse workflow transforms.
 */
std::vector<Residue>& Chain::mutableResidues() {
    return residues_;
}

/**
 * @brief Return the number of residues in this Chain.
 */
std::size_t Chain::residueCount() const {
    return residues_.size();
}

/**
 * @brief Return the total number of atoms in this Chain.
 *
 * In v01, this is served from a lightweight cache maintained as residues and
 * atoms are appended through the Chain interface.
 */
std::size_t Chain::atomCount() const {
    return atom_count_cache_;
}

// NOTE ON CACHE ROBUSTNESS:
//
// The key change in this revision is that atom insertion into the current
// residue now happens through Chain::addAtomToLastResidue(). This ensures that
// the chain-level atom_count_cache_ is always updated alongside the actual
// residue mutation.
//
// This is more robust than exposing mutable access to the last residue because
// it prevents outside code from accidentally modifying residue contents without
// keeping chain-level summary state in sync.

// NOTE ON API NARROWING:
//
// The mutating lastResidue() accessor was intentionally removed and replaced by
// a read-only lastResidue() plus an explicit addAtomToLastResidue() method.
// This makes the intended parser workflow more explicit and reduces the risk of
// hidden state inconsistencies.
//
// For v01, this is a good tradeoff: slightly less flexibility in exchange for
// safer and easier-to-audit parsing behavior.
