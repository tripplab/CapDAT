#ifndef CAPDAT_CHAIN_HPP
#define CAPDAT_CHAIN_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "atom.hpp"
#include "residue.hpp"

/**
 * @brief Represents one internally reconstructed capsid subunit.
 *
 * In v01, the implementation class name remains Chain for convenience, but
 * conceptually this object represents a unique structural subunit, not merely
 * the raw one-letter PDB chain label.
 *
 * This distinction matters because large capsid PDB files may reuse the same
 * chain identifier across multiple independent proteins. For that reason,
 * Chain stores:
 *
 * - a unique internal subunit identifier,
 * - the original PDB chain label as metadata,
 * - an ordered collection of residues,
 * - an atom-count cache for efficient reporting.
 */
class Chain {
public:
    /**
     * @brief Default constructor.
     *
     * Useful for placeholder initialization and container compatibility during
     * early scaffolding.
     */
    Chain() = default;

    /**
     * @brief Construct a Chain with a unique internal ID and original PDB label.
     *
     * @param internal_id Unique internal identifier assigned by CapDAT.
     * @param pdb_chain_id Original PDB chain identifier as read from the file.
     */
    Chain(std::size_t internal_id, char pdb_chain_id);

    // -------------------------------------------------------------------------
    // Basic metadata accessors
    // -------------------------------------------------------------------------

    /// @return Unique internal subunit identifier.
    [[nodiscard]] std::size_t internalId() const;

    /// @return Original PDB chain identifier.
    [[nodiscard]] char pdbChainId() const;

    // -------------------------------------------------------------------------
    // Residue storage and access
    // -------------------------------------------------------------------------

    /**
     * @brief Append a residue to this Chain.
     *
     * Residues are expected to be added in parsing order. In v01, the parser is
     * responsible for determining when a new residue should be created.
     *
     * @param residue Residue object to append.
     */
    void addResidue(Residue residue);

    /**
     * @brief Append an atom to the most recently added residue.
     *
     * This method exists to keep residue mutation and atom-count cache updates
     * synchronized in one place.
     *
     * Precondition: at least one residue must already exist.
     *
     * @param atom Atom object to append to the last residue.
     */
    void addAtomToLastResidue(Atom atom);

    /**
     * @brief Access the most recently added residue.
     *
     * This overload provides read-only access for inspection/debugging.
     *
     * Precondition: at least one residue must already exist.
     *
     * @return Const reference to the last residue in the collection.
     */
    [[nodiscard]] const Residue& lastResidue() const;

    /// @return Read-only access to all residues in this Chain.
    [[nodiscard]] const std::vector<Residue>& residues() const;

    /**
     * @brief Return mutable residues for post-parse workflow transforms.
     *
     * This supports in-place coordinate rotation after parsing. It should not be
     * used to rebuild topology or change residue identity metadata.
     */
    std::vector<Residue>& mutableResidues();

    /// @return Number of residues in this Chain.
    [[nodiscard]] std::size_t residueCount() const;

    /**
     * @brief Return the total atom count for this Chain.
     *
     * In v01, this is maintained through a lightweight cache updated whenever
     * atoms are appended through the Chain interface.
     *
     * @return Total number of atoms in this Chain.
     */
    [[nodiscard]] std::size_t atomCount() const;

private:
    // -------------------------------------------------------------------------
    // Chain identity metadata
    // -------------------------------------------------------------------------

    std::size_t internal_id_ = 0;
    char pdb_chain_id_ = ' ';

    // -------------------------------------------------------------------------
    // Contained residues
    // -------------------------------------------------------------------------

    std::vector<Residue> residues_;

    // -------------------------------------------------------------------------
    // Cached summary data
    // -------------------------------------------------------------------------

    std::size_t atom_count_cache_ = 0;
};

// NOTE ON CONTROLLED MUTATION:
//
// This version intentionally routes atom appends through
// addAtomToLastResidue() instead of exposing a mutable lastResidue() to
// external code. The reason is to keep the chain-level atom-count cache
// consistent automatically.
//
// This slightly narrows the mutating API, but it is a worthwhile tradeoff in
// v01 because it reduces the risk of subtle counter bugs during incremental
// parsing.

// NOTE ON CLASS NAMING:
//
// In v01, the class name Chain is retained for implementation convenience and
// because it maps naturally to familiar structural hierarchy terminology.
// However, this object should be understood as an internally reconstructed
// subunit, not as a guarantee of one unique raw PDB chain label.
//
// A future revision could rename this class to something like Subunit if that
// improves conceptual clarity across later analytical modules. For the
// foundation release, keeping the name Chain is a practical compromise between
// clarity, implementation simplicity, and continuity with the design document.

#endif // CAPDAT_CHAIN_HPP
