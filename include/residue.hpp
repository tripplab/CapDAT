#ifndef CAPDAT_RESIDUE_HPP
#define CAPDAT_RESIDUE_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "atom.hpp"

/**
 * @brief Represents a residue within one internally reconstructed capsid subunit.
 *
 * A Residue groups the atoms that belong to the same residue identity inside a
 * single internal subunit. In v01, residue identity is tracked using:
 *
 * - residue name,
 * - residue sequence number,
 * - insertion code,
 * - original PDB chain label as metadata,
 * - internal subunit identifier.
 *
 * The internal subunit identifier is important because the raw one-letter PDB
 * chain label is not guaranteed to be globally unique in large capsid
 * assemblies.
 */
class Residue {
public:
    /**
     * @brief Default constructor.
     *
     * Useful for placeholder initialization and container compatibility during
     * early scaffolding.
     */
    Residue() = default;

    /**
     * @brief Construct a residue with its identifying metadata.
     *
     * @param name Residue name (e.g. "ALA", "GLY").
     * @param seq_number Residue sequence number from the PDB record.
     * @param insertion_code PDB insertion code, or blank if absent.
     * @param pdb_chain_id Original PDB chain identifier as read from file.
     * @param internal_subunit_id Unique internal subunit identifier assigned by CapDAT.
     */
    Residue(std::string name,
            int seq_number,
            char insertion_code,
            char pdb_chain_id,
            std::size_t internal_subunit_id);

    // -------------------------------------------------------------------------
    // Basic metadata accessors
    // -------------------------------------------------------------------------

    /// @return Residue name.
    [[nodiscard]] const std::string& name() const;

    /// @return Residue sequence number.
    [[nodiscard]] int seqNumber() const;

    /// @return PDB insertion code, or blank if absent.
    [[nodiscard]] char insertionCode() const;

    /// @return Original PDB chain identifier.
    [[nodiscard]] char pdbChainId() const;

    /// @return Unique internal subunit identifier.
    [[nodiscard]] std::size_t internalSubunitId() const;

    // -------------------------------------------------------------------------
    // Atom storage and access
    // -------------------------------------------------------------------------

    /**
     * @brief Append an atom to this residue.
     *
     * In v01, no advanced chemical validation is performed here. The parser is
     * responsible for deciding whether the atom belongs to this residue before
     * calling addAtom().
     *
     * @param atom Atom object to append.
     */
    void addAtom(Atom atom);

    /// @return Number of atoms currently stored in this residue.
    [[nodiscard]] std::size_t atomCount() const;

    /// @return Read-only access to the atoms in this residue.
    [[nodiscard]] const std::vector<Atom>& atoms() const;

    /**
     * @brief Return mutable access to atoms for post-parse in-place workflows.
     *
     * Parser construction should continue to use addAtom(). This mutable view is
     * reserved for explicit workflow-stage coordinate transformations.
     */
    std::vector<Atom>& mutableAtoms();

    /**
     * @brief Build a compact residue key for reporting/debugging.
     *
     * The intended v01 logic is based primarily on:
     * - internal subunit identifier,
     * - residue sequence number,
     * - insertion code.
     *
     * Residue name is kept as associated metadata but is not intended to be the
     * primary grouping key.
     *
     * @return A string key representing this residue identity.
     */
    [[nodiscard]] std::string residueKey() const;

private:
    // -------------------------------------------------------------------------
    // Residue identity metadata
    // -------------------------------------------------------------------------

    std::string name_;
    int seq_number_ = 0;
    char insertion_code_ = ' ';
    char pdb_chain_id_ = ' ';
    std::size_t internal_subunit_id_ = 0;

    // -------------------------------------------------------------------------
    // Contained atoms
    // -------------------------------------------------------------------------

    std::vector<Atom> atoms_;
};

// NOTE ON CONTAINER CHOICE:
//
// In v01, Residue stores atoms in std::vector<Atom>. This is intentional:
// it keeps ownership simple, preserves insertion order naturally, and provides
// efficient contiguous storage for iteration.
//
// A future performance-oriented revision could revisit this if profiling shows
// that a different layout is beneficial for memory locality, indexing, or bulk
// geometry operations. For the foundation release, std::vector offers the best
// balance of clarity, correctness, and maintainability.

#endif // CAPDAT_RESIDUE_HPP
