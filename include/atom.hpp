#ifndef CAPDAT_ATOM_HPP
#define CAPDAT_ATOM_HPP

#include <array>
#include <string>

/**
 * @brief Represents a single atomic coordinate record from a PDB file.
 *
 * In v01, Atom is intentionally kept as a lightweight data holder.
 * Its main responsibility is to store atomic-level information parsed
 * from ATOM or HETATM records and provide read-only access to that data.
 *
 * Complex geometric or chemical logic is intentionally deferred to
 * future versions of CapDAT.
 */
class Atom {
public:
    /**
     * @brief Default constructor.
     *
     * This is useful for containers and placeholder initialization during
     * early scaffolding. The resulting object is valid but empty/defaulted.
     */
    Atom() = default;

    /**
     * @brief Construct a fully initialized Atom.
     *
     * @param serial Atom serial number from the PDB record.
     * @param name Atom name (e.g. "CA", "N", "O").
     * @param alt_loc Alternate-location indicator.
     * @param residue_name Residue name (e.g. "ALA", "GLY").
     * @param chain_id Original PDB chain identifier as read from file.
     * @param residue_seq Residue sequence number.
     * @param insertion_code PDB insertion code.
     * @param x X Cartesian coordinate.
     * @param y Y Cartesian coordinate.
     * @param z Z Cartesian coordinate.
     * @param occupancy Occupancy value, if available.
     * @param temp_factor Temperature factor (B-factor), if available.
     * @param element Element symbol, if available.
     * @param charge Formal charge field, if available.
     * @param is_hetatm True if the source record was HETATM, false if ATOM.
     */
    Atom(int serial,
         std::string name,
         char alt_loc,
         std::string residue_name,
         char chain_id,
         int residue_seq,
         char insertion_code,
         double x,
         double y,
         double z,
         double occupancy,
         double temp_factor,
         std::string element,
         std::string charge,
         bool is_hetatm);

    // -------------------------------------------------------------------------
    // Basic field accessors
    // -------------------------------------------------------------------------

    /// @return PDB atom serial number.
    [[nodiscard]] int serial() const;

    /// @return Atom name.
    [[nodiscard]] const std::string& name() const;

    /// @return Alternate-location indicator, or blank if absent.
    [[nodiscard]] char altLoc() const;

    /// @return Residue name.
    [[nodiscard]] const std::string& residueName() const;

    /// @return Original PDB chain identifier.
    [[nodiscard]] char chainId() const;

    /// @return Residue sequence number.
    [[nodiscard]] int residueSeq() const;

    /// @return PDB insertion code, or blank if absent.
    [[nodiscard]] char insertionCode() const;

    /// @return X coordinate.
    [[nodiscard]] double x() const;

    /// @return Y coordinate.
    [[nodiscard]] double y() const;

    /// @return Z coordinate.
    [[nodiscard]] double z() const;

    /// @return Occupancy value.
    [[nodiscard]] double occupancy() const;

    /// @return Temperature factor (B-factor).
    [[nodiscard]] double tempFactor() const;

    /// @return Element symbol.
    [[nodiscard]] const std::string& element() const;

    /// @return Charge string.
    [[nodiscard]] const std::string& charge() const;

    // -------------------------------------------------------------------------
    // Convenience methods
    // -------------------------------------------------------------------------

    /// @return True if this atom came from a HETATM record.
    [[nodiscard]] bool isHetatm() const;

    /**
     * @brief Indicates whether this atom has an alternate-location code.
     *
     * In standard PDB usage, a blank altLoc means no alternate location
     * is explicitly assigned.
     *
     * @return True if altLoc is not blank.
     */
    [[nodiscard]] bool hasAltLoc() const;

    /**
     * @brief Return Cartesian coordinates as a compact array.
     *
     * This is convenient for future geometric operations while keeping
     * the Atom API simple.
     *
     * @return {x, y, z}
     */
    [[nodiscard]] std::array<double, 3> position() const;

private:
    // -------------------------------------------------------------------------
    // Stored PDB-derived fields
    // -------------------------------------------------------------------------

    int serial_ = 0;
    std::string name_;
    char alt_loc_ = ' ';
    std::string residue_name_;
    char chain_id_ = ' ';
    int residue_seq_ = 0;
    char insertion_code_ = ' ';
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    double occupancy_ = 0.0;
    double temp_factor_ = 0.0;
    std::string element_;
    std::string charge_;
    bool is_hetatm_ = false;
};

#endif // CAPDAT_ATOM_HPP


// NOTE ON STRING-BASED TEXT FIELDS:
//
// In v01, fields such as name_, residue_name_, element_, and charge_ are kept
// as std::string on purpose. This favors correctness, readability, simplicity,
// and easier debugging during the foundation stage of the project.
//
// A future performance-oriented revision may revisit some of these fields if
// profiling shows that parser throughput, memory footprint, or cache behavior
// becomes a real bottleneck for very large structures or high-throughput runs.
// In that case, possible alternatives could include fixed-width character
// arrays, compact integer/enumerated encodings, or shared/interned strings.
//
// For v01, the added clarity and maintainability of std::string are considered
// more important than premature micro-optimization.



