#include "atom.hpp"

#include <utility>

/**
 * @brief Construct a fully initialized Atom.
 *
 * String-like parameters are accepted by value and then moved into the data
 * members. This keeps the public API simple while still allowing efficient
 * construction from parser-produced temporaries.
 */
Atom::Atom(int serial,
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
           bool is_hetatm)
    : serial_(serial),
      name_(std::move(name)),
      alt_loc_(alt_loc),
      residue_name_(std::move(residue_name)),
      chain_id_(chain_id),
      residue_seq_(residue_seq),
      insertion_code_(insertion_code),
      x_(x),
      y_(y),
      z_(z),
      occupancy_(occupancy),
      temp_factor_(temp_factor),
      element_(std::move(element)),
      charge_(std::move(charge)),
      is_hetatm_(is_hetatm) {
}

/**
 * @brief Return the PDB atom serial number.
 */
int Atom::serial() const {
    return serial_;
}

/**
 * @brief Return the atom name.
 */
const std::string& Atom::name() const {
    return name_;
}

/**
 * @brief Return the alternate-location indicator.
 *
 * A blank character means no alternate location was explicitly assigned.
 */
char Atom::altLoc() const {
    return alt_loc_;
}

/**
 * @brief Return the residue name.
 */
const std::string& Atom::residueName() const {
    return residue_name_;
}

/**
 * @brief Return the original PDB chain identifier.
 */
char Atom::chainId() const {
    return chain_id_;
}

/**
 * @brief Return the residue sequence number.
 */
int Atom::residueSeq() const {
    return residue_seq_;
}

/**
 * @brief Return the PDB insertion code.
 *
 * A blank character means no insertion code was explicitly assigned.
 */
char Atom::insertionCode() const {
    return insertion_code_;
}

/**
 * @brief Return the X coordinate.
 */
double Atom::x() const {
    return x_;
}

/**
 * @brief Return the Y coordinate.
 */
double Atom::y() const {
    return y_;
}

/**
 * @brief Return the Z coordinate.
 */
double Atom::z() const {
    return z_;
}

/**
 * @brief Return the occupancy value.
 */
double Atom::occupancy() const {
    return occupancy_;
}

/**
 * @brief Return the temperature factor (B-factor).
 */
double Atom::tempFactor() const {
    return temp_factor_;
}

/**
 * @brief Return the element symbol.
 */
const std::string& Atom::element() const {
    return element_;
}

/**
 * @brief Return the charge field.
 */
const std::string& Atom::charge() const {
    return charge_;
}

/**
 * @brief Return whether this atom came from a HETATM record.
 */
bool Atom::isHetatm() const {
    return is_hetatm_;
}

/**
 * @brief Return whether this atom has a non-blank alternate-location code.
 */
bool Atom::hasAltLoc() const {
    return alt_loc_ != ' ';
}

/**
 * @brief Return the Cartesian position as {x, y, z}.
 */
std::array<double, 3> Atom::position() const {
    return {x_, y_, z_};
}

// NOTE ON PASS-BY-VALUE + std::move:
//
// The constructor takes string-like fields by value and then moves them into
// the data members. This is a practical v01 choice: it keeps the public API
// simple while still allowing efficient construction from temporary strings
// produced by the parser.
//
// If profiling ever shows that construction overhead matters significantly, we
// could revisit argument-passing strategy. For the foundation release, this is
// a clean and modern compromise.

// NOTE ON METHOD SIMPLICITY:
//
// The implementation intentionally keeps all accessors trivial and avoids
// embedding higher-level logic in Atom. This follows the v01 design principle
// that Atom should behave primarily as a lightweight data holder.
//
// If future analytical modules need derived geometric or chemical behavior,
// those additions should be introduced carefully to avoid turning Atom into an
// overly heavy domain class too early.
