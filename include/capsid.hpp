#ifndef CAPDAT_CAPSID_HPP
#define CAPDAT_CAPSID_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "atom.hpp"
#include "chain.hpp"

/**
 * @brief Represents the complete capsid assembly parsed from one input file.
 *
 * Capsid is the top-level domain object in the v01 hierarchy. It owns the
 * collection of internally reconstructed subunits (implemented as Chain
 * objects in v01) and stores assembly-level metadata and summary counters.
 *
 * Its main responsibilities in the foundation release are:
 *
 * - preserve the full parsed structural hierarchy,
 * - retain source-file metadata,
 * - expose whole-assembly access,
 * - provide summary counts for reporting.
 *
 * Advanced scientific analysis is intentionally out of scope for this class in
 * v01.
 */
class Capsid {
public:
    /**
     * @brief Explicit source mode used for the most recent in-place reorientation.
     */
    enum class OrientationSourceMode {
        none,
        fold,
        custom_vector
    };

    /**
     * @brief Explicit frame/orientation state for current in-memory coordinates.
     *
     * Important convention:
     * - The parser initializes Capsid in the original parsed input frame.
     * - Workflow operations may later rotate coordinates in place.
     * - Once reoriented, this object no longer represents the original input
     *   coordinate frame.
     *
     * Downstream code that cares about frame identity must consult this state
     * instead of assuming coordinates are always still in the source frame.
     */
    struct OrientationState {
        bool in_original_parsed_frame = true;
        bool reoriented_in_place = false;
        bool already_aligned_identity = false;
        std::array<std::array<double, 3>, 3> applied_rotation_matrix = {{{1.0, 0.0, 0.0},
                                                                          {0.0, 1.0, 0.0},
                                                                          {0.0, 0.0, 1.0}}};
        bool has_rotation_axis = false;
        std::array<double, 3> rotation_axis = {0.0, 0.0, 1.0};
        bool has_rotation_angle = false;
        double rotation_angle_radians = 0.0;
        OrientationSourceMode source_mode = OrientationSourceMode::none;
        std::string source_description;
        std::array<double, 3> source_direction = {0.0, 0.0, 1.0};
        char requested_target_axis = 'z';
        std::array<double, 3> target_direction = {0.0, 0.0, 1.0};
    };

    /**
     * @brief Default constructor.
     */
    Capsid() = default;

    /**
     * @brief Construct a Capsid associated with a source file path.
     *
     * @param source_path Path to the input structure file.
     */
    explicit Capsid(std::string source_path);

    // -------------------------------------------------------------------------
    // Source metadata
    // -------------------------------------------------------------------------

    /**
     * @brief Set the source path associated with this Capsid.
     *
     * @param source_path Input file path.
     */
    void setSourcePath(std::string source_path);

    /// @return Source file path associated with this Capsid.
    [[nodiscard]] const std::string& sourcePath() const;

    // -------------------------------------------------------------------------
    // Structural hierarchy ownership
    // -------------------------------------------------------------------------

    /**
     * @brief Append a reconstructed subunit to the Capsid.
     *
     * @param chain Chain object to append.
     */
    void addChain(Chain chain);

    /// @return Read-only access to all reconstructed subunits.
    [[nodiscard]] const std::vector<Chain>& chains() const;

    /// @return Mutable access for workflow-stage in-place coordinate updates.
    std::vector<Chain>& mutableChains();

    /**
     * @brief Return mutable access to the most recently appended Chain.
     *
     * This method exists primarily so the parser can continue single-pass
     * hierarchy construction without breaking const-correctness.
     *
     * Precondition: at least one Chain must already exist.
     *
     * @return Mutable reference to the last stored Chain.
     */
    Chain& lastChain();

    /**
     * @brief Append an atom to the last residue of the last Chain.
     *
     * This is a convenience pass-through used by the parser to keep mutation
     * localized and cache updates consistent through the Chain interface.
     *
     * Precondition:
     * - at least one Chain exists,
     * - the last Chain already contains at least one Residue.
     *
     * @param atom Atom object to append.
     */
    void addAtomToLastChainResidue(Atom atom);

    // -------------------------------------------------------------------------
    // Global summary counters
    // -------------------------------------------------------------------------

    /// @return Total number of reconstructed internal subunits.
    [[nodiscard]] std::size_t subunitCount() const;

    /// @return Total number of atoms in the Capsid.
    [[nodiscard]] std::size_t atomCount() const;

    /// @return Total number of residues in the Capsid.
    [[nodiscard]] std::size_t residueCount() const;

    /// @return Total number of accepted HETATM records.
    [[nodiscard]] std::size_t hetatmCount() const;

    /// @return Total number of alternate-location atoms encountered.
    [[nodiscard]] std::size_t altLocCount() const;

    /// @return Total number of skipped or malformed records.
    [[nodiscard]] std::size_t skippedRecordCount() const;

    // -------------------------------------------------------------------------
    // Incremental counter updates
    // -------------------------------------------------------------------------

    /// @brief Increment total atom count by one.
    void incrementAtomCount();

    /// @brief Increment total residue count by one.
    void incrementResidueCount();

    /// @brief Increment total subunit count by one.
    void incrementSubunitCount();

    /// @brief Increment total HETATM count by one.
    void incrementHetatmCount();

    /// @brief Increment total alternate-location atom count by one.
    void incrementAltLocCount();

    /// @brief Increment total skipped/malformed record count by one.
    void incrementSkippedRecordCount();

    // -------------------------------------------------------------------------
    // Finalization
    // -------------------------------------------------------------------------

    /**
     * @brief Recompute or validate summary counters from the stored hierarchy.
     *
     * In v01, counts may either be maintained incrementally during parsing or
     * recomputed in a finalization pass. This method exists so the final policy
     * remains explicit and easy to revise.
     */
    void finalizeCounts();

    /**
     * @brief Read current authoritative orientation/frame state.
     */
    [[nodiscard]] const OrientationState& orientationState() const;

    /**
     * @brief Overwrite orientation/frame state after a successful workflow step.
     */
    void setOrientationState(const OrientationState& state);

private:
    // -------------------------------------------------------------------------
    // Source metadata
    // -------------------------------------------------------------------------

    std::string source_path_;

    // -------------------------------------------------------------------------
    // Structural hierarchy
    // -------------------------------------------------------------------------

    std::vector<Chain> chains_;
    OrientationState orientation_state_{};

    // -------------------------------------------------------------------------
    // Summary counters
    // -------------------------------------------------------------------------

    std::size_t total_atoms_ = 0;
    std::size_t total_residues_ = 0;
    std::size_t total_subunits_ = 0;
    std::size_t total_hetatm_ = 0;
    std::size_t total_altloc_ = 0;
    std::size_t total_skipped_records_ = 0;
};

// NOTE ON PARSER MUTATION ACCESS:
//
// The parser needs a clean way to continue building the currently active
// subunit after it has already been inserted into the Capsid. Providing
// lastChain() and addAtomToLastChainResidue() is a small but important design
// improvement over relying on const_cast or exposing overly broad mutable
// access patterns.
//
// For v01, this keeps the single-pass parser simple while preserving clearer
// ownership and const-correctness boundaries.

// NOTE ON COUNT OWNERSHIP:
//
// In v01, Capsid stores assembly-level counters directly because summary
// reporting is one of the core deliverables and these values are needed often.
// This keeps reporting simple and avoids scattering top-level statistics across
// multiple modules.
//
// A future revision could centralize some of these values in a dedicated report
// or statistics layer if the architecture grows more complex. For the
// foundation release, keeping the counters in Capsid is the most direct and
// maintainable choice.

// NOTE ON FINALIZATION POLICY:
//
// The design document allows either incremental counting during parsing or a
// dedicated finalization pass that recomputes totals from the stored hierarchy.
// The finalizeCounts() method is included so that policy remains explicit.
//
// For v01, this is useful because it keeps the implementation flexible while
// also giving us a clean place to validate counter consistency if needed later.

#endif // CAPDAT_CAPSID_HPP
