# CapDAT Technical Specification and Development Guide

## Purpose

This document consolidates the current CapDAT technical specification and a practical development guide into a single repository-ready reference. It is intended to support onboarding, architectural understanding, and the implementation of new functionality in a way that remains consistent with the current codebase.

The document is based on the present `main` branch implementation of **CapDAT v0.1.0**, including the build system, command-line interface, parser, core domain model, geometry and symmetry utilities, reorientation workflow, export subsystem, structural summary module, and current tests.

---

## 1. Project Overview

**CapDAT (Capsid Data Analysis Toolkit)** is a C++17 command-line application for the structural analysis of viral capsids represented in fixed-column PDB-like coordinate files. The current release is intentionally a foundation release. Its main goal is to establish a trustworthy and extensible structural core rather than a large analytical feature set.

# CapDAT Geometry Branch — Updated Technical Specification Report

## Scope and basis

This report is based on the **`geometry` branch** of `tripplab/CapDAT`, not on `main`. It reflects the current code-visible state of the branch as exposed by the repository sources, especially:

- `CMakeLists.txt`  
  <https://github.com/tripplab/CapDAT/blob/geometry/CMakeLists.txt>
- `src/main.cpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/src/main.cpp>
- `include/geometry_analysis.hpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_analysis.hpp>
- `src/geometry_analysis.cpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/src/geometry_analysis.cpp>
- `include/geometry_symmetry.hpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_symmetry.hpp>
- `include/reorientation_workflow.hpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/include/reorientation_workflow.hpp>
- `include/capsid.hpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/include/capsid.hpp>
- `include/export_capsid.hpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/include/export_capsid.hpp>

This is not a generic project summary. It is a code-driven technical specification intended to support **new functionality design and implementation** on top of the actual geometry branch architecture.

---

## 1. Executive technical state

The `geometry` branch has evolved CapDAT from a parser/export/reorientation foundation into a **multi-stage local geometry-analysis pipeline** centered on **fold-aligned cylindrical patch analysis**.

At code level, the branch now exposes a complete geometry workflow through `runFoldPatchGeometryAnalysis(...)` in `geometry_analysis`, with the following realized stages:

1. **Stage 1 — fold resolution and in-place reorientation to +Z**
2. **Stage 2 — cylindrical patch selection**
3. **Stage 3 — analytical patch normalization with vdW assignment**
4. **Stage 4 — raw outer/inner sheet detection by vertical line–sphere contacts**
5. **Stage 5 — surface preparation and domain masks**
6. **Stage 6 — smooth field reconstruction over the patch grid**
7. **Stage 7 — optional smoothing / regularization**, with two method families:
   - `smooth`
   - `thin_plate_grid_fit`

This means the branch is no longer just “preparing geometry”; it already contains a real numerical surface pipeline, mesh export logic, provenance-aware traceability, and runtime artifact generation.

---

## 2. Build and repository state

### 2.1 Build system

The branch remains a **CMake-based C++17** project. The executable is still `capsid_analyzer`. Tests are enabled through `CAPDAT_BUILD_TESTS`, and warnings through `CAPDAT_ENABLE_WARNINGS`. The build graph explicitly includes `src/geometry_analysis.cpp`, which confirms geometry is now a first-class compiled subsystem, not an experimental side file.

Source: `CMakeLists.txt`  
<https://github.com/tripplab/CapDAT/blob/geometry/CMakeLists.txt>

### 2.2 Primary compiled modules

The branch’s main executable is built from these modules:

- `main.cpp`
- core structure/domain:
  - `atom.cpp`
  - `residue.cpp`
  - `chain.cpp`
  - `capsid.cpp`
- parser:
  - `pdb_parser.cpp`
- geometry foundation:
  - `geometry_symmetry.cpp`
- workflow:
  - `reorientation_workflow.cpp`
- export:
  - `export_capsid.cpp`
- reporting:
  - `structural_summary.cpp`
  - `summary_reporter.cpp`
- infrastructure:
  - `logger.cpp`
  - `timer.cpp`
- geometry pipeline:
  - `geometry_analysis.cpp`

This matters because new functionality should follow the existing separation instead of bloating `main.cpp` or piggybacking on parser internals.

### 2.3 Test coverage visible in build

The branch defines dedicated test executables for:

- structural summary
- geometry symmetry
- reorientation workflow
- geometry analysis
- export

This is a strong signal that the intended development style is **modular feature-specific tests**, not one monolithic integration-only test harness.

Source: `CMakeLists.txt`  
<https://github.com/tripplab/CapDAT/blob/geometry/CMakeLists.txt>

---

## 3. High-level architecture

The current architecture is best understood as seven layers.

### 3.1 CLI/orchestration layer

`src/main.cpp` parses flags, enforces top-level invariants, builds config structs, executes the parser, prints structural summary, invokes reorientation if requested, invokes geometry analysis if enabled, and optionally exports the final current coordinates.

Source: `src/main.cpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/src/main.cpp>

### 3.2 Parsing layer

`PdbParser` remains the authoritative fixed-column input interpreter. Geometry does **not** parse raw file text again. That separation is still intact.

### 3.3 Persistent structural domain

The persistent object graph remains:

- `Atom`
- `Residue`
- `Chain`
- `Capsid`

The important branch-specific point is that `Capsid` now carries explicit **orientation/frame state**, so downstream geometry can tell whether coordinates are still in original parsed frame or already transformed.

Source: `include/capsid.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/capsid.hpp>

### 3.4 Canonical geometry/symmetry layer

`geometry_symmetry` provides the reusable icosahedral reference math:

- canonical folds
- vector/matrix primitives
- normalization, dot, cross, angle utilities
- rotation construction and application
- IAU definition and classification helpers
- proper-rotation validation

This module is intentionally parser-independent and is the single correct place for reusable fold/axis math.

Source: `include/geometry_symmetry.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_symmetry.hpp>

### 3.5 Workflow mutation layer

`reorientation_workflow` resolves user-facing rotation requests and applies in-place coordinate changes to the `Capsid`. Geometry Stage 1 reuses this workflow instead of duplicating rotation logic.

Source: `include/reorientation_workflow.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/reorientation_workflow.hpp>

### 3.6 Export/serialization layer

`export_capsid` writes the **current in-memory accepted structure**, including transformed coordinates when present. It also supports subset export via `atom_subset`, which is essential for patch/contact exports in geometry stages.

Source: `include/export_capsid.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/export_capsid.hpp>

### 3.7 Derived geometry-analysis layer

`geometry_analysis` is now a large orchestration-plus-algorithm module implementing the Stage 1–7 patch pipeline and its artifact generation.

Source: `include/geometry_analysis.hpp`, `src/geometry_analysis.cpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_analysis.hpp>  
<https://github.com/tripplab/CapDAT/blob/geometry/src/geometry_analysis.cpp>

---

## 4. CLI specification as implemented in the geometry branch

The geometry branch CLI materially extends the older interface.

### 4.1 General options

- `--input`
- `--log`
- `--verbose`
- `--quiet`
- `--include-hetatm`
- `--export-final`

### 4.2 Reorientation options

- `--reorient`
- `--align-fold`
- `--align-vector`
- `--align-axis`

### 4.3 Geometry-analysis enablement

- `--geometry-analysis`
- `--debug`

### 4.4 Core geometry parameters

- `--geometry_fold_type`
- `--geometry_fold_index`
- `--geometry_cylinder_radius`
- `--dvdW`
- `--geometry_grid_spacing`
- `--geometry_min_atoms_in_patch`
- `--geometry_boundary_margin`
- `--geometry_support_radius`
- `--geometry_min_support_nodes`
- `--geometry_reliable_radius`
- `--geometry_out_prefix`
- `--export_mesh_format`
- `--split_in_out_mesh`
- `--surf_min_separation`

### 4.5 Stage 7 smoothing / regularization options

- `--geometry_smooth_enabled`
- `--geometry_smooth_weight`
- `--geometry_smooth_max_iterations`
- `--geometry_smooth_convergence_tolerance`
- `--geometry_smooth_pin_seed`
- `--geometry_smooth_enforce_non_crossing`
- `--geometry_smooth_min_separation`
- `--geometry_smooth_export_meshes`
- `--geometry_s7_method`
- `--geometry_s7_lambda`
- `--geometry_s7_data_weight_seed`
- `--geometry_s7_data_weight_interp`
- `--geometry_s7_use_reliable_core_only_for_fit`
- `--geometry_s7_boundary_condition_mode`
- `--geometry_s7_solver_max_iterations`
- `--geometry_s7_solver_tolerance`
- `--geometry_s7s6_deltas`

Source: `src/main.cpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/src/main.cpp>

### 4.6 CLI invariants enforced in code

`main.cpp` explicitly rejects:

- missing `--input`
- `--align-fold` and `--align-vector` used together
- `--reorient` without exactly one source
- source flags without `--reorient`
- simultaneous `--geometry-analysis` and `--reorient`
- invalid Stage 7 lambda / weights / solver params

This is not cosmetic. It tells you where the branch draws the line between:
- **top-level mode consistency checks** in `main.cpp`
- **feature-specific deep validation** inside the geometry subsystem

---

## 5. Core domain model relevant to geometry development

### 5.1 `Capsid` is the authoritative coordinate owner

Geometry stages traverse the accepted atoms already stored in the `Capsid`. They do not create a second full structural object. Stage 1 may rotate coordinates **in place**, and later stages operate on that updated state.

### 5.2 Orientation state is central

`Capsid::OrientationState` records:

- whether coordinates are still in original parsed frame
- whether reorientation was applied in place
- whether the transform was identity
- rotation matrix
- rotation axis and angle
- source mode / description / direction
- requested target axis
- target direction

Source: `include/capsid.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/capsid.hpp>

This is the branch’s key architectural defense against frame ambiguity. Any new geometry feature that assumes canonical fold alignment must consult this state or be explicit about when it runs.

### 5.3 Patch-level derived structures are separate from persistent domain

`geometry_analysis.hpp` introduces analysis-layer structs rather than polluting `Atom` or `Capsid`:

- `CylinderMembership`
- `PatchAtom`
- `AnalyticalPatch`
- stage-specific result structs
- top-level `GeometryAnalysisResult`

This is correct architecture. Persistent biology stays in the domain model; analytical state lives in derived result objects.

Source: `include/geometry_analysis.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_analysis.hpp>

---

## 6. `geometry_symmetry` as the canonical math source

### 6.1 Public types

- `Vector3`
- `Matrix3`
- `FoldDefinition`
- `RotationDefinition`
- `IauDefinition`

### 6.2 Public capabilities

- canonical fold registry
- fold lookup by name and by `(type,index)`
- fold reference and unit vectors
- norm / normalize / dot / cross
- angle and similarity helpers
- nearest-fold lookup
- Euclidean and point-to-axis distance
- rotation alignment between directions
- align fold or arbitrary direction to +Z
- apply rotation to direction / point / point set
- IAU margins and classification
- proper rotation matrix validation

Source: `include/geometry_symmetry.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_symmetry.hpp>

### 6.3 Canonical folds actually exposed

The module provides these canonical names:

- `2_0`
- `2_1`
- `3_0`
- `3_1`
- `5_0`

This is the authoritative branch vocabulary for fold resolution. New functionality must reuse these definitions instead of hardcoding its own fold map.

---

## 7. Reorientation workflow contract

### 7.1 Request model

`ReorientationRequest` carries:

- whether reorientation is requested
- source mode (`fold` or `custom_vector`)
- fold name or vector text
- target axis
- export request metadata
- verbosity

### 7.2 Result model

`ReorientationResult` returns:

- status
- resolved source description
- source direction
- target direction
- rotation matrix
- rotation axis
- rotation angle
- whether coordinates were modified
- export info
- messages

Source: `include/reorientation_workflow.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/reorientation_workflow.hpp>

### 7.3 Geometry Stage 1 depends on it

Stage 1 in `geometry_analysis.cpp` calls `applyReorientationWorkflow(...)` with a fold-derived request targeting `+Z`. That is the branch’s intended reuse pattern: **workflow layer owns transformation semantics; analysis layer consumes it**.

---

## 8. Geometry-analysis public API

The geometry branch exposes the following core geometry functions:

- `classifyPatchCylinder(...)`
- `normalizeElementSymbol(...)`
- `vdwRadius(...)`
- `makePatchAtom(...)`
- `intersectVerticalLineWithSphere(...)`
- `detectRawFirstContactAtNode(...)`
- `buildStage4RegularGrid(...)`
- `prepareGeometryAnalysisStage1(...)`
- `runGeometryAnalysisStage2PatchSelection(...)`
- `runGeometryAnalysisStage3PatchNormalization(...)`
- `runGeometryAnalysisStage4RawSheetDetection(...)`
- `runGeometryAnalysisStage5SurfacePreparation(...)`
- `runGeometryAnalysisStage6SurfaceReconstruction(...)`
- `runGeometryAnalysisStage7SurfaceSmoothing(...)`
- `runFoldPatchGeometryAnalysis(...)`

Source: `include/geometry_analysis.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_analysis.hpp>

This is the real feature surface of the branch. Any new geometry functionality should either:
- extend one of these stage families,
- or add a clean new stage/module instead of hiding logic in `main.cpp`.

---

## 9. Stage-by-stage technical specification

## 9.1 Stage 1 — preparation and canonical fold alignment

### Purpose

Resolve a requested canonical fold, reorient the in-memory `Capsid` so that fold aligns with `+Z`, and verify the working frame.

### Inputs

- mutable `Capsid&`
- `FoldPatchAnalysisConfig`
- `ParserConfig`
- `Logger*`

### Behavior

- validates fold type and radius-level config
- resolves fold by `(fold_type, fold_index)` through `geometry_symmetry::foldByTypeIndex(...)`
- constructs a fold-based `ReorientationRequest`
- applies reorientation in place
- checks `Capsid::orientationState()` after transform
- optionally exports a rotated full capsid PDB when `export_rotated_capsid` is enabled

### Outputs

`GeometryPreparationResult` includes:

- requested and resolved fold info
- reference and unit vectors
- identity/non-identity transform state
- rotation matrix / axis / angle
- whether coordinates changed
- final frame description
- optional rotated capsid export path
- nested `ReorientationResult`
- messages

### Key design consequence

Stage 1 mutates coordinates **once** and then downstream stages assume the analysis frame is already established. New geometry stages should not secretly rotate again.

Source: `src/geometry_analysis.cpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/src/geometry_analysis.cpp>

---

## 9.2 Stage 2 — cylindrical patch selection

### Purpose

Extract the fold-centered local patch as the subset of atoms satisfying the cylindrical selection rule in the Stage 1 working frame.

### Selection rule

For each accepted atom:

- compute `radial_xy = sqrt(x^2 + y^2)`
- require `z > 0`
- require `radial_xy <= cylinder_radius`

The selection is encoded through `CylinderMembership`.

### Outputs

`GeometryPatchSelectionResult` stores:

- counts of examined and selected atoms
- rejection counters
- radius and min-atom threshold
- exported patch PDB path
- `patch_atoms`
- `selected_atom_refs`
- messages

### Export semantics

Stage 2 always exports patch atoms as a PDB subset using `ExportCapsidWriter` with `atom_subset`.

### Architectural note

The stage preserves stable pointers back to original atoms, which is critical for traceability and later exports.

---

## 9.3 Stage 3 — analytical patch normalization

### Purpose

Convert Stage 2 selected atoms into a normalized analytical patch representation suitable for geometric calculations.

### Behavior

For each Stage 2 selected atom:

- re-derive a normalized `PatchAtom`
- normalize element symbol
- assign vdW radius
- apply optional `delta_vdw`
- preserve selection-membership facts
- preserve pointer to original atom

### vdW resolution policy

Resolution is hierarchical:

1. explicit normalized element if recognized
2. inferred element from atom name if possible
3. fallback radius `1.70 Å`

Lookup table currently covers:

- H
- C
- N
- O
- S
- P
- SE

### Outputs

`GeometryPatchNormalizationResult` contains an `AnalyticalPatch` with:

- normalized atoms
- original atom refs
- cylinder radius
- atom count
- counts of explicit / inferred / fallback vdW assignments
- export path inherited from Stage 2

### Important meaning

This stage is not a vague “prep” step. It is the **normalization boundary** between raw selected structure and analytically usable patch data.

---

## 9.4 Stage 4 — raw outer/inner sheet detection

### Purpose

Construct raw outer and inner envelopes over a regular XY disk grid by vertical line–sphere intersection against patch atoms.

### Grid construction

`buildStage4RegularGrid(...)` creates a regular grid spanning:

- `x ∈ [-R, +R]`
- `y ∈ [-R, +R]`
- spacing = `grid_spacing`

### Node-level contact rule

For each inside-disk grid node `(x,y)`:

- intersect the vertical line through `(x,y)` with every patch-atom sphere
- for all intersected spheres:
  - outer raw contact = maximum `z_high`
  - inner raw contact = minimum `z_low`
- node is valid if both are finite and `z_inner <= z_outer + tol`

This is implemented by:

- `intersectVerticalLineWithSphere(...)`
- `detectRawFirstContactAtNode(...)`

### Outputs

`GeometryStage4RawSheetResult` includes:

- grid definition
- raw outer/inner scalar fields
- inside-disk and valid masks
- candidate counts
- serial-number provenance per node
- patch-atom-index provenance per node
- raw contact records
- aggregate counts:
  - valid/invalid nodes
  - zero-thickness nodes
  - negative-thickness nodes
  - unique contact atoms
  - unique contact serials
- optional CSV artifact paths
- contact-atom PDB path
- timestamps and runtime
- messages

### Debug artifacts

When `debug` is enabled, Stage 4 can emit:

- Stage 3 normalized atoms CSV
- raw outer CSV
- raw inner CSV
- valid-mask CSV
- outer-only mask CSV
- inner-only mask CSV
- negative-thickness mask CSV
- Stage 4 summary CSV

It always exports a PDB containing the union of contact atoms used in the raw envelope.

### Development significance

Stage 4 is the first real geometric measurement stage. It converts the atom cloud into gridded envelope evidence, with strong provenance tracking.

---

## 9.5 Stage 5 — surface preparation and domain definition

### Purpose

Transform Stage 4 raw envelope output into the masks and seed domains needed for Stage 6 reconstruction and Stage 7 metric evaluation.

### Derived parameters

If not explicitly set, Stage 5 derives:

- `boundary_margin = 2 * grid_spacing`
- `support_radius = 2.5 * grid_spacing`
- `reliable_radius = cylinder_radius - boundary_margin`

### Main constructs

Stage 5 builds:

- `z_outer_seed`
- `z_inner_seed`
- `outer_seed_mask`
- `inner_seed_mask`
- `paired_seed_mask`
- `boundary_exclusion_mask`
- `interp_allowed_outer_mask`
- `interp_allowed_inner_mask`
- `paired_interp_allowed_mask`
- `hard_invalid_mask`
- `reliable_core_mask`

### Logic

- every Stage 4 valid node becomes a paired seed
- boundary-excluded nodes are marked near the edge
- non-seed nodes are checked for enough nearby support within `support_radius`
- nodes without enough support become hard-invalid
- reliable-core nodes are those within `reliable_radius` and not invalid

### Provenance

Stage 5 also aggregates:

- unique outer seed atom serials
- unique inner seed atom serials
- unique seed patch-atom indices
- union counts

### Outputs

`GeometryStage5SurfacePrepResult` stores all masks, seed fields, derived radii, provenance counts, optional artifact paths, and messages.

### Technical importance

This stage defines the branch’s **trust geometry**:
- edge exclusion
- interpolation domain
- hard-invalid holes
- reliable metric core

That is the branch’s real quality-control boundary for downstream metrics.

---

## 9.6 Stage 6 — surface reconstruction

### Purpose

Reconstruct continuous outer and inner scalar fields over the patch grid from Stage 5 seeds and interpolation-allowed nodes.

### Numerical method

Stage 6 implements a deterministic grid relaxation scheme:

1. keep paired seeds as anchors
2. initialize interpolation-allowed nodes from finite 4-neighbor averages when available
3. iterate Jacobi-like relaxation with smoothing factor
4. optionally enforce non-crossing between outer and inner surfaces

The core internal routine is `runStage6FieldReconstruction(...)`.

### Neighborhood model

Stage 6 uses a **4-neighbor stencil**, not diagonal neighbors.

### Non-crossing policy

If enabled and `z_outer < z_inner + min_separation`, the stage symmetrically adjusts both values around the midpoint.

### Outputs

`GeometryStage6SurfaceReconstructionResult` includes:

- reconstructed outer/inner fields
- masks:
  - reconstructed
  - final_valid_analysis
  - non_crossing_adjustment
  - obj_vertex
- counts:
  - reconstructed nodes
  - seed nodes
  - interpolation nodes
  - unresolved nodes
- iteration stats
- separation stats
- optional CSV paths
- mesh export paths and counts
- messages

### Mesh export

Stage 6 can export:

- STL or OBJ
- combined inner+outer mesh file or split files
- using `metric_domain_mask`-like scalar node availability logic for face generation

The implementation includes explicit OBJ and ASCII STL writers.

### Practical meaning

Stage 6 is the first stage that yields a surface-like reconstruction suitable for later smoothing or direct visualization.

---

## 9.7 Stage 7 — smoothing / regularization

### Purpose

Produce the final smoothed or regularized patch surfaces used for reliable metric-domain analysis.

### Method selector

`FoldPatchAnalysisConfig::Stage7Method` supports:

- `smooth`
- `thin_plate_grid_fit`

### Method A: `smooth`

A legacy local smoothing pass over Stage 6 reconstructed fields:

- iterates over reconstructed nodes
- optionally preserves seed values
- uses 4-neighbor averaging relaxation
- supports convergence and iteration control

Implemented in `runStage7FieldSmoothingLegacy(...)`.

### Method B: `thin_plate_grid_fit`

A discrete regularized fit over the grid:

- defines fit domain from reconstructed nodes, optionally restricted to reliable core
- builds node-specific fidelity weights
- constructs boundary policy:
  - `free`
  - `fixed_to_stage6`
  - `soft_to_stage6`
- solves a regularized system using a conjugate-gradient-like iterative approach
- computes residuals and a discrete bending-energy proxy

Implemented through:
- `buildStage7FitDomainMask(...)`
- `buildStage7FidelityWeights(...)`
- `buildStage7BoundaryMask(...)`
- `applyStage7Laplacian(...)`
- `applyStage7System(...)`
- `runStage7FieldThinPlateGridFit(...)`

### Post-fit smoothing constraints

After either method:

- optional non-crossing enforcement is applied
- metric domain is restricted to reliable core
- smooth separation statistics are computed

### Outputs

`GeometryStage7SmoothedSurfaceResult` contains:

- smoothed outer/inner fields
- masks:
  - reconstructed
  - reliable_core
  - smooth_valid
  - metric_domain
  - smooth_non_crossing_adjustment
- iteration / solver stats
- fit residuals and bending energies
- separation stats
- Stage 7 vs Stage 6 deltas if enabled
- optional CSV paths
- optional mesh paths and counts
- messages

### Artifact generation

Depending on config, Stage 7 can emit:

- outer smooth CSV
- inner smooth CSV
- smooth-valid mask CSV
- metric-domain mask CSV
- non-crossing-adjustment mask CSV
- summary CSV
- Stage 7 vs Stage 6 outer/inner/thickness delta CSVs
- smoothed meshes in OBJ or STL

### What Stage 7 is really doing

This stage is the branch’s current approximation to a **final analysis-ready smooth local shell representation**. It is not yet curvature/thickness metrics proper, but it already defines the final metric domain and a numerical regularization layer.

---

## 10. Geometry configuration model

The geometry branch centralizes geometry-analysis runtime behavior in `FoldPatchAnalysisConfig`.

Source: `include/geometry_analysis.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_analysis.hpp>

### 10.1 Main configuration groups

#### Global geometry controls

- enabled
- debug
- fold_type
- fold_index
- cylinder_radius
- delta_vdw
- grid_spacing
- min_atoms_in_patch
- output_prefix
- export_rotated_capsid

#### Stage 5 controls

- stage5_boundary_margin
- stage5_support_radius
- stage5_min_support_nodes
- stage5_reliable_radius

#### Stage 6 controls

- stage6_smoothing_weight
- stage6_max_iterations
- stage6_convergence_tolerance
- stage6_enforce_non_crossing
- stage6 mesh export controls
- stage6 split-mesh controls
- stage6 minimum separation control used by the implementation path

#### Stage 7 controls

- enabled
- method
- smoothing weight
- iteration limits
- convergence tolerance
- preserve seed values
- thin-plate lambda
- fidelity weights
- reliable-core-only fit switch
- boundary condition mode
- solver iteration/tolerance
- Stage 7 vs Stage 6 delta export
- non-crossing
- minimum separation
- mesh export

### 10.2 Design implication

Any new geometry feature that changes algorithmic behavior should be added here first rather than hidden in hardcoded constants.

---

## 11. Artifact/export model generated by geometry analysis

The geometry branch is heavily artifact-oriented. That is one of its most important practical characteristics.

### 11.1 PDB exports

The pipeline can emit:

- rotated full capsid PDB from Stage 1
- Stage 2 patch atom subset PDB
- Stage 4 raw contact atom subset PDB
- final full exported current capsid via `--export-final`

### 11.2 CSV exports

The pipeline can emit stage-specific CSVs for:

- normalized patch atoms
- raw outer/inner sheets
- valid and anomaly masks
- Stage 5 seed and domain masks
- Stage 6 reconstructed fields and masks
- Stage 7 smoothed fields and masks
- Stage 7 vs Stage 6 deltas
- stage summary tables

### 11.3 Mesh exports

The pipeline can emit:

- Stage 6 reconstructed surfaces
- Stage 7 smoothed or thin-plate surfaces

in:

- `obj`
- `stl`

and either:

- combined inner+outer file
- split inner and outer files

### 11.4 Run-summary export

The pipeline writes a top-level run-summary JSON:
- `<output_prefix>_run_summary.json`

This captures geometry and parser config state and is useful for reproducibility.

Source: `src/geometry_analysis.cpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/src/geometry_analysis.cpp>

---

## 12. Data provenance model

One of the strongest parts of this branch is that provenance is not hand-wavy.

### 12.1 Atom traceability

`PatchAtom` stores:
- normalized working coordinates
- element/vdW data
- cylinder-membership facts
- pointer to original atom

### 12.2 Node traceability

Stage 4 and Stage 5 carry:
- contact serial numbers
- contact patch-atom indices
- unique provenance sets

### 12.3 Export traceability

Subset PDB exports are driven from original atom pointers, not lossy reconstructed pseudo-records.

This is good engineering. It reduces the risk of analysis outputs drifting away from biological input identity.

---

## 13. Numerical design characteristics visible in code

This matters if you are about to extend the branch.

### 13.1 Grid-based analysis, not mesh-first analysis

The current geometry pipeline is fundamentally **regular-grid scalar-field based** over XY, not triangulation-native from atoms.

### 13.2 Locality is defined in cylindrical coordinates around +Z

Everything downstream assumes Stage 1 established the working frame and that the selected fold axis is now `+Z`.

### 13.3 Outer/inner shell are modeled as separate scalar surfaces

The shell is not yet handled as a volumetric signed-distance field or explicit manifold object. It is represented as two scalar height maps:
- outer: `z_outer(...)`
- inner: `z_inner(...)`

### 13.4 Thickness is currently vertical-gap based in Stage 7 deltas/metadata

The metadata explicitly labels local thickness as `"vertical_z_difference"`, not normal-distance thickness.

That is a major conceptual limitation if future work wants physically meaningful thickness.

### 13.5 Stage 6 and legacy Stage 7 both use local 4-neighbor relaxation

This means the present numerical core is simple, deterministic, and easy to reason about, but also limited in expressiveness and isotropy.

### 13.6 Thin-plate Stage 7 is discrete-grid regularization, not a full continuous TPS implementation

The method is closer to a **discrete Laplacian/Laplacian-squared regularized grid fit** than a classical radial-basis-function thin-plate spline formulation over arbitrary points.

That distinction matters. The naming is directionally fair, but mathematically it is a grid-regularized discrete solver.

---

## 14. What is already stable enough to build on

The following parts look structurally solid enough to use as foundations for new development:

### 14.1 Canonical fold and rotation layer

`geometry_symmetry` and reorientation workflow are clearly separated and reusable.

### 14.2 Patch selection and normalization boundary

Stage 2 plus Stage 3 form a clean and useful upstream interface for any new local analysis.

### 14.3 Provenance-aware subset export

The writer and pointer-based subset mechanism are already good enough to support many new geometry/debug outputs.

### 14.4 Stage-structured result objects

Each stage returns a plain result struct with explicit fields and messages. This is exactly the right pattern for extending the branch.

### 14.5 Artifact-first debugging model

The branch already assumes complex stages need inspectable outputs. That is the correct mindset for numerical geometry development.

---

## 15. Current limitations and blind spots

Here is the part people usually avoid saying clearly.

### 15.1 Geometry analysis is still monolithic in one `.cpp`

`src/geometry_analysis.cpp` is now very large and mixes:

- config validation
- algorithmic kernels
- CSV writers
- mesh writers
- provenance bookkeeping
- stage orchestration

That is manageable for now, but it will become a drag fast. New functionality should not continue piling into the same file forever.

### 15.2 No final metric layer yet

Despite Stage 7’s sophistication, the branch still does **not** expose final scientific metrics such as:

- thickness by local normals
- mean curvature
- Gaussian curvature
- area metrics
- robust QC summary metrics over the reliable core

The pipeline has prepared the surfaces, but the final analysis stage is still missing.

### 15.3 Thickness is not yet physically rigorous

The code’s own metadata says thickness is a vertical `z` difference. That is only acceptable if the local patch is sufficiently aligned and gently curved. It is not a general shell-thickness definition.

### 15.4 Surface quality control is heuristic

Boundary margin, support radius, hard-invalid logic, and reliable-core restriction are sensible, but still heuristic. They are not yet grounded in a formal error model.

### 15.5 Mesh export is topology-naive

The mesh writers generate faces from valid grid quads. They do not do advanced topology repair, hole handling, normal consistency enforcement, or manifold validation.

### 15.6 Parser/domain limitations remain upstream

The geometry branch still inherits upstream limitations such as:
- heuristic protein-only residue acceptance
- simplistic internal subunit reconstruction assumptions
- dependence on accepted parsed coordinates already being meaningful

### 15.7 Metrics and smoothing live in the same conceptual subsystem

The branch currently conflates “surface generation”, “regularization”, and “analysis-preparation” inside one geometry module. That is functional, but not yet ideal architecture for long-term growth.

---

## 16. Extension points for new functionality

These are the most natural places to add new work.

### 16.1 Post-Stage 7 metric computation

The cleanest near-term addition is a new stage or module consuming:
- `GeometryStage7SmoothedSurfaceResult`
- `GeometryStage5SurfacePrepResult`
- maybe `GeometryStage6SurfaceReconstructionResult`

to compute:
- thickness
- curvature
- patch area
- reliability/QC summaries

### 16.2 Alternate Stage 7 methods

The branch already has a method enum. That means adding another regularizer is straightforward **if** it obeys the same result contract.

### 16.3 Better Stage 6/7 solvers

You could replace or augment the current relaxation schemes with:
- weighted Laplacian
- biharmonic smoother
- sparse linear solve backend
- anisotropic smoothing
- masked finite-difference operators

### 16.4 Better metric-domain definition

The current reliable-core logic is radial and boundary-based. Future work could define:
- confidence scores
- local coverage density
- uncertainty propagation
- gradient-based artifact detection

### 16.5 Separate surface/mesh utilities module

This branch is overdue for extraction of:
- CSV writing
- grid utilities
- mesh writers
- solver helpers

into smaller internal modules.

---

## 17. Recommended refactor boundaries before major new features

If development is about to accelerate, the following refactor is worth doing early.

### 17.1 Split `geometry_analysis.cpp`

Break it into at least:

- `geometry_patch_selection.cpp`
- `geometry_patch_normalization.cpp`
- `geometry_raw_sheet.cpp`
- `geometry_surface_prep.cpp`
- `geometry_surface_reconstruction.cpp`
- `geometry_surface_smoothing.cpp`
- `geometry_mesh_export.cpp`
- `geometry_csv_export.cpp`

You do not need to expose all of them publicly, but keeping everything in one implementation file is asking for avoidable coupling.

### 17.2 Add a dedicated metric layer

Create something like:

- `include/geometry_metrics.hpp`
- `src/geometry_metrics.cpp`

and keep it read-only over Stage 7 outputs.

### 17.3 Separate algorithm kernels from artifact writers

Right now too much stage code is mixed with output-writing code. That makes algorithm changes harder to test cleanly.

---

## 18. Practical engineering conclusions

### 18.1 What this branch is, technically

The `geometry` branch is currently a **fold-aligned local shell surface reconstruction pipeline** built on top of the CapDAT structural core.

### 18.2 What it is not yet

It is **not yet** a finished geometry-analysis product. It lacks final scientific metrics, more rigorous shell thickness definitions, stronger numerical abstraction, and cleaner module boundaries.

### 18.3 What it is ready for

It is already strong enough to support development of:

- patch thickness metrics
- curvature metrics
- QC metrics
- alternate smoothers
- improved surface reconstruction
- mesh/report export enhancements

### 18.4 The single most important development rule

Do **not** bypass the existing branch architecture.

That means:

- do not duplicate fold math outside `geometry_symmetry`
- do not parse raw coordinate text inside geometry code
- do not create hidden frame assumptions
- do not store analytical scratch state inside persistent domain classes unless it is truly authoritative
- do not keep enlarging `main.cpp`
- do not collapse new numerical work into one giant undifferentiated stage

---

## 19. Suggested mental model for future developers

The most accurate mental model for the current geometry branch is:

- `Capsid` = authoritative accepted structure and current coordinates
- `geometry_symmetry` = canonical icosahedral math and fold vocabulary
- `reorientation_workflow` = explicit coordinate-frame mutation
- `geometry_analysis` = staged local patch pipeline over the current capsid
- `export_capsid` = authoritative serializer of current accepted state
- stage result structs = the real analytical contract

That is the backbone you should preserve when adding new functionality.


# CapDAT Geometry Branch — Synthesized Development Guide for Implementing New Functionality

## 1. First: classify the feature correctly

Before touching code, decide what the feature actually is. Most bad changes come from misclassification.

### Put it in the parser only if

- you need new fields from input records
- you need to change acceptance policy
- you need to change how structure is reconstructed from file records
- the information belongs to the persistent biological structure itself

### Put it in the persistent domain only if

- the data must live as authoritative state across workflows/modules
- multiple later modules need it as part of the structural truth
- it is not merely a derived temporary analysis artifact

### Put it in `geometry_symmetry` if

- it is reusable icosahedral math
- it involves fold lookup, directions, angles, axis distances, canonical reference vectors, IAU logic, or generic rotations
- more than one future module could need it

### Put it in `reorientation_workflow` if

- it is a user-triggered coordinate transform
- it mutates the current frame
- it needs request/result semantics
- it should update `Capsid::OrientationState`

### Put it in `geometry_analysis` or a sibling geometry module if

- it operates on accepted parsed coordinates
- it is patch/surface/metric logic
- it is derived analysis, not raw parsing
- it should consume stage results or `const Capsid&`

### Put it in export code only if

- it is purely serialization
- it turns current state or derived geometry results into files
- it should not contain the scientific algorithm itself

---

## 2. Follow the existing branch workflow, not your own improvised one

The geometry branch already has a functional staged pipeline. Reuse it.

### Current geometry workflow contract

1. parse input into `Capsid`
2. compute general structural summary
3. optionally apply standalone reorientation workflow
4. or run geometry-analysis pipeline:
   - Stage 1: align requested fold to +Z
   - Stage 2: select cylindrical patch
   - Stage 3: normalize analytical patch
   - Stage 4: detect raw outer/inner contacts
   - Stage 5: prepare seeds and masks
   - Stage 6: reconstruct surfaces
   - Stage 7: smooth / regularize surfaces
5. optionally export current full capsid
6. write geometry artifacts

Source files to study before implementing anything:

- `src/main.cpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/src/main.cpp>
- `include/geometry_analysis.hpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_analysis.hpp>
- `src/geometry_analysis.cpp`  
  <https://github.com/tripplab/CapDAT/blob/geometry/src/geometry_analysis.cpp>

### Non-negotiable rule

If your new feature depends on local fold-centered geometry, do **not** reinvent preprocessing. Start from the existing stages unless you have a very strong reason not to.

---

## 3. The safest implementation entry points

## 3.1 If you want to add final geometry metrics

This is the cleanest next step.

### Best insertion point

Consume:

- `GeometryStage7SmoothedSurfaceResult`
- `GeometryStage5SurfacePrepResult`
- optionally `GeometryStage6SurfaceReconstructionResult`

### Why here

Because by Stage 7 you already have:

- final smoothed outer/inner fields
- reliable-core mask
- metric-domain mask
- non-crossing already handled
- optional deltas against Stage 6
- mesh-ready scalar surfaces

### Recommended new module

- `include/geometry_metrics.hpp`
- `src/geometry_metrics.cpp`
- `tests/geometry_metrics_tests.cpp`

### Suggested public API shape

Use plain structs like the existing branch style:

- `GeometryMetricsConfig`
- `GeometryMetricsResult`
- `runGeometryMetrics(...)`

Do not dump metric outputs into Stage 7 result unless they are inseparable from Stage 7 itself. Usually they are not.

---

## 3.2 If you want to add a new Stage 7 method

The branch is already designed for this.

### Existing extension hook

`FoldPatchAnalysisConfig::Stage7Method`

Current methods:

- `smooth`
- `thin_plate_grid_fit`

### Implementation pattern

1. extend the enum
2. add CLI option parsing in `main.cpp`
3. add method label handling in `geometry_analysis.cpp`
4. implement a new internal solver function
5. plug it into `runGeometryAnalysisStage7SurfaceSmoothing(...)`
6. populate the same result struct fields so the rest of the pipeline stays stable
7. add targeted tests

### Requirement

Do not return a method-specific ad hoc result shape. The whole point is to preserve Stage 7 as a stable contract.

---

## 3.3 If you want to add a new surface-reconstruction strategy

You have two options.

### Option A: extend Stage 6

Do this if the new method still conceptually belongs to “reconstruction from Stage 5 seeds and interpolation masks.”

### Option B: add a parallel stage/module

Do this if the method changes the conceptual representation too much, for example:

- moving from grid relaxation to mesh-native reconstruction
- moving from scalar sheets to volumetric representation
- introducing confidence/uncertainty modeling as first-class outputs

Rule: if you must distort existing Stage 6 semantics to fit your method, do not force it. Add a sibling module instead.

---

## 4. Respect frame semantics or you will create garbage

This is one of the easiest ways to break scientific meaning.

### What to check

Always know whether your code assumes:

- original parsed frame
- arbitrary current frame
- canonical fold-aligned frame
- Stage 1 working frame

### The authoritative source

`Capsid::orientationState()`

Source: `include/capsid.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/capsid.hpp>

### Rule

If your feature assumes the selected fold is aligned to +Z, make that precondition explicit. Do not silently hope the current coordinates are in the expected frame.

### Practical pattern

For stage-based local geometry code, require successful Stage 1 result or later stage results that already guarantee the working frame.

---

## 5. Reuse the branch’s data contracts instead of inventing new ad hoc containers

The branch already gives you the useful derived-state boundaries.

### Use these when possible

- `PatchAtom`
- `AnalyticalPatch`
- `GeometryPatchNormalizationResult`
- `GeometryStage4RawSheetResult`
- `GeometryStage5SurfacePrepResult`
- `GeometryStage6SurfaceReconstructionResult`
- `GeometryStage7SmoothedSurfaceResult`

### What not to do

Do not create random duplicate containers for:
- patch atoms
- grid nodes
- mask sets
- reconstructed fields
- surface provenance

unless the existing result objects are genuinely insufficient.

---

## 6. Where new code should physically go

## 6.1 Good path for a new read-only analysis feature

Use this shape:

- `include/<feature>.hpp`
- `src/<feature>.cpp`
- `tests/<feature>_tests.cpp`

Example:

- `include/geometry_metrics.hpp`
- `src/geometry_metrics.cpp`
- `tests/geometry_metrics_tests.cpp`

## 6.2 Good path for a geometry helper extraction

If you keep touching the same helper categories repeatedly, split them out from `geometry_analysis.cpp`.

Useful future splits:

- `geometry_grid_utils.cpp`
- `geometry_surface_io.cpp`
- `geometry_mesh_export.cpp`
- `geometry_stage6_solver.cpp`
- `geometry_stage7_solver.cpp`

## 6.3 What not to do

Do not keep adding large unrelated helper blocks to `geometry_analysis.cpp` forever. That file is already too broad.

---

## 7. Coding pattern the branch clearly prefers

The branch has a visible style. Follow it.

### The pattern is

- plain config structs
- plain result structs
- explicit stage entry functions
- explicit validation
- deterministic control flow
- logger messages collected as strings
- artifact paths carried in result structs
- provenance retained wherever possible

### Copy this style

Bad additions are the ones that:
- hide behavior in constructors
- mutate domain objects casually
- use magic globals
- bury state in anonymous side effects
- skip structured result reporting

---

## 8. How to design a new metric module correctly

Here is the pattern worth using.

## 8.1 Inputs

Take the minimum necessary upstream stage outputs.

For example, a metric module might take:

- `const GeometryStage7SmoothedSurfaceResult& stage7`
- `const GeometryStage5SurfacePrepResult& stage5`
- `const FoldPatchAnalysisConfig& config`
- `Logger* logger`

## 8.2 Config

Create a dedicated config only if you actually need new knobs. Do not shove every future metric parameter into `FoldPatchAnalysisConfig` by reflex.

Add to `FoldPatchAnalysisConfig` only when:
- the parameter is core to geometry pipeline orchestration
- the CLI needs to control it directly as part of geometry-analysis mode

Otherwise, separate config is cleaner.

## 8.3 Result

Return a plain result struct containing:

- success flag
- computed scalar metrics
- optional per-node fields
- optional export paths
- messages

## 8.4 Exports

If you export per-node maps, keep them in dedicated writer helpers or a small sibling exporter. Do not mix metric formulas and CSV formatting too aggressively.

---

## 9. Specific guidance for the likely next feature classes

## 9.1 Thickness metrics

### Bad approach

Thickness = `z_outer - z_inner` everywhere forever.

That is easy, but conceptually weak.

### Better approach

Compute both:

- `vertical_thickness`
- and later `normal_projected_thickness`

Then be explicit about which one is reported and why.

### Suggested first deliverable

Implement a metric module that reports:

- mean vertical thickness over metric domain
- median vertical thickness
- min/max
- robust spread
- coverage stats

Then later extend to normal-based thickness.

---

## 9.2 Curvature metrics

### Best starting point

Use Stage 7 smoothed scalar fields over the reliable metric domain.

### Practical route

From `z(x,y)` surfaces, compute finite-difference derivatives:

- first derivatives
- second derivatives
- then:
  - mean curvature
  - Gaussian curvature
  - area element if needed

### Warning

Do not compute curvature on Stage 4 raw fields. That is asking for noise-driven nonsense.

---

## 9.3 Patch area metrics

### Best route

Use Stage 7 metric domain and local surface area element over the scalar grid.

### Important distinction

You may want both:

- projected XY area
- reconstructed surface area

Do not conflate them.

---

## 9.4 Reliability/QC metrics

The branch already has the raw ingredients.

Good QC candidates:

- fraction of inside-disk nodes that survive into metric domain
- fraction of nodes requiring non-crossing adjustment
- Stage 7 vs Stage 6 delta magnitudes
- seed density
- interpolation fraction
- hard-invalid fraction

This is low-hanging fruit and scientifically useful.

---

## 10. Use `geometry_symmetry` aggressively

If you need any of the following, use the existing module:

- fold lookup
- fold direction
- axis-angle comparisons
- point-to-axis distances
- nearest fold
- canonical IAU logic
- rotation construction/application

Source: `include/geometry_symmetry.hpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/include/geometry_symmetry.hpp>

### Rule

If you find yourself hardcoding fold vectors or writing custom direction-alignment math in a new feature, you are doing it wrong.

---

## 11. How to wire a new CLI-controlled feature

Follow the branch’s existing pattern in `main.cpp`.

### Correct steps

1. add local variables near the other geometry options
2. parse the new flag in the argument loop
3. validate top-level invariants if needed
4. map CLI values into config structs
5. keep actual feature logic out of `main.cpp`
6. call a dedicated module function after prerequisites are satisfied

### Wrong pattern

Do not drop real algorithmic logic directly inside `main.cpp`.

Source: `src/main.cpp`  
<https://github.com/tripplab/CapDAT/blob/geometry/src/main.cpp>

---

## 12. How to test new functionality without making a mess

The branch already signals the expected testing strategy.

## 12.1 Preferred test style

- deterministic small synthetic inputs
- explicit expected results
- small feature-specific test executable
- minimal cross-feature coupling

## 12.2 For geometry features, use layers of tests

### Layer 1: tiny synthetic grid/result tests

Good for:
- mask logic
- metric formulas
- solver invariants
- separation enforcement

### Layer 2: small synthetic patch tests

Good for:
- Stage contracts
- provenance behavior
- reconstruction outcomes

### Layer 3: integration tests on real capsid input

Good for:
- artifact generation
- CLI wiring
- runtime sanity
- nonzero-output checks

## 12.3 What not to do

Do not rely only on one giant end-to-end test using `1cwp_full.vdb`. That makes debugging slower and less precise.

---

## 13. When to refactor instead of extending

Be honest here. Sometimes adding the feature “directly” is the wrong move.

### Refactor first if

- you need to touch the same helper category in multiple places
- you are adding a fourth or fifth unrelated writer function to one stage block
- you cannot test the new logic without going through huge orchestration
- you are adding multiple new method branches inside Stage 6 or Stage 7
- `geometry_analysis.cpp` becomes harder to navigate than the feature is worth

### The blunt truth

The current branch is already past the point where dumping everything into one implementation file is a good long-term idea.

---

## 14. A concrete roadmap for the next useful development cycle

If the goal is “add scientifically meaningful new functionality fast without wrecking the branch,” do this:

## Step 1

Create a dedicated metric module:
- `geometry_metrics`

## Step 2

Implement robust Stage-7-based scalar outputs:
- mean thickness
- min/max thickness
- patch coverage
- fraction of interpolated nodes
- fraction of adjusted nodes

## Step 3

Export:
- metrics summary CSV/JSON
- optional per-node thickness CSV

## Step 4

Add curvature on Stage 7 smoothed surfaces:
- mean curvature
- Gaussian curvature
- domain-restricted summary statistics

## Step 5

Only then consider deeper solver upgrades.

Reason: metrics unlock scientific use sooner than solver perfectionism.

---

## 15. Minimal checklist before coding

Use this every time.

1. What layer does this feature belong to?
2. Does it depend on frame/orientation assumptions?
3. Can it consume existing stage results instead of rebuilding state?
4. Can it reuse `geometry_symmetry`?
5. Does it need a new config struct?
6. Does it need a new result struct?
7. Does it need its own test target?
8. Does it require CLI wiring?
9. Does it need artifact export?
10. Is `geometry_analysis.cpp` still the right home, or should this be a new module?

---

## 16. Bottom-line development advice

The branch already gives you a strong backbone. Use it.

### Do this

- build on Stage 5–7 outputs
- keep frame semantics explicit
- preserve provenance
- return plain result structs
- keep algorithm code separate from CLI glue
- add small deterministic tests
- extract modules when logic starts repeating

### Do not do this

- duplicate fold math
- parse raw file text inside geometry work
- assume original frame silently
- mix metric formulas with file writing everywhere
- keep stuffing everything into `main.cpp`
- keep inflating `geometry_analysis.cpp` indefinitely

### The practical best next move

Add a **post-Stage-7 metrics module** and keep it read-only. That gives you real scientific output fast while staying aligned with the branch’s existing architecture.
