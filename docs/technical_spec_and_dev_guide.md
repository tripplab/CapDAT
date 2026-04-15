# CapDAT Geometry Branch — Updated Technical Specification Report

## Scope and basis

This report is based on direct inspection of the **`geometry` branch** codebase of `tripplab/CapDAT`, not on `main`, and not on the repository documentation alone. The in-repo technical guide is already partially stale relative to the current headers and CLI wiring, so this report is normalized against the actual source files.

Primary files inspected:

- `CMakeLists.txt`
- `src/main.cpp`
- `include/capsid.hpp`
- `include/geometry_symmetry.hpp`
- `include/reorientation_workflow.hpp`
- `include/export_capsid.hpp`
- `include/geometry_analysis.hpp`

This document is meant to be a **developer-facing technical map** of the current branch so new functionality can be designed against the real architecture instead of assumptions.

---

## 1. Executive technical state

The `geometry` branch is no longer just a parser plus a simple reorientation/export layer. It is now a **staged local shell-geometry analysis pipeline** for icosahedral capsids, organized around a fold-centered cylindrical patch extracted after canonical frame alignment.

The current pipeline is implemented through `runFoldPatchGeometryAnalysis(...)` and, as exposed by `include/geometry_analysis.hpp`, now spans **nine stages**:

1. **Stage 1** — preparation and fold alignment to `+Z`
2. **Stage 2** — cylindrical patch selection
3. **Stage 3** — analytical patch normalization with vdW assignment
4. **Stage 4** — raw inner/outer sheet detection by vertical line–sphere contact
5. **Stage 5** — surface-preparation masks, seeds, and trust domain setup
6. **Stage 6** — surface reconstruction
7. **Stage 7** — smoothing / regularization
8. **Stage 8** — local derivative estimation
9. **Stage 9** — curvature computation plus QC flagging

That matters because the geometry branch has already crossed from “preparation infrastructure” into a real numerical analysis subsystem. Any new work should treat it that way.

---

## 2. Build, executable, and test layout

## 2.1 Build system

The repository remains a **CMake-based C++17** project.

- project name: `CapDAT`
- version: `0.1.0`
- primary executable: `capsid_analyzer`
- tests controlled by `CAPDAT_BUILD_TESTS`
- warnings controlled by `CAPDAT_ENABLE_WARNINGS`

The geometry subsystem is first-class in the build graph: `src/geometry_analysis.cpp` is part of the main executable, not a side experiment.

## 2.2 Compiled source composition

The main executable is built from:

- orchestration:
  - `src/main.cpp`
- domain:
  - `src/atom.cpp`
  - `src/residue.cpp`
  - `src/chain.cpp`
  - `src/capsid.cpp`
- parsing:
  - `src/pdb_parser.cpp`
- symmetry / geometry foundation:
  - `src/geometry_symmetry.cpp`
- workflow:
  - `src/reorientation_workflow.cpp`
- export:
  - `src/export_capsid.cpp`
- reporting:
  - `src/structural_summary.cpp`
  - `src/summary_reporter.cpp`
- infrastructure:
  - `src/logger.cpp`
  - `src/timer.cpp`
- geometry pipeline:
  - `src/geometry_analysis.cpp`

## 2.3 Test targets visible in build

The branch defines dedicated test executables for:

- `capdat_structural_summary_tests`
- `capdat_geometry_symmetry_tests`
- `capdat_reorientation_workflow_tests`
- `capdat_geometry_analysis_tests`
- `capdat_export_capsid_tests`

There are also CLI-level tests for export and reorientation combinations.

The implication is clear: the intended development style is **modular, feature-local testing**, not only giant end-to-end runs.

---

## 3. High-level architecture

The branch is easiest to understand as six interacting layers.

## 3.1 CLI/orchestration layer

`src/main.cpp` is responsible for:

- parsing CLI arguments
- validating top-level mode combinations
- configuring parser behavior
- parsing the input file
- computing and printing structural summary
- optionally invoking standalone reorientation
- optionally invoking geometry analysis
- optionally exporting the final current capsid state

It is a coordinator, not the scientific core.

## 3.2 Persistent structural domain

The persistent biological structure remains the `Atom -> Residue -> Chain -> Capsid` hierarchy.

This is still the authoritative owner of accepted coordinates.

Important current rule: **geometry stages operate on the coordinates already stored in `Capsid`**, and Stage 1 may rotate them **in place**.

## 3.3 Canonical geometry/symmetry layer

`geometry_symmetry` is the repository’s reusable icosahedral math layer. It is parser-independent and should remain the only authoritative place for:

- canonical fold definitions
- fold lookup by name or `(type,index)`
- vector and matrix primitives
- direction alignment and rotation construction
- IAU boundary logic
- proper-rotation validation

## 3.4 Workflow mutation layer

`reorientation_workflow` handles user-facing coordinate reorientation as an explicit workflow with request/result models. It owns:

- source resolution
- target axis mapping
- transform computation
- in-place coordinate mutation
- orientation-state update semantics

The geometry pipeline reuses this workflow instead of duplicating rotation logic.

## 3.5 Export/serialization layer

`export_capsid` writes the **current in-memory accepted structure state**, not a verbatim copy of the source file. If the capsid was reoriented in place, export writes the transformed coordinates.

The writer also supports subset export through `atom_subset`, which is critical for patch- and contact-level geometry artifacts.

## 3.6 Derived geometry-analysis layer

`geometry_analysis` is the staged derived-analysis subsystem. It does not redefine the biological structure. Instead, it builds analysis-layer objects and result structs on top of the persistent capsid.

---

## 4. Core domain model relevant to geometry development

## 4.1 `Capsid` is the authoritative coordinate owner

`Capsid` owns the accepted parsed structure. Geometry stages traverse it directly. The branch does **not** create a second full capsid representation for rotated or analyzed states.

That is an intentional design choice and should be preserved unless there is an overwhelming reason to change it.

## 4.2 Explicit orientation/frame state exists and is central

`Capsid` now includes an embedded `OrientationState` with fields such as:

- `in_original_parsed_frame`
- `reoriented_in_place`
- `already_aligned_identity`
- `applied_rotation_matrix`
- `rotation_axis`
- `rotation_angle_radians`
- `source_mode`
- `source_description`
- `source_direction`
- `requested_target_axis`
- `target_direction`

This is one of the branch’s most important architectural safeguards. Any new feature that depends on frame semantics must respect it.

## 4.3 Persistent vs derived data is cleanly separated

The geometry pipeline introduces analysis-only structs instead of polluting `Atom` or `Capsid`:

- `CylinderMembership`
- `PatchAtom`
- `AnalyticalPatch`
- stage-specific result structs
- `GeometryAnalysisResult`

That separation is correct. Derived numerical state belongs in analysis-layer result models, not in the persistent biology objects.

---

## 5. `geometry_symmetry` as the canonical mathematical source

## 5.1 Public types

The module exposes:

- `Vector3`
- `Matrix3`
- `FoldDefinition`
- `RotationDefinition`
- `IauDefinition`

It also defines `RotationStatus` and `IauClassification`.

## 5.2 Canonical fold vocabulary

The geometry branch relies on the following canonical fold names:

- `2_0`
- `2_1`
- `3_0`
- `3_1`
- `5_0`

They are exposed through stable lookup functions and should be treated as the canonical branch vocabulary.

## 5.3 Functional scope

The module covers:

- fold lookup by exact name
- fold lookup by `(fold_type, fold_index)`
- exact reference vectors and normalized unit vectors
- norm, normalize, dot, cross
- angular relations and cosine similarity
- nearest-fold lookup
- Euclidean and point-to-axis distances
- direction-to-direction alignment
- alignment to `+Z`
- rotation of directions and points
- canonical IAU boundary margin logic
- proper-rotation validation

## 5.4 Development consequence

Do not hardcode fold vectors, axis logic, or custom direction-alignment math anywhere else. If a new feature needs reusable icosahedral math, it belongs here.

---

## 6. Reorientation workflow contract

## 6.1 Request model

`ReorientationRequest` contains:

- `request_reorientation`
- `source_mode`
- `fold_name`
- `custom_vector_text`
- `target_axis`
- `request_export`
- `export_path`
- `verbose`

## 6.2 Result model

`ReorientationResult` returns:

- `status`
- `resolved_source_description`
- `source_direction`
- `requested_target_axis`
- `target_direction`
- `rotation_matrix`
- `rotation_axis`
- `rotation_angle_radians`
- `coordinates_modified_in_place`
- `export_requested`
- `export_path`
- `messages`

## 6.3 Architectural role

This workflow owns the semantics of “take this user request and mutate the capsid frame accordingly.” The geometry-analysis layer should continue reusing it instead of re-implementing frame mutation.

---

## 7. Export subsystem contract

## 7.1 Export semantics

`ExportCapsidWriter` writes the **current accepted capsid coordinates** exactly as stored in memory.

That means:

- original parsed coordinates if no in-place transform was applied
- transformed coordinates if reorientation or geometry Stage 1 changed them

## 7.2 Configuration

`ExportCapsidConfig` includes:

- `output_path`
- `emit_header_comments`
- `emit_ter_records`
- `emit_end_record`
- `preserve_atom_serial_numbers`
- `coordinate_records_only`
- `atom_subset`

## 7.3 Why this matters

The branch already has the right foundation for exporting:

- full current capsid
- Stage 2 patch subsets
- Stage 4 contact-atom subsets
- future derived subsets

New functionality should reuse this path instead of building ad hoc PDB writers.

---

## 8. CLI surface as actually implemented

The real CLI is broader than the older geometry documentation suggests.

## 8.1 General options

- `--input`
- `--log`
- `--verbose`
- `--quiet`
- `--include-hetatm`
- `--export-final`

## 8.2 Reorientation options

- `--reorient`
- `--align-fold`
- `--align-vector`
- `--align-axis`

## 8.3 Geometry enablement and Stage 1–6 controls

- `--geometry-analysis`
- `--debug`
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

## 8.4 Stage 7 controls

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

## 8.5 Stage 8 controls

- `--geometry_s8_enabled`
- `--geometry_s8_fit_radius`
- `--geometry_s8_min_points`
- `--geometry_s8_max_points`
- `--geometry_s8_max_rms_residual`
- `--geometry_s8_max_abs_residual`
- `--geometry_s8_max_condition_indicator`
- `--geometry_s8_require_centered_support`
- `--geometry_s8_min_directional_span`
- `--geometry_s8_export_csv`

## 8.6 Stage 9 / curvature QC controls

- `--qc_n_tail`
- `--qc_n_spike`

## 8.7 Important top-level invariants enforced in `main.cpp`

The CLI rejects:

- missing `--input`
- simultaneous `--align-fold` and `--align-vector`
- `--reorient` without exactly one source
- source flags without `--reorient`
- simultaneous `--geometry-analysis` and `--reorient`
- invalid Stage 7 weights and solver settings
- invalid Stage 8 fit controls
- nonpositive Stage 9 QC thresholds

The design implication is good: `main.cpp` performs **mode consistency** and basic parameter sanity, while deeper stage-specific logic remains inside the geometry subsystem.

---

## 9. Geometry configuration model

The entire geometry workflow is driven by `FoldPatchAnalysisConfig`.

This struct is now more advanced than older documentation suggests.

## 9.1 Global / shared controls

- `enabled`
- `debug`
- `fold_type`
- `fold_index`
- `cylinder_radius`
- `delta_vdw`
- `grid_spacing`
- `min_atoms_in_patch`
- `output_prefix`
- `export_rotated_capsid`

## 9.2 Stage 5 controls

- `stage5_boundary_margin`
- `stage5_support_radius`
- `stage5_min_support_nodes`
- `stage5_reliable_radius`

## 9.3 Stage 6 controls

- `stage6_smoothing_weight`
- `stage6_max_iterations`
- `stage6_convergence_tolerance`
- `stage6_enforce_non_crossing`
- `stage6_min_separation`
- `stage6_export_obj_meshes`
- `stage6_mesh_export_format`
- `stage6_split_in_out_meshes`

## 9.4 Stage 7 controls

- `stage7_enabled`
- `stage7_method`
- `stage7_smoothing_weight`
- `stage7_max_iterations`
- `stage7_convergence_tolerance`
- `stage7_preserve_seed_values`
- `stage7_lambda`
- `stage7_data_weight_seed`
- `stage7_data_weight_interp`
- `stage7_use_reliable_core_only_for_fit`
- `stage7_boundary_condition_mode`
- `stage7_solver_max_iterations`
- `stage7_solver_tolerance`
- `stage7_export_s7s6_deltas`
- `stage7_enforce_non_crossing`
- `stage7_min_separation`
- `stage7_export_meshes`

## 9.5 Stage 8 controls

- `stage8_enabled`
- `stage8_fit_radius`
- `stage8_min_points`
- `stage8_max_points`
- `stage8_max_rms_residual`
- `stage8_max_abs_residual`
- `stage8_max_condition_indicator`
- `stage8_require_centered_support`
- `stage8_min_directional_span`
- `stage8_export_csv`

## 9.6 Stage 9 controls

- `stage9_enabled`
- `stage9_export_csv`
- `stage9_qc_n_tail`
- `stage9_qc_n_spike`
- `stage9_qc_min_neighbors`
- `stage9_qc_abs_scale_floor`

## 9.7 Design consequence

Any new behavior that materially changes the geometry pipeline should be surfaced through structured config, not buried in hidden constants.

---

## 10. Public geometry-analysis API surface

The header exposes both stage entry points and reusable kernels.

## 10.1 Reusable helpers and kernels

- `classifyPatchCylinder(...)`
- `normalizeElementSymbol(...)`
- `vdwRadius(...)`
- `makePatchAtom(...)`
- `intersectVerticalLineWithSphere(...)`
- `detectRawFirstContactAtNode(...)`
- `buildStage4RegularGrid(...)`

## 10.2 Stage entry points

- `prepareGeometryAnalysisStage1(...)`
- `runGeometryAnalysisStage2PatchSelection(...)`
- `runGeometryAnalysisStage3PatchNormalization(...)`
- `runGeometryAnalysisStage4RawSheetDetection(...)`
- `runGeometryAnalysisStage5SurfacePreparation(...)`
- `runGeometryAnalysisStage6SurfaceReconstruction(...)`
- `runGeometryAnalysisStage7SurfaceSmoothing(...)`
- `runGeometryAnalysisStage8DerivativeEstimation(...)`
- `runGeometryAnalysisStage9CurvatureComputation(...)`
- `runFoldPatchGeometryAnalysis(...)`

This is the real contract surface future development must respect.

---

## 11. Stage-by-stage current technical specification

## 11.1 Stage 1 — preparation and fold alignment

### Purpose

- resolve the requested canonical fold
- align it to `+Z`
- mutate the current `Capsid` coordinates in place
- capture the new frame state

### Inputs

- mutable `Capsid&`
- `FoldPatchAnalysisConfig`
- `ParserConfig`
- `Logger*`

### Outputs

`GeometryPreparationResult` includes:

- requested fold type/index
- resolved fold name
- resolved reference vector
- resolved unit vector
- identity-rotation flag
- rotation matrix
- rotation axis
- rotation angle
- in-place mutation flag
- final frame description
- export path
- nested `ReorientationResult`
- messages

### Architectural meaning

This stage establishes the working frame for all downstream local geometry. New stages should not silently rotate the capsid again.

---

## 11.2 Stage 2 — cylindrical patch selection

### Purpose

Select the fold-centered local patch from the working-frame capsid.

### Selection logic

For each accepted atom in current coordinates:

- compute `radial_xy`
- require `z > 0`
- require `radial_xy <= cylinder_radius`

Selection facts are encoded in `CylinderMembership`.

### Outputs

`GeometryPatchSelectionResult` stores:

- total atoms examined
- selected atom count
- rejection counters
- cylinder radius
- min-atoms threshold
- export path
- `patch_atoms`
- `selected_atom_refs`
- messages

### Important feature

The stage preserves stable pointers to original atoms. That is essential for provenance and export fidelity.

---

## 11.3 Stage 3 — analytical patch normalization

### Purpose

Transform the selected patch into a normalized analysis-layer representation.

### Main actions

For each Stage 2 selected atom:

- build a normalized `PatchAtom`
- normalize the element symbol
- assign vdW radius
- apply `delta_vdw`
- preserve selection membership facts
- preserve pointer to original atom

### vdW policy

The code supports element-based vdW resolution with a fallback strategy. The lookup table includes common elements and uses a fallback radius when needed.

### Outputs

`GeometryPatchNormalizationResult` contains an `AnalyticalPatch` with:

- normalized atoms
- original atom refs
- cylinder radius
- atom count
- explicit / inferred / fallback vdW counts
- export path

### Meaning

This is the true boundary between raw structural subset and analytically usable local geometry input.

---

## 11.4 Stage 4 — raw sheet detection

### Purpose

Generate raw outer and inner envelope evidence over a regular XY grid by intersecting vertical lines with patch-atom spheres.

### Grid model

`buildStage4RegularGrid(...)` builds a regular grid over the patch disk using `grid_spacing`.

### Node-level contact logic

At each grid node `(x,y)`:

- test vertical line–sphere intersections against all patch atoms
- outer raw contact is selected from the upper intersection side
- inner raw contact is selected from the lower intersection side
- validity depends on having a consistent inner/outer pair

### Outputs

`GeometryStage4RawSheetResult` includes:

- grid descriptor
- raw outer field
- raw inner field
- inside-disk mask
- valid mask
- candidate counts
- inner/outer contact serial numbers
- contact patch-atom indices
- per-atom contact roles
- raw contact records
- counts for valid, invalid, outer-only, inner-only, both-hit, zero-thickness, negative-thickness nodes
- unique-contact provenance counts
- CSV artifact paths
- contact-atoms PDB path
- timing metadata
- messages

### Why it matters

This is where the pipeline first turns the atom cloud into structured geometric evidence.

---

## 11.5 Stage 5 — surface preparation and trust-domain definition

### Purpose

Convert Stage 4 raw evidence into seeds, interpolation masks, exclusion masks, and the reliable core used later for smooth analysis.

### Main constructs

The stage produces:

- `z_outer_seed`
- `z_inner_seed`
- `inside_disk_mask`
- `raw_valid_mask`
- `outer_seed_mask`
- `inner_seed_mask`
- `paired_seed_mask`
- `boundary_exclusion_mask`
- `interp_allowed_outer_mask`
- `interp_allowed_inner_mask`
- `paired_interp_allowed_mask`
- `hard_invalid_mask`
- `reliable_core_mask`

### Derived parameters

If not explicitly set, the stage derives sensible geometry-dependent defaults such as boundary margin, support radius, and reliable radius.

### Outputs

`GeometryStage5SurfacePrepResult` includes all masks, seed fields, provenance sets, derived radii, artifact paths, and messages.

### Meaning

This stage defines the branch’s practical quality-control geometry. It is the trust boundary for later smoothing, derivatives, and curvature.

---

## 11.6 Stage 6 — surface reconstruction

### Purpose

Reconstruct continuous outer and inner scalar surfaces over the patch grid from Stage 5 seeds and interpolation-allowed nodes.

### Representation

The shell is represented as two scalar fields over XY:

- `z_outer_reconstructed`
- `z_inner_reconstructed`

### Outputs

`GeometryStage6SurfaceReconstructionResult` includes:

- reconstructed outer/inner fields
- masks for seeds, interpolation, hard-invalid, reconstructed, final-valid-analysis, non-crossing adjustment, and OBJ vertex usage
- counts for reconstructed, seed, interpolation, valid, adjusted, and unresolved nodes
- iteration statistics
- reconstructed separation statistics
- mesh counts
- CSV paths
- mesh paths
- messages

### Important limitation

This is still a regular-grid scalar-surface formulation, not a native manifold shell model.

---

## 11.7 Stage 7 — smoothing / regularization

### Purpose

Regularize the reconstructed Stage 6 surfaces into smoother, analysis-ready fields.

### Supported methods

The code currently supports:

- `smooth`
- `thin_plate_grid_fit`

### Outputs

`GeometryStage7SmoothedSurfaceResult` includes:

- smoothed outer/inner fields
- reconstructed/reliable-core/smooth-valid/metric-domain masks
- non-crossing adjustment mask
- counts for smooth-valid and metric-domain nodes
- iteration and solver statistics
- method labels and parameter echoes
- residual and bending-energy summaries
- Stage 7 vs Stage 6 delta summaries
- separation statistics
- semantic metadata such as:
  - outer normal orientation label
  - inner normal orientation label
  - thickness definition label
  - metric surface definition label
- CSV paths
- mesh paths and counts
- messages

### Critical observation

The result explicitly labels local thickness as `vertical_z_difference`. That is honest, but it also means physically rigorous thickness is still not fully solved here.

---

## 11.8 Stage 8 — derivative estimation

### Purpose

Estimate first and second derivatives of the Stage 7 smoothed surfaces on the metric domain.

### What the stage computes

For both outer and inner surfaces it stores:

- `dz/dx`
- `dz/dy`
- `d2z/dx2`
- `d2z/dy2`
- `d2z/dxdy`

It also tracks fit-quality metadata such as:

- RMS residual
- maximum absolute residual
- condition indicator
- neighbor count
- max support radius

### Validation masks

The stage explicitly distinguishes failure reasons through masks such as:

- insufficient points
- rank deficiency
- poor conditioning
- high residual
- bad boundary-neighbor geometry

### Outputs

`GeometryStage8DerivativeEstimationResult` includes derivative fields, fit diagnostics, validity masks, counts, CSV paths, and messages.

### Why this is a big deal

This means the branch now already has a real numerical differential-geometry bridge. It is no longer just reconstructing surfaces.

---

## 11.9 Stage 9 — curvature computation and QC

### Purpose

Compute curvature fields from Stage 8 derivatives and attach quality-control flags.

### Fields computed

For both outer and inner surfaces the stage stores:

- mean curvature
- oriented mean curvature
- Gaussian curvature
- graph-normal vs radial alignment metrics
- outward-normal alignment flags
- orientation-flip-applied flags

### QC outputs

The stage also produces QC masks and counts including:

- invalid nonfinite input/output
- global tail flags for Gaussian curvature
- local spike flags
- QC warning flags
- confidence classes

### Summaries

The result includes domain-level summaries such as:

- mean / median / stddev / min / max for `H`
- mean / median / stddev / min / max for oriented `H`
- mean / median / stddev / min / max for `K`
- curvature-valid fraction of metric domain
- curvature-valid fraction of derivative-valid nodes
- QC fractions

### Outputs

`GeometryStage9CurvatureComputationResult` includes curvature fields, QC flags, summary statistics, CSV paths, and messages.

### Consequence

Older branch narratives that say “curvature is not implemented yet” are now false. Curvature is present. What is still missing is a broader post-curvature scientific metric/reporting layer.

---

## 12. Result aggregation model

The top-level `GeometryAnalysisResult` aggregates:

- `preparation`
- `stage2_patch`
- `stage3_patch`
- `stage4_raw`
- `stage5_prep`
- `stage6_surfaces`
- `stage7_smooth`
- `stage8_derivatives`
- `stage9_curvature`
- `messages`

This is the branch’s true high-level analysis contract.

---

## 13. Artifact/export model of the geometry pipeline

The geometry branch is artifact-heavy by design, which is a strength for debugging numerical work.

## 13.1 PDB artifacts

The branch can emit:

- rotated capsid export
- Stage 2 patch-atom subset PDB
- Stage 4 contact-atom subset PDB
- final current capsid export via `--export-final`

## 13.2 CSV artifacts

The stage results expose paths for many CSV outputs, including:

- normalized atoms
- raw outer/inner sheets
- seed fields
- masks
- reconstructed fields
- smoothed fields
- delta maps
- derivative fields
- curvature fields
- stage-level summaries

## 13.3 Mesh artifacts

The branch supports mesh export for reconstructed and smoothed surfaces, with format control and split inner/outer options.

## 13.4 Design meaning

The branch is built around inspectable intermediate state. That is exactly what you want in a geometry-analysis codebase.

---

## 14. Current numerical design characteristics

## 14.1 The pipeline is still grid-first

The local shell is analyzed on a regular XY grid aligned to the selected fold axis after Stage 1.

## 14.2 The shell representation is two scalar sheets

The shell is represented as:

- outer height map
- inner height map

This is simple and practical, but it is not a general-purpose shell manifold representation.

## 14.3 Derivatives and curvature are local grid-based numerical estimates

Stage 8 and Stage 9 show that the branch has moved into differential geometry, but still in the context of local polynomial/grid-based estimation, not a more general mesh-native DG framework.

## 14.4 Thickness remains semantically limited

Even with Stage 9 present, the Stage 7 metadata still says the local thickness definition is vertical Z difference. That remains a conceptual limitation.

## 14.5 QC is explicit but heuristic

The branch now has explicit masks and thresholds for trust-domain construction, derivative validity, and curvature outlier/spike detection. That is good engineering, but it is still heuristic QC, not uncertainty propagation.

---

## 15. What is structurally strong right now

The following parts are strong enough to build on directly:

## 15.1 Frame handling is explicit

The combination of `Capsid::OrientationState`, `geometry_symmetry`, and `reorientation_workflow` is solid and should be reused.

## 15.2 Stage boundaries are real

The staged result structs are explicit, data-rich, and suitable as extension points.

## 15.3 Provenance is treated seriously

The code preserves atom pointers, contact provenance, and validity masks instead of collapsing everything into anonymous numeric fields.

## 15.4 The pipeline already reaches curvature

The branch is further along than earlier specs implied. It already reconstructs surfaces, regularizes them, estimates derivatives, and computes curvature with QC.

## 15.5 Artifact-oriented debugging is built in

That is one of the most valuable architectural features for future numerical development.

---

## 16. Current limitations and technical debt

Here is the blunt version.

## 16.1 `geometry_analysis.cpp` is doing too much

The geometry subsystem is still concentrated in one very large implementation file. That file now likely mixes:

- stage orchestration
- numerical kernels
- validation
- mesh writing
- CSV writing
- provenance logic
- summary serialization

That is survivable in the short term and a liability in the medium term.

## 16.2 Final scientific metric synthesis is still incomplete

The branch now has curvature, but it still lacks a cleaner top-level post-analysis layer for results such as:

- thickness summaries beyond vertical gap semantics
- patch area metrics
- consolidated QC reports
- one coherent final analysis summary object for publication-facing outputs

## 16.3 The shell model is still local and graph-based

Everything assumes a local height-field representation relative to the selected axis. That is fine for small, well-behaved patches and limiting for more general geometry.

## 16.4 Mesh export is downstream utility, not geometric truth

The mesh outputs are useful artifacts, but the underlying truth representation remains the scalar grid fields.

## 16.5 Stage naming in external docs is behind the code

The repository docs lag the actual implementation. Future development should trust the headers and code first.

---

## 17. Natural extension points for new functionality

## 17.1 Post-Stage-9 reporting layer

The cleanest next addition is a dedicated module that consumes Stage 7/8/9 outputs and produces a final scientific summary:

- thickness summaries
- area estimates
- reliability summaries
- consolidated geometry reports

## 17.2 Better thickness definitions

A high-value next step is to add:

- normal-projected thickness
- explicit distinction between vertical and normal-based thickness
- summary comparison between both

## 17.3 Additional metric-domain QC

The branch already contains the ingredients for stronger confidence measures:

- interpolation fraction
- derivative-valid fraction
- curvature-valid fraction
- outlier/spike burden
- non-crossing adjustment burden

## 17.4 Solver and representation upgrades

Possible future directions include:

- better Stage 6/7 solvers
- extracted sparse-system backends
- anisotropic smoothing
- confidence-weighted fitting
- eventually mesh-native analysis

---

## 18. Recommended refactor boundaries before major expansion

If this branch will keep growing, the following split is worth doing early:

- `geometry_patch_selection.*`
- `geometry_patch_normalization.*`
- `geometry_raw_sheet.*`
- `geometry_surface_prep.*`
- `geometry_surface_reconstruction.*`
- `geometry_surface_smoothing.*`
- `geometry_derivatives.*`
- `geometry_curvature.*`
- `geometry_surface_io.*`
- `geometry_mesh_export.*`

The main problem is not that the current code is conceptually wrong. The problem is that continued growth in one implementation file will make correctness harder to maintain.

---

## 19. Bottom-line technical conclusion

The current `geometry` branch is best described as a **fold-aligned local shell reconstruction and differential-geometry analysis pipeline** on top of the CapDAT structural core.

It already has:

- canonical fold math
- explicit frame mutation semantics
- cylindrical patch extraction
- vdW-normalized analytical patch representation
- raw sheet detection
- trust-domain preparation
- surface reconstruction
- surface regularization
- derivative estimation
- curvature computation with QC
- artifact exports across multiple stages

It still lacks:

- a cleaner final scientific metric/report layer
- stronger thickness definitions
- cleaner module decomposition
- a less monolithic geometry implementation layout

That is the honest state. The branch is already substantial and usable as a development base, but it is also far enough along that sloppy extension work will create avoidable architecture debt fast.





# CapDAT Geometry Branch — Synthesized Development Guide for Implementing New Functionality

## 1. Start by classifying the feature correctly

Most bad implementations begin with putting the feature in the wrong layer.

## Put it in the parser only if

- the new information comes directly from raw file records
- the acceptance policy changes
- structural reconstruction rules change
- the data belongs to the persistent biological model itself

## Put it in `Capsid` or the persistent domain only if

- the state is authoritative and long-lived
- multiple workflows need it as structural truth
- it is not merely temporary analysis output

## Put it in `geometry_symmetry` if

- it is reusable icosahedral math
- it involves fold lookup, axis geometry, direction comparison, IAU logic, or generic rotations
- more than one future module could need it

## Put it in `reorientation_workflow` if

- it is a coordinate-transform workflow
- it mutates the current frame
- it needs request/result semantics
- it must update `Capsid::OrientationState`

## Put it in geometry analysis or a sibling geometry module if

- it operates on accepted coordinates after parsing
- it is patch logic, surface logic, derivative logic, curvature logic, or metric logic
- it is derived analysis rather than parsing

## Put it in export code only if

- the logic is pure serialization
- it writes already-computed state to files
- it should not decide the scientific algorithm

That classification step is not optional. If you skip it, you will place logic incorrectly and spend the rest of the work compensating for that mistake.

---

## 2. Respect the current workflow instead of inventing a parallel one

The branch already has a real staged analysis pipeline. Reuse it.

## Current geometry workflow

1. parse input file into `Capsid`
2. compute structural summary
3. optionally run standalone reorientation
4. or run geometry analysis:
   1. Stage 1 align requested fold to `+Z`
   2. Stage 2 select cylindrical patch
   3. Stage 3 normalize analytical patch
   4. Stage 4 detect raw inner/outer sheets
   5. Stage 5 prepare seeds, masks, and reliable core
   6. Stage 6 reconstruct surfaces
   7. Stage 7 smooth / regularize
   8. Stage 8 estimate derivatives
   9. Stage 9 compute curvature and QC
5. optionally export current full capsid state

## Non-negotiable rule

If your feature depends on fold-centered local geometry, do **not** rebuild preprocessing from scratch. Start from the existing stage outputs unless there is a hard technical reason not to.

---

## 3. Know the safest extension points

## 3.1 Best place for new scientific metrics

The cleanest current insertion point is **after Stage 9**.

Why:

- Stage 7 already gives smoothed surfaces
- Stage 8 gives derivatives
- Stage 9 gives curvature and QC masks
- Stage 5 still provides the reliable-core and trust-domain semantics

That means a new scientific metrics layer should usually consume:

- `GeometryStage7SmoothedSurfaceResult`
- `GeometryStage8DerivativeEstimationResult`
- `GeometryStage9CurvatureComputationResult`
- possibly `GeometryStage5SurfacePrepResult`

## Best module shape

- `include/geometry_metrics.hpp`
- `src/geometry_metrics.cpp`
- `tests/geometry_metrics_tests.cpp`

## Why this is the best next move

Because it adds real scientific value without destabilizing the established Stage 1–9 pipeline.

---

## 3.2 Best place for new smoothing / regularization methods

Add them to Stage 7 only if they still conceptually belong to “surface regularization over the existing grid fields.”

Current Stage 7 already supports:

- `smooth`
- `thin_plate_grid_fit`

If the new method still produces the same kind of outer/inner smoothed scalar fields, extend Stage 7.

If the method changes representation too much, do not force it into Stage 7. Create a sibling module.

---

## 3.3 Best place for new differential-geometry refinements

If the change affects:

- local derivative fitting
- support geometry for derivatives
- derivative quality filters
- derivative-confidence logic

then it belongs in Stage 8 or a dedicated sibling derivative module.

If the change affects:

- curvature formulas
- orientation conventions
- curvature QC logic
- confidence labeling for curvature

then it belongs in Stage 9 or a curvature sibling module.

---

## 4. Never ignore frame semantics

This branch already solved one major class of bugs by making frame state explicit. Do not reintroduce them.

## Always know which frame your feature assumes

Possible assumptions are:

- original parsed frame
- current arbitrary in-memory frame
- fold-aligned Stage 1 working frame
- downstream geometry-analysis frame after Stage 1

## The authoritative source

Use `Capsid::orientationState()`.

## Rule

If your feature assumes the selected fold is aligned to `+Z`, make that precondition explicit. Do not silently assume it.

## Practical pattern

For new local geometry work, require one of:

- successful Stage 1 result
- or later stage results that already guarantee the working frame

Anything else is asking for garbage geometry.

---

## 5. Reuse existing result contracts instead of inventing new containers

The branch already gives you the main data boundaries you need.

## Existing useful contracts

- `PatchAtom`
- `AnalyticalPatch`
- `GeometryPreparationResult`
- `GeometryPatchSelectionResult`
- `GeometryPatchNormalizationResult`
- `GeometryStage4RawSheetResult`
- `GeometryStage5SurfacePrepResult`
- `GeometryStage6SurfaceReconstructionResult`
- `GeometryStage7SmoothedSurfaceResult`
- `GeometryStage8DerivativeEstimationResult`
- `GeometryStage9CurvatureComputationResult`
- `GeometryAnalysisResult`

## Do not do this

Do not create random parallel containers for:

- patch atoms
- masks
- derivative maps
- curvature maps
- reliable-domain logic
- provenance data

unless the existing result objects truly cannot support the feature.

Most duplicate containers are just unacknowledged design drift.

---

## 6. Physical placement of new code

## Good pattern for a new read-only analysis feature

Use:

- `include/<feature>.hpp`
- `src/<feature>.cpp`
- `tests/<feature>_tests.cpp`

Example:

- `include/geometry_metrics.hpp`
- `src/geometry_metrics.cpp`
- `tests/geometry_metrics_tests.cpp`

## Good pattern for helper extraction

If repeated logic appears across geometry stages, extract helpers into focused internal modules such as:

- `geometry_grid_utils.*`
- `geometry_surface_io.*`
- `geometry_mesh_export.*`
- `geometry_derivative_utils.*`
- `geometry_curvature_utils.*`

## Bad pattern

Do not keep adding unrelated blocks to `geometry_analysis.cpp` forever. That file is already a future maintenance trap.

---

## 7. Coding style the branch clearly prefers

Match the branch style or your code will feel bolted on.

## The visible style is

- plain config structs
- plain result structs
- explicit stage entry functions
- explicit validation
- deterministic control flow
- logger-based runtime messaging
- artifact paths carried in result objects
- provenance preserved where possible

## That means

Good additions look like:

- `FeatureConfig`
- `FeatureResult`
- `runFeature(...)`

Bad additions look like:

- hidden state in constructors
- giant side effects
- implicit global behavior
- algorithm logic embedded directly in CLI parsing
- undocumented mutations to persistent domain objects

The codebase is not built around OO cleverness. It is built around explicit structured workflow.

---

## 8. Correct design pattern for a new metrics module

This is likely the next high-value extension area, so do it cleanly.

## 8.1 Inputs

Take only the upstream data you really need.

Typical read-only input set:

- `const GeometryStage7SmoothedSurfaceResult& stage7`
- `const GeometryStage8DerivativeEstimationResult& stage8`
- `const GeometryStage9CurvatureComputationResult& stage9`
- `const GeometryStage5SurfacePrepResult& stage5`
- `const FoldPatchAnalysisConfig& config`
- `Logger* logger`

## 8.2 Config

Create a dedicated config if the feature needs its own knobs.

Do **not** automatically stuff every new parameter into `FoldPatchAnalysisConfig`.

Add new options to `FoldPatchAnalysisConfig` only when:

- the parameter is core to the geometry-analysis pipeline itself
- it must be controlled directly at the top-level geometry CLI

Otherwise, a feature-local config is cleaner.

## 8.3 Result

Return a plain result struct containing:

- success flag
- scalar summaries
- optional per-node fields
- optional QC summaries
- optional artifact paths
- messages

## 8.4 Export strategy

If the feature writes CSVs or summaries, keep the computational logic separate from the file formatting logic as much as possible.

Do not mix formulas and serialization into one giant untestable block.

---

## 9. Specific advice for likely next features

## 9.1 Better thickness metrics

This is the most obvious remaining scientific gap.

## Wrong move

Treat `z_outer - z_inner` as the final thickness definition for everything.

That is easy and incomplete.

## Better move

Implement both:

- `vertical_thickness`
- `normal_projected_thickness` or an explicitly improved shell-thickness estimate

Then report them separately and say exactly what each means.

## Good first deliverable

A read-only metrics module that computes over the metric domain:

- mean vertical thickness
- median thickness
- min/max thickness
- robust spread
- valid coverage fraction
- interpolation burden
- non-crossing-adjustment burden

That is useful immediately and does not require a representation rewrite.

---

## 9.2 Patch area metrics

This is another clean addition.

Useful distinctions:

- projected XY patch area
- reconstructed outer surface area
- reconstructed inner surface area
- maybe mean-surface area if defined later

Do not conflate them. A projected disk area is not the same thing as reconstructed shell area.

---

## 9.3 Final QC synthesis

The branch already computes many partial QC signals. A good next feature is to consolidate them.

Good candidates:

- metric-domain fraction of inside-disk nodes
- derivative-valid fraction
- curvature-valid fraction
- Stage 7 non-crossing adjustment fraction
- Stage 8 fit-failure reason burden
- Stage 9 tail/spike burden

That gives users a practical “can I trust this patch?” summary instead of making them inspect many CSVs manually.

---

## 9.4 New Stage 7 methods

If you add a new Stage 7 method:

1. extend `FoldPatchAnalysisConfig::Stage7Method`
2. add CLI parsing in `main.cpp`
3. keep the same Stage 7 result contract
4. implement a dedicated internal solver function
5. wire it through `runGeometryAnalysisStage7SurfaceSmoothing(...)`
6. add targeted tests

## Non-negotiable rule

Do not invent a method-specific ad hoc output shape if it still belongs to Stage 7. Keep the stage contract stable.

---

## 9.5 New derivative estimation logic

If you improve Stage 8:

- preserve the explicit failure masks
- preserve fit diagnostics
- preserve nodewise validity semantics

The current Stage 8 contract is already good because it distinguishes why a node failed. Do not collapse that into one generic “invalid” flag.

---

## 9.6 New curvature logic

If you improve Stage 9:

- keep outward/inward normal semantics explicit
- keep orientation-flip tracking explicit
- preserve QC flag layers
- preserve summary statistics

Do not silently change curvature sign conventions without exposing that clearly in result metadata.

---

## 10. Use `geometry_symmetry` aggressively

Whenever you need:

- fold lookup
- fold direction
- angle to fold
- point-to-axis distance
- direction alignment
- IAU logic
- proper-rotation checks

use `geometry_symmetry`.

## Rule

If you are hardcoding fold vectors or writing custom fold-axis math in a new feature, you are almost certainly doing it wrong.

---

## 11. How to wire a new CLI-controlled feature cleanly

Follow the existing `main.cpp` pattern.

## Correct steps

1. add local CLI variables near related options
2. parse the new flag in the argument loop
3. validate only high-level invariants in `main.cpp`
4. map CLI values into config structs
5. keep all real feature logic out of `main.cpp`
6. call a dedicated module function after prerequisites are ready

## Wrong steps

Do not:

- put numerical kernels inside `main.cpp`
- build one-off ad hoc execution branches with hidden semantics
- bypass the stage pipeline because it seems faster

That shortcut becomes permanent technical debt.

---

## 12. Testing strategy for new features

The branch already tells you the preferred testing model.

## Use three layers

### Layer 1 — small deterministic unit tests

Best for:

- formulas
- mask logic
- derivative kernels
- curvature sign conventions
- solver invariants

### Layer 2 — small synthetic geometry tests

Best for:

- stage behavior
- provenance preservation
- failure masks
- support-geometry logic

### Layer 3 — integration tests on realistic input

Best for:

- artifact generation
- CLI wiring
- end-to-end sanity
- nonzero result checks

## What not to do

Do not rely only on one giant real-capsid test. It is slow, opaque, and terrible for debugging numerical regressions.

---

## 13. When to refactor before extending

Be honest about this. Sometimes adding the feature “directly” is the wrong move.

## Refactor first if

- you need to touch the same helper logic in multiple stages
- you are adding another large writer block to `geometry_analysis.cpp`
- you cannot test the new logic without running huge orchestration
- Stage 7 or Stage 8 branching is getting crowded
- you are about to add a second or third related feature into the same giant file section

## Blunt truth

The branch is already beyond the point where “just one more block in `geometry_analysis.cpp`” is a healthy long-term habit.

---

## 14. Practical roadmap for the next useful development cycle

If the goal is to add value fast without damaging the architecture, do this:

## Step 1

Create a dedicated post-Stage-9 metrics module.

## Step 2

Implement final scalar summary outputs such as:

- vertical thickness summaries
- patch area summaries
- derivative-valid fraction
- curvature-valid fraction
- QC burden summaries

## Step 3

Export:

- metrics summary CSV
- metrics summary JSON
- optional per-node thickness or QC maps

## Step 4

Only then add improved thickness definitions based on normals.

## Step 5

After that, consider deeper solver upgrades if they are still necessary.

Why this order matters: it gives scientifically useful outputs sooner and avoids getting trapped in solver perfectionism before the branch exposes the final metrics users actually need.

---

## 15. Minimal checklist before writing code

Use this every time.

1. What layer does this feature belong to?
2. What frame assumptions does it make?
3. Can it consume an existing stage result instead of rebuilding state?
4. Should it reuse `geometry_symmetry`?
5. Does it need its own config struct?
6. Does it need its own result struct?
7. Does it need artifact export?
8. Does it need CLI wiring?
9. Does it need a dedicated test target?
10. Is `geometry_analysis.cpp` really the right home, or should this be a new module?

If you cannot answer those clearly, you are not ready to implement the feature cleanly.

---

## 16. Bottom-line implementation advice

## Do this

- build on the existing stage outputs
- keep frame semantics explicit
- preserve provenance
- return plain result structs
- separate algorithm logic from CLI glue
- keep exports structured and inspectable
- add deterministic tests
- extract modules once a helper category repeats

## Do not do this

- duplicate fold math
- parse raw coordinate text inside geometry features
- assume the working frame silently
- hide scientific semantics in export code
- stuff real algorithm logic into `main.cpp`
- keep inflating `geometry_analysis.cpp` without restraint

## Best immediate next move

Build a **read-only post-Stage-9 metrics/reporting module**. It is the highest-value extension that matches the current architecture, uses the existing Stage 1–9 pipeline correctly, and closes the biggest remaining gap between numerical geometry and user-facing scientific output.




