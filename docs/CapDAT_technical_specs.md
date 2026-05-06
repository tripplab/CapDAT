# CapDAT `geometry` Branch — Detailed Technical Specification

## 2. What the repository is

CapDAT is a C++17 command-line toolkit for parsing and analyzing large viral capsid coordinate files in fixed-column PDB-like format, including VIPERdb-style `.vdb` inputs. The `geometry` branch is a **post-parse fold-centered local shell geometry pipeline** built on top of a structural core, with explicit frame handling, deterministic icosahedral fold math, in-place reorientation, canonical export, and a staged numerical analysis path from patch extraction through curvature and thickness estimation.

In plain terms, the repo does three big things:

1. It reconstructs an internal structural model of an icosahedral capsid from raw coordinate records.
2. It can reorient that model in place into a fold-aligned working frame.
3. It can derive a local geometric patch around a chosen fold and run a staged analysis on that patch to reconstruct surfaces, estimate derivatives, compute curvature, and compute thickness.

## 3. Build and executable structure

The build is CMake-based, uses C++17, and produces a single main executable named `capsid_analyzer` plus dedicated test executables. The main source list includes the structural domain, parser, symmetry module, export module, reorientation workflow, summary/reporting, infrastructure, and the full geometry pipeline through `src/geometry_analysis.cpp`.

The build graph makes several architectural facts obvious:

- geometry is first-class, not optional dead code
- the geometry pipeline is linked into the main executable
- tests are modularized by subsystem
- the expected development workflow includes CTest, not only manual runs

The test executables explicitly include:

- `capdat_structural_summary_tests`
- `capdat_geometry_symmetry_tests`
- `capdat_reorientation_workflow_tests`
- `capdat_geometry_analysis_tests`
- `capdat_export_capsid_tests`

and CLI-level tests for export and reorientation combinations :contentReference[oaicite:18]{index=18}.

## 4. High-level architecture

The repo is best understood as six interacting layers.

### 4.1 CLI orchestration layer

`src/main.cpp` is the coordinator. It does argument parsing, top-level mode validation, parser configuration, parse execution, structural summary printing, optional reorientation, optional geometry analysis, and optional final export :contentReference[oaicite:19]{index=19}.

It is intentionally not the place where scientific kernels live. The CLI layer maps user intent into structured configs and workflow calls.

### 4.2 Persistent structural domain layer

The persistent biological/structural hierarchy is:

**Atom → Residue → Chain → Capsid**

This is the authoritative owner of accepted coordinates. Geometry code does not create a second full structural model. It consumes the `Capsid` and derived stage outputs built from it :contentReference[oaicite:20]{index=20} :contentReference[oaicite:21]{index=21} :contentReference[oaicite:22]{index=22}.

Important nuance: `Chain` is only the implementation class name. Conceptually it represents an internally reconstructed subunit, not a unique raw one-letter PDB chain identity. This matters because large capsids can reuse raw chain IDs across multiple proteins, so the code assigns its own internal subunit identity and only preserves the original chain label as metadata :contentReference[oaicite:23]{index=23}.

### 4.3 Parsing and acceptance-policy layer

`PdbParser` is the only layer that should understand the raw fixed-column layout. It owns record extraction, temporary record validation, accept/reject policy, and hierarchical reconstruction into `Capsid` :contentReference[oaicite:24]{index=24}.

The parser configuration already encodes an important domain policy: `protein_only = true` by default, meaning the analysis basis excludes non-protein components unless policy changes :contentReference[oaicite:25]{index=25}.

### 4.4 Canonical geometry/symmetry layer

`geometry_symmetry` is the repository’s canonical icosahedral math module. It is parser-independent and exposes fold definitions, normalized fold vectors, vector algebra, distance and angle helpers, direction-alignment rotations, IAU classification, and proper-rotation checks :contentReference[oaicite:26]{index=26}.

This module is the single source of truth for fold-axis math.

### 4.5 Workflow mutation layer

`reorientation_workflow` owns the semantics of applying a user-requested transform to the current capsid frame. It translates a structured request into source resolution, axis mapping, rotation calculation, in-place coordinate updates, and workflow result reporting :contentReference[oaicite:27]{index=27}.

### 4.6 Export and artifact layer

`export_capsid` writes the **current in-memory accepted structure state**, not the original file as text. If the capsid has been reoriented in place, export writes those transformed coordinates. It also supports subset export through `atom_subset`, which is crucial for patch- and contact-level debug artifacts :contentReference[oaicite:28]{index=28}.

### 4.7 Derived geometry-analysis layer

`geometry_analysis` is a staged derived-analysis subsystem built on top of the persistent capsid. It does not mutate topology. It constructs stage outputs, scalar fields, masks, metrics, artifacts, and reportable summaries from the current accepted coordinates :contentReference[oaicite:29]{index=29}.

## 5. Core persistent model semantics

### 5.1 Atom

`Atom` is intentionally lightweight. It stores parsed PDB-derived identity and coordinate information and supports only a small mutator, `setPosition(...)`, for post-parse coordinate transforms like reorientation :contentReference[oaicite:30]{index=30}.

### 5.2 Chain

`Chain` caches atom counts, stores a unique internal subunit ID plus the original PDB chain label, owns residues, and supports controlled mutation through append-oriented methods plus mutable residue access for post-parse transformations :contentReference[oaicite:31]{index=31}.

### 5.3 Capsid

`Capsid` is the authoritative structure owner. It stores source metadata, owns all reconstructed subunits, stores summary counters, and now also stores an explicit `OrientationState` that records whether the in-memory coordinates are still in the original parsed frame or have been reoriented in place :contentReference[oaicite:32]{index=32}.

This orientation state is architecturally central. Downstream code is not supposed to assume the frame.

## 6. Frame and orientation semantics

A major strength of this branch is that frame identity is explicit instead of implicit.

`Capsid::OrientationState` tracks, among other things:

- whether the structure is still in the original parsed frame
- whether a reorientation was applied in place
- whether the transform was identity because the requested fold was already aligned
- the applied rotation matrix
- rotation axis and angle
- source mode and source description
- source direction
- requested target axis
- target direction

This means frame-aware features have an authoritative state object to consult instead of guessing from workflow history :contentReference[oaicite:33]{index=33}.

## 7. Canonical symmetry and fold model

The `geometry_symmetry` module defines lightweight geometry primitives and canonical fold definitions. The canonical fold vocabulary visible in the public API is:

- `2_0`
- `2_1`
- `3_0`
- `3_1`
- `5_0`

The module exposes lookup by exact name and by `(fold_type, fold_index)`, plus helpers for normalized directions, nearest-fold identification, distances to axes, alignment to arbitrary target directions, alignment to `+Z`, and IAU classification functions :contentReference[oaicite:34]{index=34}.

Architecturally, this means that custom fold hardcoding elsewhere is a design violation.

## 8. Reorientation workflow contract

The reorientation layer is driven by `ReorientationRequest`, which contains:

- `request_reorientation`
- `source_mode` as fold or custom vector
- `fold_name`
- `custom_vector_text`
- `target_axis`
- optional export request fields
- verbosity flag

It returns `ReorientationResult`, which contains:

- status
- resolved source description
- source direction
- target direction
- rotation matrix
- rotation axis
- rotation angle
- whether coordinates were modified in place
- export state
- messages

This is important because `main.cpp` does not perform deep source resolution itself. It only wires CLI state into this request model and delegates the real semantics to the workflow layer :contentReference[oaicite:35]{index=35} :contentReference[oaicite:36]{index=36}.

## 9. Export subsystem contract

The export writer is not a parser reverse-dumper. It serializes the current accepted in-memory coordinates and can also export subsets of atoms. The config includes:

- output path
- header/comment emission
- TER/END record control
- whether to preserve serial numbers
- coordinate-record-only mode
- optional explicit atom subset pointer

This is what makes Stage 2 patch exports and Stage 4 contact-atom exports possible without creating ad hoc writers :contentReference[oaicite:37]{index=37}.

## 10. Actual CLI surface on the `geometry` branch

The real CLI surface is broader than older branch descriptions imply.

### General
- `--input`
- `--log`
- `--verbose`
- `--quiet`
- `--include-hetatm`
- `--export-final`

### Reorientation
- `--reorient`
- `--align-fold`
- `--align-vector`
- `--align-axis`

### Geometry analysis enablement and early-stage controls
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

### Stage 7 smoothing / regularization
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

### Stage 8 derivative estimation
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

### Stage 9 curvature QC
- `--qc_n_tail`
- `--qc_n_spike`

### Stage 10 thickness
- `--geometry_thickness_enabled`
- `--geometry_thickness_method`
- `--geometry_thickness_export_csv`
- `--geometry_thickness_use_curvature_valid_domain_only`
- `--geometry_thickness_min_thickness`
- `--geometry_thickness_max_thickness`

All of that is visible directly in `printHelp(...)` and the argument parser in `main.cpp` :contentReference[oaicite:38]{index=38}.

Current mesh format options are `obj`, `stl`, and `ply`, with `ply` as the default format for Stage 6/7 mesh export in the CLI/config surface.

## 11. Top-level CLI invariants and execution rules

`main.cpp` enforces several important invariants before execution:

- `--input` is mandatory
- `--align-fold` and `--align-vector` are mutually exclusive
- `--reorient` requires exactly one source
- source flags cannot appear without `--reorient`
- `--geometry-analysis` cannot be combined with `--reorient`
- various Stage 7, 8, 9, and 10 numerical parameters must satisfy positivity or ordering constraints

That split is good engineering. High-level mode consistency lives in `main.cpp`; deeper stage semantics belong in the geometry subsystem :contentReference[oaicite:39]{index=39}.

## 12. Geometry pipeline: public API and staged contract

`include/geometry_analysis.hpp` exposes a formal staged API. This is the real contract for geometry development.

### Helper/kernel API
- `classifyPatchCylinder(...)`
- `normalizeElementSymbol(...)`
- `vdwRadius(...)`
- `makePatchAtom(...)`
- `intersectVerticalLineWithSphere(...)`
- `detectRawFirstContactAtNode(...)`
- `buildStage4RegularGrid(...)`

### Stage entry points
- `prepareGeometryAnalysisStage1(...)`
- `runGeometryAnalysisStage2PatchSelection(...)`
- `runGeometryAnalysisStage3PatchNormalization(...)`
- `runGeometryAnalysisStage4RawSheetDetection(...)`
- `runGeometryAnalysisStage5SurfacePreparation(...)`
- `runGeometryAnalysisStage6SurfaceReconstruction(...)`
- `runGeometryAnalysisStage7SurfaceSmoothing(...)`
- `runGeometryAnalysisStage8DerivativeEstimation(...)`
- `runGeometryAnalysisStage9CurvatureComputation(...)`
- `runGeometryAnalysisStage10ThicknessComputation(...)`
- `runFoldPatchGeometryAnalysis(...)`
- `buildGeometrySummaryReport(...)`

The pipeline is therefore a formalized sequence, not an informal script :contentReference[oaicite:40]{index=40}.

## 13. Stage-by-stage technical meaning

### Stage 1 — preparation and fold alignment

This stage resolves the requested canonical fold, aligns it to `+Z`, mutates the `Capsid` in place when necessary, and records the resulting frame semantics. Its result includes the fold name, reference vector, unit vector, rotation matrix, axis, angle, whether identity rotation was used, whether coordinates were modified, final frame description, and nested reorientation workflow result :contentReference[oaicite:41]{index=41}.

Tests show that identity alignment for `2_0` is treated explicitly and recorded semantically rather than being silently ignored, while non-identity alignments such as `3_0` do change coordinates in place :contentReference[oaicite:42]{index=42}.

### Stage 2 — cylindrical patch selection

This stage selects the local patch in the working frame. It evaluates each accepted atom by:

- radial distance in the XY plane
- positivity of `z`
- inclusion within the cylinder radius

The result preserves both `PatchAtom` records and stable original atom pointers, plus counters for rejected atoms and the patch export path :contentReference[oaicite:43]{index=43}.

Tests confirm that `z > 0` and `radial_xy <= cylinder_radius` are the selection rule, that edge inclusion is deterministic, and that min-patch-size thresholds are enforced as hard failures :contentReference[oaicite:44]{index=44}.

### Stage 3 — analytical patch normalization

Stage 3 turns the selected subset into an `AnalyticalPatch`. For each selected atom it normalizes element identity, assigns vdW radius, applies global `delta_vdw`, preserves membership facts, and keeps original atom provenance. Counts are tracked for explicit, inferred, and fallback vdW assignments :contentReference[oaicite:45]{index=45}.

Tests confirm explicit lookup for common elements, inference from atom names when the element field is blank, and fallback behavior when the element remains unsupported :contentReference[oaicite:46]{index=46}.

### Stage 4 — raw inner/outer sheet detection

Stage 4 introduces the first real geometric envelope construction. It builds a regular XY grid over the patch disk, then at each node intersects the vertical line with all patch-atom spheres. It selects:

- outer raw contact from the winning upper intersection
- inner raw contact from the winning lower intersection

The result includes raw scalar fields, inside-disk and valid masks, provenance arrays for winning atoms, contact roles, CSV paths, a contact-atom PDB path, timestamps, and runtime metadata :contentReference[oaicite:47]{index=47}.

Tests show that the selection is based on the winning envelope intersections, not on naive atom-center heuristics, and that provenance is preserved deterministically :contentReference[oaicite:48]{index=48}.

### Stage 5 — surface preparation and trust-domain setup

Stage 5 converts raw Stage 4 evidence into seeds, interpolation-allowed masks, exclusion masks, hard-invalid masks, and a reliable-core mask. It is effectively the pipeline’s trust-domain definition stage. It also derives or accepts parameters like boundary margin, support radius, and reliable radius :contentReference[oaicite:49]{index=49}.

Tests show that this stage explicitly distinguishes seed nodes, interpolation-eligible nodes, boundary-excluded nodes, hard-invalid nodes, and reliable-core nodes, and that provenance from Stage 4 paired seeds is propagated and deduplicated carefully :contentReference[oaicite:50]{index=50}.

### Stage 6 — surface reconstruction

Stage 6 reconstructs continuous outer and inner scalar surfaces over the grid as `z_outer_reconstructed` and `z_inner_reconstructed`. It tracks reconstructed nodes, seed and interpolation counts, final valid analysis mask, non-crossing adjustments, separation statistics, and mesh export metadata :contentReference[oaicite:51]{index=51}.

Tests show that seed values remain fixed, interpolation nodes are solved, hard-invalid nodes remain invalid, non-crossing constraints can be enforced, and deterministic mesh/CSV exports are expected :contentReference[oaicite:52]{index=52}.

### Stage 7 — smoothing / regularization

Stage 7 regularizes Stage 6 surfaces. The config exposes two methods:

- `smooth`
- `thin_plate_grid_fit`

The result includes smoothed outer/inner fields, masks, fit metadata, solver iteration metadata, Stage 7 vs Stage 6 delta summaries, smooth separation summaries, mesh paths, and semantic labels such as normal orientation, thickness-definition label, and metric-surface definition :contentReference[oaicite:53]{index=53}.

Tests demonstrate that Stage 7 supports seed-pinning or non-pinning behavior, reliable-core-restricted metric domains, boundary-condition modes, deterministic fitting, optional delta-map export, and both legacy smoothing and thin-plate grid fitting modes :contentReference[oaicite:54]{index=54}.

### Stage 8 — derivative estimation

Stage 8 estimates first and second derivatives for both outer and inner smoothed surfaces and tracks fit diagnostics and explicit failure reasons. The result includes derivative fields, RMS and absolute residuals, condition indicators, neighbor counts, and a family of masks indicating why a node failed, such as insufficient support, rank deficiency, poor conditioning, high residual, or bad boundary geometry :contentReference[oaicite:55]{index=55}.

Tests show exact recovery on synthetic quadratics, deterministic invalidation on pathological neighborhoods, and explicit export of derivative CSV artifacts :contentReference[oaicite:56]{index=56}.

### Stage 9 — curvature computation and QC

Stage 9 computes mean curvature, oriented mean curvature, Gaussian curvature, graph-normal-versus-radial alignment metrics, orientation-flip flags, and QC layers such as global-tail and local-spike flags plus confidence classes. Summary statistics are stored directly in the result object :contentReference[oaicite:57]{index=57}.

Tests show correct plane and paraboloid behavior, invalid-input and non-finite-output guards, QC flags, and CSV export schema including QC columns :contentReference[oaicite:58]{index=58}.

### Stage 10 — thickness computation

This is the newest and most important correction relative to older documentation. Stage 10 computes thickness using configurable methods:

- `vertical`
- `radial`

Its result stores both vertical and radial thickness-related fields, domain masks, invalid-reason masks, active-method summaries, domain-fraction summaries, semantic metadata, and CSV artifact paths :contentReference[oaicite:59]{index=59}.

The public API and tests reveal several important semantic decisions:

- vertical thickness is still explicitly defined from Stage 7 smooth surfaces as an outer-minus-inner graph-surface separation
- radial thickness is treated as a separate explicit mode, not as a silent reinterpretation
- Stage 10 can optionally restrict its domain to the Stage 9 curvature-valid subset
- radial mode tracks special failure reasons such as no bracket, root failure, and outside-inner-domain loss
- Stage 10 is allowed to run even when Stage 8 failed, because vertical thickness fundamentally depends on Stage 7 surfaces, not on derivatives

All of that is exercised in tests, so it is not speculative :contentReference[oaicite:60]{index=60}.

## 14. Geometry config model

The entire geometry pipeline is driven by `FoldPatchAnalysisConfig`. This config is broad and stage-aware. It includes:

- top-level enable/debug/fold/radius/grid controls
- Stage 5 boundary/support/reliable-domain controls
- Stage 6 reconstruction and mesh export controls
- Stage 7 method and regularization controls
- Stage 8 derivative-fit controls
- Stage 9 curvature QC controls
- Stage 10 thickness controls
- output prefix and rotated-capsid export flag

This matters because the branch favors explicit structured config over hidden constants :contentReference[oaicite:61]{index=61}.

## 15. Geometry result aggregation model

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
- `stage10_thickness`
- `run_summary_json_path`
- `messages`

This is the top-level analysis contract that downstream reporting and future features should consume :contentReference[oaicite:62]{index=62}.

## 16. Output and artifact model

The branch is artifact-rich by design.

### PDB-like artifacts
- final full capsid export
- rotated capsid export
- Stage 2 patch subset export
- Stage 4 contact-atom export

### CSV artifacts
- Stage 3 normalized atom tables
- Stage 4 raw sheet and mask CSVs
- Stage 5 seed/mask summaries
- Stage 6 reconstructed field CSVs
- Stage 7 smoothed fields and delta maps
- Stage 8 derivative CSVs
- Stage 9 curvature CSVs
- Stage 9 curvature-colored PLY meshes (red=negative H, white=0, blue=positive H) when mesh format is `ply`
- Stage 10 thickness CSVs

### Mesh artifacts
- Stage 6 reconstructed meshes
- Stage 7 smoothed meshes
- selectable OBJ or STL export
- split or combined inner/outer output modes

Tests repeatedly enforce existence, determinism, and schema expectations for these artifacts :contentReference[oaicite:63]{index=63}.

## 17. Numerical design character of the branch

The real numerical character of the branch, based on the API and tests, is this:

- it is **grid-first**
- it represents the local shell as two graph surfaces over XY
- it is strongly **fold-aligned** and local, not a global shell manifold framework
- it explicitly separates raw contact evidence, trust-domain construction, reconstruction, smoothing, derivatives, curvature, and thickness
- it relies heavily on explicit masks and deterministic artifact generation for debuggability

This is a practical and inspectable architecture, but it is still graph-surface based rather than a general mesh-native shell geometry framework.

## 18. What is strong in the current design

The strongest parts of the branch are:

- explicit frame/orientation semantics
- centralized symmetry math
- clean separation between persistent structure and derived analysis state
- formal stage boundaries with rich result structs
- strong provenance preservation
- explicit invalid-reason masks instead of undifferentiated failure
- artifact-oriented debugging
- modular tests with synthetic fixtures and integration runs

Those are real strengths visible directly in the current contracts and tests :contentReference[oaicite:64]{index=64} :contentReference[oaicite:65]{index=65}.

## 19. Limitations and technical debt visible from the inspected surface

Several limitations are also clear.

First, the older internal docs lag the code. The public API and tests already include Stage 10, while the stale technical guide still frames Stage 9 as the end of the implemented numerical path :contentReference[oaicite:66]{index=66} :contentReference[oaicite:67]{index=67} :contentReference[oaicite:68]{index=68}.

Second, the build links only one large geometry pipeline source file, `src/geometry_analysis.cpp`, for a very wide set of responsibilities. That strongly suggests implementation concentration in a single compilation unit, which is survivable now and likely to become maintenance debt as the pipeline keeps expanding :contentReference[oaicite:69]{index=69}.

Third, the geometric representation remains a local graph-surface model aligned to a selected fold axis. That is excellent for the intended patch analysis and naturally limiting for more general shell representations.

Fourth, some semantics remain explicitly limited even when honestly labeled. For example, vertical thickness is still graph-surface separation, and radial thickness carries explicit domain-loss limitations in its own API/tests rather than pretending to solve the problem universally :contentReference[oaicite:70]{index=70} :contentReference[oaicite:71]{index=71}.

## 20. Bottom-line technical conclusion

The `geometry` branch should be understood as a **fold-aligned local shell reconstruction and differential-geometry analysis pipeline** on top of a C++17 structural capsid core.

It already includes:

- parser and structural hierarchy
- protein-only structural acceptance policy
- internal subunit reconstruction beyond raw PDB chain semantics
- explicit in-place frame mutation semantics
- canonical icosahedral fold math
- deterministic fold-centered patch extraction
- vdW-normalized analytical patch representation
- raw envelope detection
- trust-domain preparation
- surface reconstruction
- optional smoothing and thin-plate fitting
- derivative estimation
- curvature computation with QC
- thickness computation with vertical and radial modes
- heavy intermediate artifact export
- broad geometry-focused tests

So the honest state is this: this branch is not “early scaffolding” anymore. It is already a substantial scientific-numerical subsystem. The next development work should treat it with that level of seriousness.













