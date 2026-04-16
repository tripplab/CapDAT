# CapDAT — Synthesized Development Guide for Implementing New Functionality

## 1. Start by putting the feature in the correct layer

Most extension mistakes in this repo will come from adding logic to the wrong layer.

### Put it in the parser only if
- the feature depends on raw text extraction from coordinate records
- acceptance/rejection policy changes
- hierarchical reconstruction rules change
- the new state belongs to the persistent structural truth

### Put it in `Atom`, `Residue`, `Chain`, or `Capsid` only if
- the data is authoritative and long-lived
- multiple workflows need it as structural truth
- it is not just a derived analysis artifact

### Put it in `geometry_symmetry` if
- the logic is reusable icosahedral geometry
- it involves fold definitions, canonical axes, vector alignment, IAU logic, or reusable rotation math
- more than one future geometry feature could need it

### Put it in `reorientation_workflow` if
- the feature changes the coordinate frame
- it is a user-facing transform workflow
- it must update `Capsid::OrientationState`
- it belongs to source/axis resolution and transform semantics

### Put it in geometry analysis or a sibling geometry module if
- it works on accepted coordinates after parsing
- it consumes Stage 1+ data or later stage outputs
- it is patch logic, surface logic, metrics, derivatives, curvature, thickness, QC, or reporting

### Put it in export code only if
- the logic is pure serialization
- the feature writes already-computed state
- it should not decide or alter the scientific algorithm

If you skip this classification step, you will create design drift fast.

## 2. Respect the current execution model

The repo already has a real workflow. Reuse it instead of improvising a parallel one.

### Actual high-level flow
1. parse input into `Capsid`
2. print structural summary
3. optionally run standalone reorientation
4. or run geometry analysis:
   1. Stage 1 prepare and align fold to `+Z`
   2. Stage 2 select cylindrical patch
   3. Stage 3 normalize analytical patch
   4. Stage 4 detect raw envelope sheets
   5. Stage 5 prepare seeds and trust domain
   6. Stage 6 reconstruct surfaces
   7. Stage 7 smooth / regularize
   8. Stage 8 estimate derivatives
   9. Stage 9 compute curvature and QC
   10. Stage 10 compute thickness
5. optionally export final full capsid state

### Non-negotiable rule
If your new feature depends on fold-centered local geometry, do **not** recompute preprocessing from scratch. Consume the appropriate stage results.

## 3. Always make frame assumptions explicit

This branch already solved a major class of bugs by making orientation state explicit.

### Always know which frame your feature assumes
- original parsed frame
- current arbitrary in-memory frame
- Stage 1 fold-aligned working frame
- later stage analysis frame derived from Stage 1

### The authoritative frame source
Use `capsid.orientationState()`.

### Rule
If your feature assumes the selected fold is aligned to `+Z`, require that precondition explicitly. Do not silently assume it.

A clean geometry feature should usually consume:
- `GeometryPreparationResult`, or
- a later stage result whose existence already guarantees Stage 1 success

## 4. Reuse the existing result contracts

The branch already gives you the core data boundaries.

### Existing contracts you should prefer to reuse
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
- `GeometryStage10ThicknessResult`
- `GeometryAnalysisResult`

### Do not do this
Do not create random duplicate containers for:
- patch atoms
- masks
- derivative maps
- curvature fields
- thickness fields
- reliable-domain semantics
- provenance data

unless the current result model is genuinely insufficient.

## 5. Best extension points by feature type

### 5.1 New scientific summary metrics
The cleanest place is **after Stage 10** now, not after Stage 9. The code already computes surfaces, derivatives, curvature, and thickness. A new scientific summary module should usually consume:

- `GeometryStage5SurfacePrepResult`
- `GeometryStage7SmoothedSurfaceResult`
- `GeometryStage8DerivativeEstimationResult`
- `GeometryStage9CurvatureComputationResult`
- `GeometryStage10ThicknessResult`

This is the highest-value current extension point because the branch now has enough numerical material to synthesize meaningful scientific outputs.

### 5.2 New smoothing or regularization methods
Add them to Stage 7 only if they still conceptually produce the same kind of outer/inner smoothed scalar fields over the same grid domain.

If the new method changes representation too radically, do not force it into Stage 7. Create a sibling module.

### 5.3 New derivative logic
If the feature changes support geometry, fit models, residual filters, or derivative-confidence logic, Stage 8 is the correct insertion point.

### 5.4 New curvature logic
If the feature changes sign conventions, outward-normal treatment, curvature formulas, or QC rules, Stage 9 is the correct insertion point.

### 5.5 New thickness logic
If the feature changes how shell thickness is computed, Stage 10 is now the formal location. Do not bolt thickness alternatives into Stage 7 metadata or ad hoc report code.

## 6. Use `geometry_symmetry` aggressively

Whenever a feature needs:
- fold lookup
- fold vectors
- canonical directions
- angle to an axis
- point-to-axis distance
- direction alignment
- IAU membership logic
- proper-rotation validation

use `geometry_symmetry`.

If you are hardcoding fold vectors or duplicating alignment math elsewhere, you are doing it wrong.

## 7. Preserve persistent/derived separation

The current architecture gets one thing very right: persistent structural truth and derived analytical state are separate.

### Persistent data belongs in
- `Atom`
- `Residue`
- `Chain`
- `Capsid`

### Derived data belongs in
- stage result structs
- analytical patch structs
- metrics or reporting structs
- export-only artifact models

Do not contaminate `Atom` or `Capsid` with stage-local numeric fields just because it is convenient.

## 8. Suggested file/module shape for new work

### For a new read-only analysis feature
Use:
- `include/<feature>.hpp`
- `src/<feature>.cpp`
- `tests/<feature>_tests.cpp`

### Good examples
- `include/geometry_metrics.hpp`
- `src/geometry_metrics.cpp`
- `tests/geometry_metrics_tests.cpp`

### For reusable helpers
Extract focused helpers such as:
- `geometry_grid_utils.*`
- `geometry_mask_utils.*`
- `geometry_surface_io.*`
- `geometry_mesh_export.*`
- `geometry_derivative_utils.*`
- `geometry_curvature_utils.*`
- `geometry_thickness_utils.*`

### Bad habit to avoid
Do not keep dumping unrelated logic into `geometry_analysis.cpp` forever. The build surface already suggests that this file is the central geometry compilation unit, and continued growth there will make maintenance worse, not better.

## 9. Coding style the branch already prefers

Match the repo’s style instead of fighting it.

### The visible branch style is
- plain config structs
- plain result structs
- explicit stage entry functions
- explicit validation and hard precondition checks
- deterministic behavior
- logger-based messages
- artifact paths carried in results
- provenance preserved where possible
- tests that assert schemas, masks, counts, and determinism

### Good addition pattern
- `FeatureConfig`
- `FeatureResult`
- `runFeature(...)`

### Bad addition pattern
- hidden stateful classes with unclear lifecycle
- numerical logic embedded in CLI parsing
- silent mutation of persistent objects
- export code deciding scientific semantics
- global behavior hidden behind conditionals in `main.cpp`

This codebase is explicit, not “clever”.

## 10. How to wire new CLI-controlled functionality

Follow the existing `main.cpp` discipline.

### Correct sequence
1. add local CLI variables near related controls
2. parse the new option in the argument loop
3. keep only high-level mode and numeric sanity checks in `main.cpp`
4. map values into config structs
5. call a dedicated module function after prerequisites are ready
6. keep actual numerical work outside `main.cpp`

### Wrong sequence
- writing formulas directly inside `main.cpp`
- branching around the stage pipeline with one-off shortcuts
- bypassing config structs and smuggling hidden constants

That is how repos become unmaintainable.

## 11. How to design a new metrics/reporting module

This is probably the best next move for this branch.

### Recommended inputs
A post-analysis module should usually take only the stage results it actually needs, for example:

- `const GeometryStage5SurfacePrepResult&`
- `const GeometryStage7SmoothedSurfaceResult&`
- `const GeometryStage8DerivativeEstimationResult&`
- `const GeometryStage9CurvatureComputationResult&`
- `const GeometryStage10ThicknessResult&`
- `const FoldPatchAnalysisConfig&`
- `Logger*`

### Recommended outputs
Return a plain result struct with:
- success flag
- scalar summaries
- optional nodewise maps
- QC summaries
- artifact paths
- messages

### Recommended file names
- `include/geometry_metrics.hpp`
- `src/geometry_metrics.cpp`
- `tests/geometry_metrics_tests.cpp`

### Recommended content
A first useful metrics layer could compute:
- patch-level thickness summaries from Stage 10 active method
- curvature-valid coverage
- derivative-valid coverage
- radial-thickness sampleability loss burden
- Stage 7 non-crossing-adjustment burden
- Stage 9 QC warning burden
- integrated quality summary for “trust this patch / distrust this patch”

That adds scientific value without destabilizing the current stage pipeline.

## 12. Guidance for likely next feature categories

### 12.1 Better final reporting
This is the clearest current gap. The branch computes a lot, but users still need a cleaner final scientific summary layer.

### 12.2 Stronger thickness analysis
Stage 10 already exists, so new thickness work should extend it, not reinvent it somewhere else. Good extensions include:
- clearer vertical-versus-radial comparison summaries
- coverage and failure-mode summaries
- confidence weighting by Stage 9 or Stage 8 validity burdens
- future normal-based thickness methods if the representation can support them honestly

### 12.3 Area metrics
A dedicated metrics layer can cleanly add:
- projected patch area
- reconstructed outer-surface area
- reconstructed inner-surface area
- active metric-domain area
- curvature-valid subdomain area

Do not conflate those. They are not the same quantity.

### 12.4 QC synthesis
The branch already has raw ingredients for a good trust score:
- reliable-core fraction
- reconstructed fraction
- smooth-valid fraction
- derivative-valid fraction
- curvature-valid fraction
- Stage 7 adjustment burden
- Stage 9 QC warnings
- Stage 10 radial domain loss

A synthesis module should summarize those instead of making users inspect many CSVs manually.

## 13. How to add a new Stage 7 method cleanly

If the new method still belongs to Stage 7:

1. extend `FoldPatchAnalysisConfig::Stage7Method`
2. add CLI parsing in `main.cpp`
3. preserve the existing Stage 7 result contract
4. implement the method in a dedicated internal function
5. wire dispatch inside `runGeometryAnalysisStage7SurfaceSmoothing(...)`
6. add focused unit tests and integration tests

Do not invent a method-specific output shape if it still belongs to Stage 7. Keep the stage contract stable.

## 14. How to add a new Stage 8 or Stage 9 method cleanly

### For Stage 8
Preserve:
- explicit failure masks
- per-node diagnostics
- derivative-valid semantics
- deterministic rejection behavior

### For Stage 9
Preserve:
- outward/inward normal semantics
- orientation-flip tracking
- QC layers
- summary statistics
- explicit CSV schema expectations

Do not silently change sign conventions or node-validity semantics.

## 15. How to add or modify thickness functionality cleanly

Because Stage 10 already formalizes thickness, new thickness work should follow this path:

1. decide whether the new method is truly a Stage 10 method
2. extend `Stage10ThicknessMethod` if appropriate
3. add CLI mapping
4. implement a dedicated method branch in Stage 10
5. preserve the existing result object fields where possible
6. add new fields only when materially necessary
7. add tests for:
   - valid cases
   - invalid domain cases
   - threshold failures
   - method metadata
   - summary field semantics
   - CSV schema

Do not hide a new thickness method behind reused labels. Be explicit.

## 16. Testing strategy the branch expects

The current repo already tells you how to test new code.

### Layer 1 — small deterministic unit tests
Use these for:
- algebra
- masks
- formula correctness
- edge conditions
- sign conventions
- threshold logic

### Layer 2 — synthetic stage fixtures
Use these for:
- stage contracts
- provenance rules
- node-validity logic
- failure-reason logic
- exact recovery on synthetic fields
- pathological counterexamples

### Layer 3 — integration tests
Use these for:
- multi-stage workflows
- artifact generation
- CLI wiring
- deterministic output
- realistic data sanity

### What not to do
Do not rely only on one big real-capsid run. That is slow, opaque, and a terrible debugging tool.

## 17. When you should refactor before extending

Be honest about this. Sometimes the right move is not “add the feature now”.

### Refactor first if
- the new feature touches helper logic that is already repeated across stages
- you are adding another big unrelated block into `geometry_analysis.cpp`
- you cannot test the new feature without running a huge multi-stage orchestration
- Stage 7, 8, 9, or 10 branching is becoming crowded
- export code and numerical code are getting intertwined

### Practical early refactor targets
- grid utilities
- mask and domain helpers
- CSV artifact writers
- mesh export helpers
- derivative fit helpers
- curvature helpers
- thickness ray/intersection helpers

## 18. A sensible next-step roadmap

If the goal is maximum value with minimum architectural damage, the best development sequence is:

1. build a dedicated post-Stage-10 metrics/reporting module
2. add integrated QC and trust summaries
3. add area metrics and active-domain coverage metrics
4. improve thickness comparison/reporting between vertical and radial methods
5. only then consider representation-level upgrades or heavier solver changes

That order matters. It exposes value sooner and prevents getting stuck in solver tinkering before the branch has a coherent scientific reporting layer.

## 19. Minimal checklist before writing code

Use this every time:

1. Which layer does this feature belong to?
2. What frame does it assume?
3. Can it consume an existing stage result?
4. Should it reuse `geometry_symmetry`?
5. Does it need a dedicated config struct?
6. Does it need a dedicated result struct?
7. Does it need export artifacts?
8. Does it need CLI wiring?
9. Which exact tests prove it works?
10. Is `geometry_analysis.cpp` really the right home, or should this be a new module?

If you cannot answer those clearly, the implementation is not ready.

## 20. Bottom-line implementation advice

### Do this
- build on existing stage outputs
- keep frame semantics explicit
- preserve provenance
- use plain config/result structs
- keep numerical code out of `main.cpp`
- keep export code serialization-only
- add deterministic tests
- extract helper modules as soon as logic repeats

### Do not do this
- duplicate fold math
- assume `+Z` alignment silently
- reparse raw coordinate text inside geometry features
- store stage-local numeric state in `Capsid`
- mix scientific formulas with artifact formatting
- keep growing one giant geometry source file without restraint

### Best immediate next move
Implement a **read-only post-Stage-10 metrics and reporting module**. That matches the actual current architecture, uses the full Stage 1–10 pipeline correctly, and closes the biggest remaining gap between numerical computation and user-facing scientific output.
