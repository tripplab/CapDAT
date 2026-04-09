# CapDAT Technical Specification and Development Guide

## Purpose

This document consolidates the current CapDAT technical specification and a practical development guide into a single repository-ready reference. It is intended to support onboarding, architectural understanding, and the implementation of new functionality in a way that remains consistent with the current codebase.

The document is based on the present `main` branch implementation of **CapDAT v0.1.0**, including the build system, command-line interface, parser, core domain model, geometry and symmetry utilities, reorientation workflow, export subsystem, structural summary module, and current tests.

---

## 1. Project Overview

**CapDAT (Capsid Data Analysis Toolkit)** is a C++17 command-line application for the structural analysis of viral capsids represented in fixed-column PDB-like coordinate files. The current release is intentionally a foundation release. Its main goal is to establish a trustworthy and extensible structural core rather than a large analytical feature set.

At present, the program can:

- read one input structure file,
- parse `ATOM` records and optionally `HETATM` records,
- retain only the accepted internal structure according to parser policy,
- reconstruct a capsid hierarchy in memory,
- compute a compact extended structural summary,
- optionally apply a post-parse in-place reorientation workflow,
- optionally export the current accepted in-memory capsid structure.

The project currently builds a single executable named `capsid_analyzer` and includes optional test targets when `CAPDAT_BUILD_TESTS` is enabled.

---

## 2. Build and Repository Structure

### 2.1 Build system

CapDAT uses **CMake** as its primary build system and requires **C++17**.

Current top-level build characteristics:

- minimum CMake version: `3.16`
- project version: `0.1.0`
- language: `CXX`
- standard: `C++17`
- optional warnings: enabled by default
- optional tests: disabled by default unless `CAPDAT_BUILD_TESTS=ON`

### 2.2 Main compiled source modules

The current source list reflects the intended architecture:

- `src/main.cpp`
- `src/atom.cpp`
- `src/residue.cpp`
- `src/chain.cpp`
- `src/capsid.cpp`
- `src/pdb_parser.cpp`
- `src/geometry_symmetry.cpp`
- `src/export_capsid.cpp`
- `src/reorientation_workflow.cpp`
- `src/logger.cpp`
- `src/structural_summary.cpp`
- `src/summary_reporter.cpp`
- `src/timer.cpp`

### 2.3 Current tests

When tests are enabled, the repository currently builds dedicated test executables for:

- structural summary,
- geometry and symmetry,
- reorientation workflow,
- export behavior,
- selected CLI-level success and failure cases.

This is already enough to support disciplined incremental development.

---

## 3. Architectural Model

The codebase currently follows a clear layered structure.

### 3.1 Layers

1. **Application shell**  
   `main.cpp` handles CLI parsing, config wiring, high-level orchestration, final reporting, and top-level error handling.

2. **Parsing layer**  
   `PdbParser` is responsible for interpreting fixed-column PDB-like records and building the in-memory structure.

3. **Domain model**  
   `Atom`, `Residue`, `Chain`, and `Capsid` store the accepted structure and its authoritative state.

4. **Geometry and symmetry utilities**  
   `geometry_symmetry` provides reusable post-parse geometric and canonical icosahedral reference logic.

5. **Workflow layer**  
   `reorientation_workflow` resolves user-facing reorientation requests and applies in-place coordinate transforms.

6. **Export layer**  
   `export_capsid` serializes the current accepted in-memory state to a PDB-like file.

7. **Derived analysis and reporting layer**  
   `structural_summary` computes assembly-level geometric descriptors, while `summary_reporter` formats them.

8. **Infrastructure layer**  
   `logger` and `timer` provide lightweight runtime support.

### 3.2 Architectural principle

The most important current architectural principle is this:

- **parsing owns text interpretation,**
- **domain classes own persistent structure and authoritative state,**
- **geometry utilities own reusable post-parse math,**
- **workflows own explicit user-driven mutations,**
- **writers own output serialization,**
- **analysis modules should compute derived values from the accepted structure.**

Future development should preserve that separation whenever possible.

---

## 4. Command-Line Interface Specification

The CLI currently supports the following user-facing options:

### Required

- `-i, --input <file>`  
  Input PDB-like coordinate file.

### Optional

- `-l, --log <file>`  
  Write log output to a file.

- `-v, --verbose`  
  Increase terminal verbosity.

- `--quiet`  
  Reduce terminal output.

- `--include-hetatm`  
  Include `HETATM` records during parsing.

- `--export-final <file>`  
  Export the current accepted capsid coordinates to a PDB-like file.

- `--reorient`  
  Enable the post-parse in-place reorientation workflow.

- `--align-fold <name>`  
  Reorientation source specified as a canonical fold. Supported names are:
  - `2_0`
  - `2_1`
  - `3_0`
  - `3_1`
  - `5_0`

- `--align-vector <x,y,z>`  
  Reorientation source specified as a custom origin-based direction vector.

- `--align-axis <x|y|z>`  
  Target alignment axis. Defaults to `z` if omitted.

- `-h, --help`
- `--version`

### Current CLI invariants

The current implementation explicitly enforces:

- `--input` is required,
- `--align-fold` and `--align-vector` are mutually exclusive,
- either source requires `--reorient`,
- `--reorient` requires exactly one source,
- `--align-axis` accepts only lowercase `x`, `y`, or `z`,
- deprecated `--write-clean-pdb` is rejected in favor of `--export-final`.

This is the expected style for future user-facing features: argument parsing and high-level invariants in `main.cpp`, detailed interpretation and action in a dedicated module.

---

## 5. Core Domain Model

### 5.1 Atom

`Atom` is the smallest persistent structural unit. It stores atom-level metadata parsed from a coordinate record, including:

- serial number,
- atom name,
- alternate-location indicator,
- residue name,
- original PDB chain identifier,
- residue sequence number,
- insertion code,
- Cartesian coordinates,
- occupancy,
- B-factor,
- element,
- charge,
- record type (`ATOM` vs `HETATM`).

Its API is intentionally lightweight. The one notable mutator is `setPosition()`, which exists specifically to support post-parse coordinate workflows such as reorientation.

### 5.2 Residue

`Residue` groups atoms that belong to the same residue identity within one internally reconstructed subunit. It stores:

- residue name,
- residue sequence number,
- insertion code,
- original PDB chain label,
- internal subunit identifier,
- contained atoms.

Mutable atom access is available, but it is intended for explicit post-parse workflows rather than arbitrary topology changes.

### 5.3 Chain

`Chain` is the current implementation name for what is conceptually an **internally reconstructed subunit**.

This distinction is important because large capsid structures may reuse the same one-letter PDB chain label across multiple distinct proteins. `Chain` therefore stores:

- a unique internal subunit identifier,
- the original PDB chain label as metadata,
- an ordered list of residues,
- a cached atom count for efficient summary access.

The current API intentionally narrows mutation through `addAtomToLastResidue()` so the internal atom-count cache remains synchronized.

### 5.4 Capsid

`Capsid` is the top-level domain object. It owns the full assembly and stores:

- source file path,
- all reconstructed internal subunits,
- top-level counters,
- authoritative orientation/frame state.

Its responsibilities in the current release are intentionally limited to storing accepted structure, summary counters, and whole-assembly metadata.

---

## 6. Orientation and Frame State

One of the most important recent architectural features is `Capsid::OrientationState`.

This state records whether the current in-memory coordinates still correspond to the original parsed input frame or to a derived frame produced by an explicit workflow.

It currently stores:

- whether the capsid is still in the original parsed frame,
- whether in-place reorientation occurred,
- whether the request resolved as identity,
- the applied rotation matrix,
- rotation axis and angle when available,
- source mode,
- source description,
- source direction,
- requested target axis,
- target direction.

### Design consequence

Any future feature that depends on frame semantics must consult this state instead of assuming coordinates are still in the original file orientation.

That matters especially for:

- fold-based spatial classification,
- canonical-frame measurements,
- exports,
- future local patch analyses,
- asymmetric-unit calculations.

---

## 7. Parsing Subsystem Specification

### 7.1 Parser role

`PdbParser` is the only module that should understand the fixed-column coordinate format. Its responsibilities are:

- opening the file,
- scanning it line by line,
- recognizing supported coordinate records,
- extracting fields using fixed-column slicing,
- applying validation and policy,
- reconstructing the in-memory hierarchy,
- collecting parsing statistics.

### 7.2 Parser configuration

The parser is configured through `ParserConfig`, which currently includes:

- `include_hetatm`
- `strict_mode`
- `keep_altloc_all`
- `ignore_blank_chain_id`
- `verbose_warnings`
- `protein_only`

Not all of these are fully exploited yet, but they indicate intended extension points.

### 7.3 Parse statistics

`ParseStats` currently tracks:

- total lines read,
- total coordinate records detected,
- total accepted atoms,
- total accepted HETATM records,
- malformed records,
- skipped non-coordinate records,
- warnings,
- alternate-location records,
- internal subunits created,
- residues created.

### 7.4 Current parsing flow

The implementation in `src/pdb_parser.cpp` currently performs this sequence:

1. reset parser state and statistics,
2. open the input file,
3. scan line by line,
4. reject non-coordinate lines,
5. parse candidate coordinate records into a temporary `PdbRecord`,
6. validate the temporary record,
7. apply acceptance policy,
8. append accepted records into the `Capsid` hierarchy in a single pass,
9. finalize hierarchy-derived counts,
10. return the completed `Capsid`.

### 7.5 Temporary `PdbRecord`

The temporary `PdbRecord` type is an important architectural detail. It separates:

- raw text extraction,
- temporary field storage,
- record validation,
- later construction of persistent domain objects.

This separation should be preserved if parser complexity increases or additional coordinate formats are introduced later.

### 7.6 Current validation model

Validation is intentionally lightweight in the current release. The parser currently checks basic record identity, line length, and the presence of key coordinate fields, while numeric parsing mostly uses fallback-based parsing.

This is acceptable for the foundation release, but it is not yet a deeply strict PDB validator.

### 7.7 Current acceptance policy

The present `protein_only` mode uses a residue-name heuristic based on a fixed list of amino-acid residue names and a few broader placeholders such as `UNK`.

This means the parser currently aims to exclude non-protein components such as water, ions, nucleic acids, and ligands, but the mechanism is still heuristic rather than chemically sophisticated.

### 7.8 Current hierarchy reconstruction rule

The current internal subunit reconstruction logic is simple by design:

- the first accepted record starts the first internal subunit,
- a change in the raw PDB chain label starts a new internal subunit,
- a change in residue sequence number or insertion code starts a new residue.

This is one of the clearest known limitations in the current codebase and one of the most important likely targets for future refinement.

---

## 8. Geometry and Symmetry Module

`geometry_symmetry` is a standalone post-parse utility namespace for canonical icosahedral geometry.

### 8.1 Current responsibilities

It currently provides:

- canonical fold definitions,
- vector and matrix types,
- vector norms and normalization,
- dot and cross products,
- angular comparisons,
- fold-axis angle and similarity helpers,
- nearest-fold assignment among candidate folds,
- Cartesian and point-to-axis distance helpers,
- proper rotation construction,
- rotation application to directions and points,
- canonical IAU definition,
- IAU half-space classification helpers,
- proper-rotation matrix checks.

### 8.2 Canonical folds

The current fold registry includes five canonical folds:

- `2_0`
- `2_1`
- `3_0`
- `3_1`
- `5_0`

Their exact reference vectors are hardcoded and treated as the single source of truth for the current VIPERdb-like canonical frame. The implementation also enforces the invariant that `2_0` must align with `+Z`.

### 8.3 Development significance

This module should be the preferred place for future reusable geometric logic, including but not limited to:

- fold-relative local coordinate systems,
- patch direction classification,
- asymmetric-unit membership,
- angle-based fold assignment,
- canonical-frame orientation utilities,
- future symmetry-aware reporting.

New geometry logic should generally be added here rather than reimplemented in workflows or analytical modules.

---

## 9. Reorientation Workflow Specification

The current codebase includes a dedicated workflow layer in `reorientation_workflow`.

### 9.1 Request model

`ReorientationRequest` carries the user request in a structured way, including:

- whether reorientation is requested,
- source mode,
- fold name or custom vector text,
- target axis,
- export request flags,
- output path,
- verbosity.

### 9.2 Result model

`ReorientationResult` reports:

- status,
- resolved source description,
- source direction,
- target axis,
- target direction,
- rotation matrix,
- rotation axis,
- rotation angle,
- whether coordinates changed,
- messages for logging or inspection.

### 9.3 Current workflow behavior

The workflow currently:

1. resolves the source direction from either a canonical fold or a custom vector,
2. resolves the target axis to `+X`, `+Y`, or `+Z`,
3. computes a pure proper rotation about the origin,
4. applies the rotation in place to all atom coordinates if needed,
5. writes authoritative orientation state back into the `Capsid`,
6. returns a structured result.

### 9.4 Design implication

This module is the template for future user-driven workflow features. Anything that is multi-step, requires validation, and may mutate the current structure should probably follow this pattern rather than expanding `main.cpp` directly.

---

## 10. Export Subsystem Specification

The export layer writes the **current accepted in-memory structure**, not a byte-faithful copy of the original input file.

That distinction is central to the current design.

### 10.1 Export semantics

The writer can currently:

- emit header remarks,
- emit `TER` records,
- emit an `END` record,
- preserve original atom serial numbers or renumber them,
- write either coordinate-only output or richer PDB-like output.

### 10.2 Orientation-aware export

If the current coordinates have been reoriented in place, the export remarks document that fact and include the corresponding alignment metadata.

### 10.3 Design implication

The export layer should remain the canonical mechanism for materializing the current accepted structural state. Future writers, even for other formats, should follow the same principle:

- serialize what the `Capsid` currently represents,
- do not silently assume original-frame coordinates,
- record transform metadata when relevant.

---

## 11. Structural Summary Module

`structural_summary` is currently the principal example of a post-parse analytical module.

### 11.1 Current computed metrics

It currently computes assembly-level descriptors from accepted atoms only:

- accepted atom count,
- geometric center,
- coordinate min/max bounds,
- axis spans,
- radial min/max/mean/stddev relative to the geometric center,
- shell-thickness estimate defined as `r_max - r_min`,
- bounding-box center,
- per-subunit atom count min/max/mean,
- per-subunit residue count min/max/mean,
- unique original label count,
- sorted set of unique original labels.

### 11.2 Reporting separation

The derived values are stored in a plain `StructuralSummary` struct and formatted separately by `summary_reporter`.

### 11.3 Development significance

This is a useful pattern for future analytical features:

- one read-only input object,
- one plain result struct,
- one compute function,
- one optional formatting layer.

That pattern should be reused whenever possible.

---

## 12. Logging and Runtime Infrastructure

The logger is intentionally lightweight and concrete.

### 12.1 Log levels

Current supported levels are:

- `ERROR`
- `WARNING`
- `INFO`
- `DEBUG`

### 12.2 Current behavior

The logger supports:

- verbosity threshold control,
- terminal messages,
- optional log-file output,
- lightweight timestamped formatting.

### 12.3 Development guidance

Future modules should use the logger rather than writing ad hoc terminal output, except for deliberate final summaries already centralized in reporter functions or `main.cpp`.

---

## 13. Current Tests and Their Implications

The present tests are not just correctness checks; they also reveal intended usage patterns.

### 13.1 Structural summary tests demonstrate

- how to build tiny in-memory capsids for deterministic checks,
- how to validate whole-assembly metrics,
- how to test parser integration with temporary PDB-like files,
- how parser policy affects accepted geometry.

### 13.2 Geometry and symmetry tests demonstrate

- canonical fold registry expectations,
- deterministic lookup behavior,
- proper rotation properties,
- IAU classification expectations.

### 13.3 Reorientation tests demonstrate

- in-place coordinate changes,
- persisted orientation state,
- custom-vector handling,
- rejection of invalid zero-vector requests.

### 13.4 Export tests demonstrate

- original-frame vs transformed-frame remarks,
- orientation-aware output semantics,
- export of current in-memory coordinates.

### 13.5 Development implication

New features should follow the same testing style:

- start with small deterministic constructed objects,
- then add narrow integration tests if needed,
- keep test targets independent and purpose-specific.

---

## 14. Known Technical Limitations

The codebase is already well-structured, but several limitations are currently explicit and important.

### 14.1 Internal subunit reconstruction remains simplistic

The current rule based on contiguous changes in raw chain label is not yet sufficient to fully solve repeated chain-label reuse in large capsids.

### 14.2 Protein-only filtering is heuristic

Residue-name matching is useful as a foundation policy, but it is not yet a comprehensive treatment of non-standard or modified residues.

### 14.3 Validation is intentionally lightweight

The parser is not yet a strict format validator and still relies on fallback-based parsing in several places.

### 14.4 Feature set is still intentionally narrow

The current analysis layer stops at assembly-level structural summary. There is not yet a general local geometry, surface, curvature, thickness, anisotropy, or batch-processing subsystem.

### 14.5 Single-file workflow

The current application flow is centered on processing one input file per run.

---

## 15. Development Guide

This section translates the current architecture into practical implementation guidance.

### 15.1 First question: what kind of feature is this?

Before writing code, classify the feature.

#### Add to the parser when

- you need new fields from the input file,
- you must change acceptance policy,
- you must improve hierarchy reconstruction,
- the persistent structural meaning of the stored object changes.

#### Add to the domain model when

- the information must remain part of the authoritative in-memory structure,
- multiple later modules must read it,
- it represents actual persistent state rather than a temporary analysis.

#### Add a new analysis module when

- the feature computes derived measurements,
- it should not mutate topology or persistent structural identity,
- it can consume `const Capsid&` and return a result struct.

#### Add a workflow module when

- the feature is user-triggered,
- it requires option resolution or validation,
- it may mutate coordinates or authoritative state,
- it benefits from explicit request/result models.

#### Add a writer/export module when

- the feature serializes current accepted state,
- output formatting deserves separation from analysis code,
- frame-aware metadata may need to be emitted.

### 15.2 Preferred implementation pattern for a new analysis module

A good new analysis feature should usually follow this pattern:

1. add a new header in `include/`,
2. define one or more plain result structs,
3. expose a read-only compute function,
4. implement the logic in `src/`,
5. optionally create a dedicated reporter/formatter,
6. add the source to `CMakeLists.txt`,
7. add focused tests in `tests/`,
8. call the feature from `main.cpp` only after parsing and any requested workflows are complete.

#### Recommended shape

- `include/<feature>.hpp`
- `src/<feature>.cpp`
- optionally `include/<feature>_reporter.hpp`
- optionally `src/<feature>_reporter.cpp`
- `tests/<feature>_tests.cpp`

### 15.3 Preferred implementation pattern for a new workflow module

A workflow should follow the reorientation model:

- one request struct,
- one status enum,
- one result struct,
- one top-level application function,
- explicit logging where useful,
- state changes applied only after successful resolution.

### 15.4 Keep frame semantics explicit

If a future feature depends on the canonical fold frame, original parsed orientation, or user-requested alignment, it must consult `Capsid::orientationState()`.

Do not silently assume the structure is still in the original input frame.

### 15.5 Prefer read-only traversal unless mutation is essential

For most analytical features, use the standard traversal:

- iterate through `capsid.chains()`
- iterate through `chain.residues()`
- iterate through `residue.atoms()`

Use mutable access only for explicit workflow operations.

### 15.6 Centralize reusable geometry

If a new feature needs:

- fold vectors,
- fold-axis angles,
- nearest-fold assignment,
- IAU classification,
- proper rotations,
- canonical-frame directional comparisons,

it should use `geometry_symmetry` rather than redefining those calculations elsewhere.

### 15.7 Keep parsing concerns inside the parser

New scientific features should not parse raw coordinate lines directly.

If extra information must be read from the file, the parser should be extended cleanly through:

- `PdbRecord` field extraction,
- `validateRecord()`,
- `shouldAcceptRecord()`,
- `appendRecordToCapsid()`.

### 15.8 Use tests as design tools

The current test style is especially useful and should be preserved:

- small deterministic in-memory capsids for exact expectations,
- tiny temporary PDB-like files for parser integration,
- focused test executables rather than one giant test binary.

---

## 16. Recommended Near-Term Development Priorities

Based on the current codebase, the most natural next functionality families are the following.

### 16.1 Stronger internal subunit reconstruction

This is likely the highest-priority structural-core improvement because many downstream geometric analyses will depend on correct subunit identity.

### 16.2 Fold-centered local analysis

Examples:

- atom or residue classification by nearest fold,
- angular patch extraction around canonical folds,
- IAU-aware partitioning.

### 16.3 Radial shell profiling

A natural extension of the current structural summary would include:

- shell distributions,
- per-subunit radial ranges,
- inner/outer shell statistics,
- radial occupancy summaries.

### 16.4 Local geometry descriptors

Examples:

- local thickness estimates,
- local curvature scaffolding,
- patch-level surface descriptors,
- orientation-dependent local metrics.

### 16.5 New machine-readable exports

Examples:

- CSV metric tables,
- JSON summaries,
- local-frame export bundles,
- future mesh or field output once surface modules exist.

### 16.6 Batch-processing support

Examples:

- processing multiple files per run,
- standardized output directory layouts,
- machine-readable summaries for later pipelines.

---

## 17. Recommended Contribution Style

To remain consistent with the repository’s current direction, new code should generally favor:

- explicit small APIs,
- plain result/config structs,
- focused `.cpp` implementations,
- narrow mutation points,
- deterministic tests,
- minimal hidden behavior,
- clear separation of parsing, analysis, workflows, and output.

That style is already visible across the parser, geometry utilities, reorientation workflow, export layer, and structural summary module.

---

## 18. Minimal Checklist Before Implementing a New Feature

Before starting implementation, check the following:

1. Is this parser behavior, domain state, analysis, workflow, or export?
2. Does it depend on original vs reoriented frame?
3. Should it mutate the `Capsid`, or only read it?
4. Can it reuse `geometry_symmetry` utilities?
5. Does it need a config/request struct?
6. Does it need a result struct?
7. Where will the tests live?
8. Does `CMakeLists.txt` need a new source or test target?
9. Does the CLI need new user-facing options?
10. Does output need to record transform or provenance metadata?

---

## 19. Summary

CapDAT already has a strong architectural spine for future work.

Its current design is especially promising because:

- parsing is isolated,
- structure ownership is explicit,
- geometry logic is centralized,
- user-triggered mutation is handled through workflows,
- export is authoritative and frame-aware,
- tests already encode intended behavior.

For future development, the safest and most scalable strategy is to extend those layers rather than collapse them together.

In practice, that means:

- keep raw text interpretation in the parser,
- keep structural truth in the domain model,
- keep reusable geometry in `geometry_symmetry`,
- keep multi-step user actions in dedicated workflows,
- keep derived measurements in separate analysis modules,
- keep serialization in dedicated writers,
- keep tests small, deterministic, and feature-specific.

That approach should allow CapDAT to grow into more advanced capsid analysis functionality without losing the clarity established by the current v0.1 foundation release.
