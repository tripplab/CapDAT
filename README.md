# Capsid Data Analysis Toolkit (CapDAT)
High-performance software tool for the structural analysis of viral capsids from atomic coordinate data in Protein Data Bank (PDB) format

by trippm@tripplab.com [04/2026]

geometry branch

**CapDAT (Capsid Data Analysis Toolkit)** is a C++17 command-line project for parsing, filtering, organizing, reorienting, exporting, and geometrically analyzing large icosahedral capsid structures from fixed-column PDB-like coordinate files, including VIPERdb-style `.vdb` inputs. The current `geometry` branch includes the structural core, a centralized icosahedral symmetry module, an in-place reorientation workflow, canonical final export, and a staged fold-centered geometry-analysis pipeline under active development.

CapDAT is designed around capsid-specific realities that standard PDB tooling handles poorly. In particular, large capsids may contain hundreds of protein subunits, which means the original one-letter chain identifier cannot be treated as a globally unique molecular key. CapDAT therefore preserves the raw chain label only as metadata while reconstructing its own internal hierarchy as **Atom -> Residue -> internal subunit (`Chain`) -> Capsid**. This makes the internal model suitable for large quaternary assemblies and for downstream geometry workflows that must remain deterministic and traceable.

At parse time, CapDAT reads coordinate records, extracts relevant structural fields, and filters the input so that only the atoms belonging to the capsid protein are retained internally by default. Non-protein entities such as nucleic acids, waters, ions, ligands, and other unrelated components are excluded from the analysis basis unless explicitly allowed by parser/CLI policy. The program also tracks malformed records, skipped records, alternate locations, accepted hetero atoms, and reconstructed subunits, and reports an extended structural summary after parsing.

The `geometry` branch also includes a dedicated `geometry_symmetry` module that centralizes canonical icosahedral fold definitions, deterministic fold lookup, reusable angle/distance utilities, and direction-alignment rotation math. On top of that, CapDAT provides an opt-in post-parse in-place reorientation workflow and a fold-centered geometry-analysis workflow driven from the CLI. The geometry-analysis path currently exposes a multi-stage pipeline with configuration for fold selection, cylindrical patch extraction, vdW-radius normalization, grid generation, surface reconstruction support, boundary handling, mesh export, and smoothing / thin-plate-style fitting controls.

The executable currently supports:
- input-file selection and standard help/version output,
- optional log-file output and verbosity control,
- optional inclusion of `HETATM` records,
- optional final export of the current accepted coordinates through `--export-final`,
- optional post-parse reorientation through `--reorient`,
- optional geometry-analysis execution through `--geometry-analysis`,
- geometry-stage tuning parameters for patch selection, surface reconstruction, and smoothing.

The current code base on this branch builds:
- the main executable `capsid_analyzer`,
- structural summary tests,
- geometry/symmetry tests,
- reorientation workflow tests,
- geometry-analysis tests,
- export tests,
- CLI integration tests registered through CTest.

## Tech Requirements

CapDAT currently requires the following development stack:

- **C++ standard:** C++17
- **Build system:** CMake (minimum 3.16)
- **Recommended generator:** Ninja
- **Test framework:** CTest through CMake test registration
- **Supported toolchain:** a modern C++ compiler compatible with C++17
- **Optional environment management:** micromamba

Current build options exposed by the project include:

- `CAPDAT_BUILD_TESTS=ON|OFF` to enable or disable test targets
- `CAPDAT_ENABLE_WARNINGS=ON|OFF` to enable compiler warnings
- `CMAKE_BUILD_TYPE=Debug|Release|...` for standard CMake build-type selection

The current project configuration explicitly sets:
- `CMAKE_CXX_STANDARD 17`
- `CMAKE_CXX_STANDARD_REQUIRED ON`
- `CMAKE_CXX_EXTENSIONS OFF`

## How to Clone CapDAT from GitHub

To clone the **`geometry` branch specifically**, use one of the following commands.

For **SSH**:

`git clone --branch geometry --single-branch git@github.com:tripplab/CapDAT.git`

For **HTTPS**:

`git clone --branch geometry --single-branch https://github.com/tripplab/CapDAT.git`

Then enter the repository:

`cd CapDAT`

Verify that you are really on the correct branch:

`git branch --show-current`

The command should print:

`geometry`

If you already cloned the repository previously and only need to switch to the branch, use:

`git fetch origin`
`git switch geometry`

or, if the branch does not yet exist locally:

`git fetch origin geometry`
`git switch --track origin/geometry`

## How to Define an Environment for CapDAT

If your machine already has a working C++17 compiler, CMake, and Ninja, you do **not** need a dedicated environment. If you want a clean and reproducible toolchain, define one explicitly.

A practical micromamba environment is:

`micromamba create -n capdat -c conda-forge python=3.11 cmake ninja cxx-compiler -y`

Activate it with:

`micromamba activate capdat`

If micromamba shell integration is not initialized yet, run:

`eval "$(micromamba shell hook --shell bash)"`
`micromamba activate capdat`

Then verify the toolchain:

`cmake --version`
`ninja --version`
`c++ --version`

The environment is sufficient for configuring, building, and testing the current `geometry` branch.

## How to Build CapDAT

From the repository root, configure the project in **Debug** mode and explicitly build the test targets:

`cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCAPDAT_BUILD_TESTS=ON`

Then compile:

`cmake --build build`

This should generate the main executable:

`build/capsid_analyzer`

It should also build the current test executables registered by the project, including the structural-summary, geometry-symmetry, reorientation-workflow, geometry-analysis, and export test binaries.

If you want to confirm the configured project summary immediately after configuration, look for output of the form:

- `Configuring CapDAT v0.1.0`
- `Build type: Debug`
- `C++ standard: 17`
- `Build tests: ON`
- `Warnings enabled: ON`

To confirm that the executable was built correctly, run:

`./build/capsid_analyzer --help`

## How to Test the CapDAT Build

Once the build completes, run the full registered test suite with CTest:

`ctest --test-dir build --output-on-failure`

This is the correct test entry point for the current branch because the project registers both unit-style test binaries and CLI integration tests directly in CMake.

A successful test pass should include coverage for:
- `capdat_structural_summary_tests`
- `capdat_geometry_symmetry_tests`
- `capdat_reorientation_workflow_tests`
- `capdat_geometry_analysis_tests`
- `capdat_export_capsid_tests`
- CLI export/reorientation checks such as:
  - `capdat_cli_export_only`
  - `capdat_cli_fold_to_z_export`
  - `capdat_cli_fold_to_x_export`
  - `capdat_cli_vector_to_y_export`
  - `capdat_cli_invalid_combo`

You should also do a direct executable sanity check:

`./build/capsid_analyzer --version`
`./build/capsid_analyzer --help`

And a real sample-data run, for example:

`./build/capsid_analyzer --input data/1cwp_full.vdb --verbose --log run.log`

For a geometry workflow sanity check, for example:

`./build/capsid_analyzer --input data/1cwp_full.vdb --geometry-analysis --debug --geometry_out_prefix geometry`

These commands validate different parts of the current development branch:
- basic CLI wiring,
- parser execution,
- structural summary reporting,
- optional logging,
- geometry pipeline entry,
- debug-artifact generation path.

## Example Output from the CapDAT Build Test

A **successful configuration** should produce output similar to:

`-- Configuring CapDAT v0.1.0`
`-- Build type: Debug`
`-- C++ standard: 17`
`-- Build tests: ON`
`-- Warnings enabled: ON`
`-- Configuring done`
`-- Generating done`
`-- Build files have been written to: /path/to/CapDAT/build`

A **successful build** should end with output similar to:

`[... ] Linking CXX executable capsid_analyzer`
`[... ] Linking CXX executable capdat_structural_summary_tests`
`[... ] Linking CXX executable capdat_geometry_symmetry_tests`
`[... ] Linking CXX executable capdat_reorientation_workflow_tests`
`[... ] Linking CXX executable capdat_geometry_analysis_tests`
`[... ] Linking CXX executable capdat_export_capsid_tests`

A **successful CTest run** should look similar to:

`Internal ctest changing into directory: /path/to/CapDAT/build`
`Test project /path/to/CapDAT/build`
`    Start 1: capdat_structural_summary_tests`
`1/10 Test #1: capdat_structural_summary_tests ...........   Passed`
`    Start 2: capdat_geometry_symmetry_tests`
`2/10 Test #2: capdat_geometry_symmetry_tests ............   Passed`
`    Start 3: capdat_reorientation_workflow_tests`
`3/10 Test #3: capdat_reorientation_workflow_tests .......   Passed`
`    Start 4: capdat_geometry_analysis_tests`
`4/10 Test #4: capdat_geometry_analysis_tests ............   Passed`
`    Start 5: capdat_export_capsid_tests`
`5/10 Test #5: capdat_export_capsid_tests ................   Passed`
`    Start 6: capdat_cli_export_only`
`6/10 Test #6: capdat_cli_export_only ....................   Passed`
`    Start 7: capdat_cli_fold_to_z_export`
`7/10 Test #7: capdat_cli_fold_to_z_export ...............   Passed`
`    Start 8: capdat_cli_fold_to_x_export`
`8/10 Test #8: capdat_cli_fold_to_x_export ...............   Passed`
`    Start 9: capdat_cli_vector_to_y_export`
`9/10 Test #9: capdat_cli_vector_to_y_export .............   Passed`
`    Start 10: capdat_cli_invalid_combo`
`10/10 Test #10: capdat_cli_invalid_combo ................ Passed`
`100% tests passed, 0 tests failed out of 10`

A **successful sample parser run** should look similar to:

`[YYYY-MM-DD HH:MM:SS]`
`[INFO] Starting CapDAT [YYYY-MM-DD HH:MM:SS]`
`[YYYY-MM-DD HH:MM:SS] [INFO] Opening input file: data/1cwp_full.vdb`
`[YYYY-MM-DD HH:MM:SS] [INFO] Parsing completed successfully`
`[YYYY-MM-DD HH:MM:SS] [INFO] Starting extended structural summary geometry`
`Input file:              data/1cwp_full.vdb`
`Total lines read:        ...`
`Coordinate records seen: ...`
`Accepted atoms:          ...`
`Accepted residues:       ...`
`Internal subunits:       ...`
`Accepted HETATM:         ...`
`Alternate locations:     ...`
`Malformed records:       ...`
`Skipped records:         ...`
`Geometric center (accepted atoms): (...)`
`Coordinate bounds:`
`  x: [...]`
`  y: [...]`
`  z: [...]`
`Axis spans: x=... y=... z=...`
`Runtime (s):             ...`
`[INFO] Run completed successfully [YYYY-MM-DD HH:MM:SS]`
`[YYYY-MM-DD HH:MM:SS]`

## Notes for This Branch

This README intentionally reflects the **current development state of the `geometry` branch**, not a minimal parser-only foundation. That means:
- geometry-analysis options are part of the documented CLI surface,
- CTest is part of the expected development workflow,
- the build instructions explicitly enable tests and Debug mode,
- branch-specific cloning instructions point directly to `geometry`,
- the old dedicated section `Reorientation and final export semantics (v01)` is intentionally removed.

The authoritative references for future README updates on this branch should be:
- `CMakeLists.txt` for build/test targets and options,
- `src/main.cpp` for the real CLI surface,
- the geometry-analysis implementation and tests for stage-level behavior as the pipeline evolves.


For more details look into the `docs` directory:

<https://github.com/tripplab/CapDAT/blob/geometry/docs/technical_spec_and_dev_guide.md>
