# Capsid Data Analysis Toolkit (CapDAT)
High-performance software tool for the structural analysis of viral capsids from atomic coordinate data in Protein Data Bank (PDB) format

**CapDAT (Capsid Data Analysis Toolkit)** is a C++ command-line software project for the structural analysis of viral capsids from atomic coordinate files in PDB-like format. The current **v01 foundation release** focuses on building a reliable and extensible structural core rather than advanced scientific analysis. Its main purpose is to read an input structure file, parse atomic coordinate records, filter the content so that only capsid protein atoms are retained internally, reconstruct the molecular hierarchy, and report summary statistics about the parsed assembly.

At this stage, CapDAT implements a lightweight object-oriented architecture centered on the classes **Atom**, **Residue**, **Chain**, **Capsid**, **PdbParser**, **Logger**, and **Timer**. The internal hierarchy is reconstructed as **atoms grouped into residues, residues grouped into internally reconstructed subunits, and subunits grouped into the full capsid object**. A key design decision is that the original one-letter PDB chain label is preserved only as metadata and is **not** treated as a globally unique identifier, since large capsid structures may reuse chain labels across many independent proteins. This makes the software more suitable as a long-term base for capsid-specific structural analysis.

The program currently provides a command-line interface with support for input-file selection, optional log-file output, verbosity control, optional inclusion of `HETATM` records, help output, and version reporting. During execution, it produces informative runtime messages and a final terminal summary that includes the input file name, total lines read, coordinate records detected, accepted atoms, accepted residues, internally reconstructed subunits, alternate-location counts, skipped or malformed records, and total runtime. Logging and timing are already integrated from the first version so runs remain traceable and easy to debug.

The main engineering goal of CapDAT v01 is to establish a trustworthy, maintainable, and performance-aware parsing foundation that can support future analytical modules. Planned later extensions may include geometric descriptors, inter-subunit distance analysis, symmetry-related grouping, shell thickness and radius distributions, local curvature and anisotropy metrics, batch processing, and OpenMP-enabled acceleration. In that sense, the present release is intentionally conservative: it prioritizes **correct parsing, clear hierarchy reconstruction, modular software organization, and future extensibility** over premature scientific complexity.

# Requirements Description

CapDAT v01 is intended as the foundation release of a C++ command-line software project for the structural analysis of viral capsids from atomic coordinate files in standard fixed-column PDB-like format. The software must be able to read a single input structure file from the command line, identify supported coordinate records, and parse the atomic information needed to reconstruct the molecular hierarchy of the assembly. At minimum, the program must support `ATOM` records and may optionally support `HETATM` records under explicit configuration. It must extract core structural fields such as atom identifiers, residue identifiers, original chain labels, Cartesian coordinates, occupancy, temperature factor, element, and charge when available.

A central functional requirement is that the software must internally reconstruct the hierarchy as **Atom -> Residue -> internal subunit (`Chain`) -> Capsid**. This hierarchy is essential because the original one-letter PDB chain identifier cannot be assumed to be globally unique in large capsid assemblies. The program must therefore preserve the raw chain label as metadata while assigning its own internal subunit identifiers for grouping and reporting. Residues must be grouped consistently within each reconstructed subunit, and the complete assembly must be represented through a top-level capsid object that owns the full parsed structure.

The parser must be robust enough to process large structures while handling malformed or unsupported content safely. It must distinguish fatal errors from recoverable warnings, fail clearly when the input file is missing, unreadable, empty, or contains no valid coordinate records, and avoid silent failure. The software must also support the internal filtering of non-capsid components, so that atoms belonging to molecules such as DNA, RNA, water, ions, ligands, or other non-protein entities can be ignored during v01 parsing under parser or CLI policy. Basic support for alternate locations and blank chain identifiers must also be present without causing crashes.

CapDAT must provide a usable and predictable command-line interface from the first version. At minimum, the executable must accept an input file path, provide help and version messages, and optionally support a log-file path, increased terminal verbosity, quiet mode, and inclusion of `HETATM` records. Runtime messaging is a required feature: the program must report key execution stages such as startup, input-file opening, parsing progress, warning conditions, completion, and total runtime. Logging support must be available so that runs can be traced later, with message levels such as INFO, WARNING, ERROR, and DEBUG.

The output requirements for v01 are intentionally modest but precise. The program must print a final terminal summary that includes, at minimum, the input file name, total lines read, total coordinate records detected, accepted atoms, accepted residues, internally reconstructed subunits, accepted hetero atoms when enabled, alternate-location counts when relevant, skipped or malformed records, and total runtime. These outputs are required both for user transparency and for basic validation of parser correctness. Internally, the software must remain modular, lightweight, and extensible so that later analytical features such as geometric measurements, symmetry analysis, anisotropy metrics, and high-throughput processing can be added without redesigning the structural core.

# How to Clone CapDAT from GitHub

To obtain the CapDAT source code from GitHub, open a terminal and move to the directory where you want the project folder to be created. Then run `git clone git@github.com:tripplab/CapDAT.git` if you use SSH authentication with GitHub, or `git clone https://github.com/tripplab/CapDAT.git` if you prefer HTTPS. This command will create a local directory named `CapDAT` containing the full repository history and all project files. After the clone finishes, enter the project directory with `cd CapDAT`. You can verify that the repository was cloned correctly by running `git status`, which should report that you are on branch `main` with a clean working tree. If you want to confirm the configured remote, run `git remote -v`, which should show the GitHub repository URL for both fetch and push operations.

# How to Define an Environment for CapDAT

To define a reproducible environment for CapDAT, create a dedicated **micromamba** environment containing the minimal tools needed to configure, compile, and test the software. A practical choice is to create an environment named `capdat` with CMake, a C++ compiler toolchain, `make`, and `ninja`. This can be done with the command `micromamba create -n capdat -c conda-forge python=3.11 cmake cxx-compiler make ninja -y`. After creation, activate the environment with `micromamba activate capdat`. If shell activation is not yet configured, first initialize the shell hook with `eval "$(micromamba shell hook --shell bash)"` and then activate the environment. Once activated, verify that the required tools are available by checking `cmake --version`, `c++ --version`, `make --version`, and `ninja --version`. Using a dedicated environment is recommended because it isolates the CapDAT toolchain from system-level software, improves reproducibility across nodes and machines, and ensures that the same build tools are used consistently during development and testing.

# How to Build CapDAT

To build CapDAT, first make sure your environment is active or that the required build tools are available, especially `cmake` and a C++ compiler. From the root directory of the cloned repository, configure the project with CMake by running `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug` if you want to use the Ninja backend, or `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` if you prefer the default Makefile-based generator. This step creates the `build/` directory and generates all necessary build files. Once configuration completes successfully, compile the program with `cmake --build build`. If the build succeeds, the executable `capsid_analyzer` will be created inside the `build/` directory. You can confirm that the build worked by running `./build/capsid_analyzer --help`, which should display the command-line usage message for the program.

# How to Test the CapDAT Build

To test that the CapDAT build completed successfully, first confirm that the executable exists in the `build/` directory after compilation. A simple functional check is to run `./build/capsid_analyzer --help`, which should print the program help message, and `./build/capsid_analyzer --version`, which should report the current software version. After these interface checks, run the program on a sample structure file included in the repository with `./build/capsid_analyzer -i data/1cwp_full.vdb --verbose --log run.log`. A successful test should produce runtime log messages, generate a final summary in the terminal, and optionally write a log file if requested. The summary should include parsed atom, residue, and internal subunit counts, together with skipped-record statistics and total runtime. If these steps complete without errors, the build can be considered functionally validated at the v01 foundation level.

# Example Output from the CapDAT Build Test

When the build was tested with the command `./build/capsid_analyzer -i data/1cwp_full.vdb --verbose --log run.log`, the program produced the following output:

`[2026-04-07 22:04:20] [INFO] Starting CapDAT`  
`[2026-04-07 22:04:20] [INFO] Input file: data/1cwp_full.vdb`  
`[2026-04-07 22:04:20] [INFO] Log file: run.log`  
`[2026-04-07 22:04:20] [INFO] Include HETATM: false`  
`[2026-04-07 22:04:20] [INFO] Protein only: true`  
`[2026-04-07 22:04:20] [INFO] Opening input file: data/1cwp_full.vdb`  
`[2026-04-07 22:04:21] [INFO] Parsing completed successfully`  
`[2026-04-07 22:04:21] [INFO] Accepted atoms: 214440`  
`[2026-04-07 22:04:21] [INFO] Accepted residues: 28620`  
`[2026-04-07 22:04:21] [INFO] Internal subunits: 180`  

`=== CapDAT Summary ===`  
`Input file:              data/1cwp_full.vdb`  
`Total lines read:        227592`  
`Coordinate records seen: 227040`  
`Accepted atoms:          214440`  
`Accepted residues:       28620`  
`Internal subunits:       180`  
`Accepted HETATM:         0`  
`Alternate locations:     0`  
`Malformed records:       0`  
`Skipped records:         12600`  
`Runtime (s):             0.384428`  
`[2026-04-07 22:04:21] [INFO] Run completed successfully`

This output indicates that the CapDAT v01 build completed and executed correctly on the sample input file, successfully parsing the structure, reconstructing the internal hierarchy, producing summary statistics, and completing the run without fatal errors.





