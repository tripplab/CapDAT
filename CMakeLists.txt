cmake_minimum_required(VERSION 3.16)

project(
    CapDAT
    VERSION 0.1.0
    DESCRIPTION "Capsid Data Analysis Toolkit"
    LANGUAGES CXX
)

# ------------------------------------------------------------------------------
# Global build settings
# ------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CAPDAT_BUILD_TESTS "Build CapDAT tests" ON)
option(CAPDAT_ENABLE_WARNINGS "Enable compiler warnings" ON)

# Default build type if not specified by user
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ------------------------------------------------------------------------------
# Source files
# ------------------------------------------------------------------------------

set(CAPDAT_SOURCES
    src/main.cpp
    src/atom.cpp
    src/residue.cpp
    src/chain.cpp
    src/capsid.cpp
    src/pdb_parser.cpp
    src/geometry_symmetry.cpp
    src/export_capsid.cpp
    src/reorientation_workflow.cpp
    src/logger.cpp
    src/structural_summary.cpp
    src/summary_reporter.cpp
    src/timer.cpp
    src/geometry_analysis.cpp
)

# ------------------------------------------------------------------------------
# Executable
# ------------------------------------------------------------------------------

add_executable(capsid_analyzer ${CAPDAT_SOURCES})

target_include_directories(capsid_analyzer
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

# ------------------------------------------------------------------------------
# Compiler warnings
# ------------------------------------------------------------------------------

if(CAPDAT_ENABLE_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(capsid_analyzer PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(capsid_analyzer PRIVATE
            /W4
            /permissive-
        )
    endif()
endif()

# ------------------------------------------------------------------------------
# Future-friendly compile definitions
# ------------------------------------------------------------------------------

target_compile_definitions(capsid_analyzer PRIVATE
    CAPDAT_VERSION="${PROJECT_VERSION}"
)

# ------------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------------

install(TARGETS capsid_analyzer
    RUNTIME DESTINATION bin
)

# ------------------------------------------------------------------------------
# Tests (placeholder for future v01 test wiring)
# ------------------------------------------------------------------------------

if(CAPDAT_BUILD_TESTS)
    enable_testing()

    add_executable(capdat_structural_summary_tests
        tests/structural_summary_tests.cpp
        src/atom.cpp
        src/residue.cpp
        src/chain.cpp
        src/capsid.cpp
        src/pdb_parser.cpp
        src/logger.cpp
        src/structural_summary.cpp
        src/summary_reporter.cpp
    )

    target_include_directories(capdat_structural_summary_tests
        PRIVATE
            ${PROJECT_SOURCE_DIR}/include
    )

    if(CAPDAT_ENABLE_WARNINGS)
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(capdat_structural_summary_tests PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
            )
        elseif(MSVC)
            target_compile_options(capdat_structural_summary_tests PRIVATE
                /W4
                /permissive-
            )
        endif()
    endif()

    add_test(NAME capdat_structural_summary_tests COMMAND capdat_structural_summary_tests)

    add_executable(capdat_geometry_symmetry_tests
        tests/geometry_symmetry_tests.cpp
        src/geometry_symmetry.cpp
    )

    target_include_directories(capdat_geometry_symmetry_tests
        PRIVATE
            ${PROJECT_SOURCE_DIR}/include
    )

    if(CAPDAT_ENABLE_WARNINGS)
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(capdat_geometry_symmetry_tests PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
            )
        elseif(MSVC)
            target_compile_options(capdat_geometry_symmetry_tests PRIVATE
                /W4
                /permissive-
            )
        endif()
    endif()

    add_test(NAME capdat_geometry_symmetry_tests COMMAND capdat_geometry_symmetry_tests)

    add_executable(capdat_reorientation_workflow_tests
        tests/reorientation_workflow_tests.cpp
        src/atom.cpp
        src/residue.cpp
        src/chain.cpp
        src/capsid.cpp
        src/geometry_symmetry.cpp
        src/reorientation_workflow.cpp
        src/logger.cpp
    )
    target_include_directories(capdat_reorientation_workflow_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
    add_test(NAME capdat_reorientation_workflow_tests COMMAND capdat_reorientation_workflow_tests)


    add_executable(capdat_geometry_analysis_tests
        tests/geometry_analysis_tests.cpp
        src/atom.cpp
        src/residue.cpp
        src/chain.cpp
        src/capsid.cpp
        src/geometry_symmetry.cpp
        src/reorientation_workflow.cpp
        src/export_capsid.cpp
        src/logger.cpp
        src/timer.cpp
        src/geometry_analysis.cpp
    )
    target_include_directories(capdat_geometry_analysis_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_compile_definitions(capdat_geometry_analysis_tests PRIVATE CAPDAT_VERSION="${PROJECT_VERSION}")
    add_test(NAME capdat_geometry_analysis_tests COMMAND capdat_geometry_analysis_tests)

    add_executable(capdat_export_capsid_tests
        tests/export_capsid_tests.cpp
        src/atom.cpp
        src/residue.cpp
        src/chain.cpp
        src/capsid.cpp
        src/export_capsid.cpp
        src/logger.cpp
    )
    target_include_directories(capdat_export_capsid_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_compile_definitions(capdat_export_capsid_tests PRIVATE CAPDAT_VERSION="${PROJECT_VERSION}")
    add_test(NAME capdat_export_capsid_tests COMMAND capdat_export_capsid_tests)

    add_test(NAME capdat_cli_export_only
             COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}
                     $<TARGET_FILE:capsid_analyzer> --input 1cwp --export-final cli_export_only.pdb)
    add_test(NAME capdat_cli_fold_to_z_export
             COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}
                     $<TARGET_FILE:capsid_analyzer> --input 1cwp --reorient --align-fold 5_0 --align-axis z --export-final cli_fold_z.pdb)
    add_test(NAME capdat_cli_fold_to_x_export
             COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}
                     $<TARGET_FILE:capsid_analyzer> --input 1cwp --reorient --align-fold 5_0 --align-axis x --export-final cli_fold_x.pdb)
    add_test(NAME capdat_cli_vector_to_y_export
             COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}
                     $<TARGET_FILE:capsid_analyzer> --input 1cwp --reorient --align-vector 1,1,1 --align-axis y --export-final cli_vec_y.pdb)
    add_test(NAME capdat_cli_invalid_combo
             COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}
                     $<TARGET_FILE:capsid_analyzer> --input 1cwp --align-fold 5_0)
    set_tests_properties(capdat_cli_invalid_combo PROPERTIES WILL_FAIL TRUE)
    add_test(NAME capdat_cli_resolve_input_success
             COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}
                     $<TARGET_FILE:capsid_analyzer> --input 1cwp)
    set_tests_properties(capdat_cli_resolve_input_success PROPERTIES
                         PASS_REGULAR_EXPRESSION "data/1cwp_full\\.vdb")

    add_test(NAME capdat_cli_missing_data_dir
             COMMAND bash -lc "cd '${PROJECT_BINARY_DIR}' && output=$('$<TARGET_FILE:capsid_analyzer>' --input 1cwp 2>&1); status=$?; [ $status -ne 0 ] && grep -Fq \"required data directory not found: data/\" <<< \"$output\"")

    add_test(NAME capdat_cli_missing_pdb_file
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && output=$('$<TARGET_FILE:capsid_analyzer>' --input 9zzz 2>&1); status=$?; [ $status -ne 0 ] && grep -Fq \"input file not found for PDBID '9zzz': data/9zzz_full.vdb\" <<< \"$output\"")

    add_test(NAME capdat_cli_reject_path_input
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && output=$('$<TARGET_FILE:capsid_analyzer>' --input data/1cwp_full.vdb 2>&1); status=$?; [ $status -ne 0 ] && grep -Fq -- \"--input expects a PDB identifier only\" <<< \"$output\"")

    add_test(NAME capdat_cli_reject_filename_input
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && output=$('$<TARGET_FILE:capsid_analyzer>' --input 1cwp.vdb 2>&1); status=$?; [ $status -ne 0 ] && grep -Fq -- \"--input expects a PDB identifier only\" <<< \"$output\"")

    add_test(NAME capdat_cli_reject_malformed_input
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && output=$('$<TARGET_FILE:capsid_analyzer>' --input .vdb 2>&1); status=$?; [ $status -ne 0 ] && grep -Fq -- \"--input expects a PDB identifier only\" <<< \"$output\"")

    add_test(NAME capdat_cli_geometry_canonical_root
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && rm -rf results/1cwp/2_0 && '$<TARGET_FILE:capsid_analyzer>' --input 1cwp --geometry-analysis >/tmp/capdat_geom_root.log 2>&1 && test -d results && test -d results/1cwp && test -d results/1cwp/2_0 && test -f results/1cwp/2_0/stage_run_summary.json")
    add_test(NAME capdat_cli_geometry_fold_name_success
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && rm -rf results/1cwp/5_0 && '$<TARGET_FILE:capsid_analyzer>' --input 1cwp --geometry-analysis --geometry_fold_name 5_0 >/tmp/capdat_geom_foldname.log 2>&1 && test -d results/1cwp/5_0 && test -f results/1cwp/5_0/stage_run_summary.json")
    add_test(NAME capdat_cli_geometry_fold_name_invalid
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && output=$('$<TARGET_FILE:capsid_analyzer>' --input 1cwp --geometry-analysis --geometry_fold_name bad_fold 2>&1); status=$?; [ $status -ne 0 ] && grep -Fq \"Invalid canonical fold name\" <<< \"$output\"")
    add_test(NAME capdat_cli_geometry_fold_name_mutual_exclusion
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && output=$('$<TARGET_FILE:capsid_analyzer>' --input 1cwp --geometry-analysis --geometry_fold_name 5_0 --geometry_fold_type 5 --geometry_fold_index 0 2>&1); status=$?; [ $status -ne 0 ] && grep -Fq \"mutually exclusive\" <<< \"$output\"")
    add_test(NAME capdat_cli_geometry_clean_results_false
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && mkdir -p results/1cwp/5_0 && echo keep > results/1cwp/5_0/keep.txt && '$<TARGET_FILE:capsid_analyzer>' --input 1cwp --geometry-analysis --geometry_fold_name 5_0 --clean-results-dir false >/tmp/capdat_geom_clean_false.log 2>&1 && test -f results/1cwp/5_0/keep.txt")
    add_test(NAME capdat_cli_geometry_clean_results_true
             COMMAND bash -lc "cd '${PROJECT_SOURCE_DIR}' && mkdir -p results/1cwp/5_0/old/nested && echo stale > results/1cwp/5_0/old/nested/stale.txt && '$<TARGET_FILE:capsid_analyzer>' --input 1cwp --geometry-analysis --geometry_fold_name 5_0 --clean-results-dir true >/tmp/capdat_geom_clean_true.log 2>&1 && test -d results && test -d results/1cwp && test ! -e results/1cwp/5_0/old && test -f results/1cwp/5_0/stage_run_summary.json")

endif()

# ------------------------------------------------------------------------------
# Build summary
# ------------------------------------------------------------------------------

message(STATUS "Configuring ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build tests: ${CAPDAT_BUILD_TESTS}")
message(STATUS "Warnings enabled: ${CAPDAT_ENABLE_WARNINGS}")
