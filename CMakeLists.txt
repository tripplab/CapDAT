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

option(CAPDAT_BUILD_TESTS "Build CapDAT tests" OFF)
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
    src/fold_patch_analysis.cpp
    src/fold_patch_frame.cpp
    src/fold_patch_extraction.cpp
    src/fold_patch_contact.cpp
    src/fold_patch_surface_fit.cpp
    src/fold_patch_metrics.cpp
    src/fold_patch_reporter.cpp
    src/fold_patch_export.cpp
    src/timer.cpp
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

    add_executable(capdat_fold_patch_stage2_tests
        tests/fold_patch_stage2_tests.cpp
        src/atom.cpp
        src/residue.cpp
        src/chain.cpp
        src/capsid.cpp
        src/pdb_parser.cpp
        src/logger.cpp
        src/geometry_symmetry.cpp
        src/fold_patch_analysis.cpp
        src/fold_patch_frame.cpp
        src/fold_patch_extraction.cpp
        src/fold_patch_contact.cpp
        src/fold_patch_surface_fit.cpp
        src/fold_patch_metrics.cpp
        src/fold_patch_reporter.cpp
        src/fold_patch_export.cpp
    )
    target_include_directories(capdat_fold_patch_stage2_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_compile_definitions(capdat_fold_patch_stage2_tests PRIVATE CAPDAT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
    add_test(NAME capdat_fold_patch_stage2_tests COMMAND capdat_fold_patch_stage2_tests)

    add_test(NAME capdat_cli_export_only
             COMMAND capsid_analyzer --input ${PROJECT_SOURCE_DIR}/data/1cwp_full.vdb --export-final cli_export_only.pdb)
    add_test(NAME capdat_cli_fold_to_z_export
             COMMAND capsid_analyzer --input ${PROJECT_SOURCE_DIR}/data/1cwp_full.vdb --reorient --align-fold 5_0 --align-axis z --export-final cli_fold_z.pdb)
    add_test(NAME capdat_cli_fold_to_x_export
             COMMAND capsid_analyzer --input ${PROJECT_SOURCE_DIR}/data/1cwp_full.vdb --reorient --align-fold 5_0 --align-axis x --export-final cli_fold_x.pdb)
    add_test(NAME capdat_cli_vector_to_y_export
             COMMAND capsid_analyzer --input ${PROJECT_SOURCE_DIR}/data/1cwp_full.vdb --reorient --align-vector 1,1,1 --align-axis y --export-final cli_vec_y.pdb)
    add_test(NAME capdat_cli_invalid_combo
             COMMAND capsid_analyzer --input ${PROJECT_SOURCE_DIR}/data/1cwp_full.vdb --align-fold 5_0)
    set_tests_properties(capdat_cli_invalid_combo PROPERTIES WILL_FAIL TRUE)
endif()

# ------------------------------------------------------------------------------
# Build summary
# ------------------------------------------------------------------------------

message(STATUS "Configuring ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build tests: ${CAPDAT_BUILD_TESTS}")
message(STATUS "Warnings enabled: ${CAPDAT_ENABLE_WARNINGS}")
