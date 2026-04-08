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
    src/accepted_structure_writer.cpp
    src/logger.cpp
    src/structural_summary.cpp
    src/summary_reporter.cpp
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
endif()

# ------------------------------------------------------------------------------
# Build summary
# ------------------------------------------------------------------------------

message(STATUS "Configuring ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build tests: ${CAPDAT_BUILD_TESTS}")
message(STATUS "Warnings enabled: ${CAPDAT_ENABLE_WARNINGS}")
