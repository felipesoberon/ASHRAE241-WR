cmake_minimum_required(VERSION 3.16)
project(ecai_cpp2 LANGUAGES CXX)

file(GLOB ECAI_SOURCE_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${ECAI_SOURCE_FILES})
set(ECAI_SOURCE_HASH_INPUT "")
foreach(source IN LISTS ECAI_SOURCE_FILES)
    file(SHA256 "${source}" source_hash)
    get_filename_component(source_name "${source}" NAME)
    string(APPEND ECAI_SOURCE_HASH_INPUT "${source_name}:${source_hash}\n")
endforeach()
string(SHA256 ECAI_SOURCE_CONTENT_HASH "${ECAI_SOURCE_HASH_INPUT}")

execute_process(COMMAND git -C ${CMAKE_CURRENT_SOURCE_DIR} rev-parse HEAD
    OUTPUT_VARIABLE ECAI_SOURCE_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT ECAI_SOURCE_REVISION)
    set(ECAI_SOURCE_REVISION "unknown")
endif()
execute_process(COMMAND git -C ${CMAKE_CURRENT_SOURCE_DIR} diff --quiet HEAD -- .
    RESULT_VARIABLE ECAI_TRACKED_DIFF_RESULT ERROR_QUIET)
if(ECAI_TRACKED_DIFF_RESULT EQUAL 0)
    set(ECAI_SOURCE_TRACKED_DIRTY "false")
else()
    set(ECAI_SOURCE_TRACKED_DIRTY "true_or_unknown")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to an optimized build with debug info when the caller does not
# choose one. The reproducible floating-point flags below apply to every
# configuration; no -Ofast / -ffast-math / -march=native is ever added.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# Warnings and reproducible floating-point flags, per docs/MODEL_SPECIFICATION.md
# section 11.2 (no fast-math, no FP contraction / FMA fusing).
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -fno-fast-math
    -ffp-contract=off
)

# CLI-only current C++ implementation. Must not reference or link against the
# legacy C++_legacy tree (see IMPLEMENTATION_PLAN.md section 1). From Stage 12
# the entry point resolves a run and drives the Stage 3-11 pipeline, so it
# links every implementation source.
add_executable(ecai_cpp2
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/model.cpp
    src/simulation.cpp
    src/statistics.cpp
    src/cli.cpp
    src/qer_direct.cpp
    src/output.cpp
    src/paired_output.cpp
    src/main.cpp
)
target_include_directories(ecai_cpp2 PRIVATE include)
target_compile_definitions(ecai_cpp2 PRIVATE
    ECAI_CPP2_PROTECTED_LEGACY_TREE="${CMAKE_CURRENT_SOURCE_DIR}/../C++_legacy"
    ECAI_SOURCE_REVISION="${ECAI_SOURCE_REVISION}"
    ECAI_SOURCE_TRACKED_DIRTY="${ECAI_SOURCE_TRACKED_DIRTY}"
    ECAI_SOURCE_CONTENT_HASH="${ECAI_SOURCE_CONTENT_HASH}"
    ECAI_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}"
)

enable_testing()

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME ecai_cpp2_paired_reader_tests
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_paired_reader.py
                ${CMAKE_BINARY_DIR})
endif()

add_executable(ecai_cpp2_tests tests/test_main.cpp)
target_compile_definitions(ecai_cpp2_tests PRIVATE
    ECAI_CPP2_EXECUTABLE="$<TARGET_FILE:ecai_cpp2>"
)

add_test(NAME ecai_cpp2_main_tests COMMAND ecai_cpp2_tests)

# Stage 2: category data (separate CTest target, independent of the CLI
# executable test above).
add_executable(ecai_cpp2_category_data_tests
    src/category_data.cpp
    tests/test_category_data.cpp
)
target_include_directories(ecai_cpp2_category_data_tests PRIVATE include)

add_test(NAME ecai_cpp2_category_data_tests COMMAND ecai_cpp2_category_data_tests)

# Stage 3: exact-size one-dimensional LHS engine (separate CTest target,
# independent of distribution transforms).
add_executable(ecai_cpp2_lhs_tests
    src/lhs.cpp
    tests/test_lhs.cpp
)
target_include_directories(ecai_cpp2_lhs_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_lhs_tests PRIVATE
    ECAI_CPP2_LHS_SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/src/lhs.cpp"
)

add_test(NAME ecai_cpp2_lhs_tests COMMAND ecai_cpp2_lhs_tests)

# Stage 4: deterministic per-replicate/per-stream seed derivation (separate
# CTest target, independent of LHS/category data).
add_executable(ecai_cpp2_seed_manager_tests
    src/seed_manager.cpp
    tests/test_seed_manager.cpp
)
target_include_directories(ecai_cpp2_seed_manager_tests PRIVATE include)

add_test(NAME ecai_cpp2_seed_manager_tests COMMAND ecai_cpp2_seed_manager_tests)

# Stage 5: inverse CDFs / Jones parameter distribution transforms (separate
# CTest target, independent of LHS/seed/category data).
add_executable(ecai_cpp2_distributions_tests
    src/distributions.cpp
    tests/test_distributions.cpp
)
target_include_directories(ecai_cpp2_distributions_tests PRIVATE include)

add_test(NAME ecai_cpp2_distributions_tests COMMAND ecai_cpp2_distributions_tests)

# Stage 6: binomial pre-sampling and exact draw-count plan (separate CTest
# target; depends on lhs, seed_manager, distributions, and qer_direct --
# sampling_plan.cpp's package 3a family-aware path calls
# qer_direct::sample_direct_qer_log10normal).
add_executable(ecai_cpp2_sampling_plan_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/sampling_plan.cpp
    src/qer_direct.cpp
    tests/test_sampling_plan.cpp
)
target_include_directories(ecai_cpp2_sampling_plan_tests PRIVATE include)

add_test(NAME ecai_cpp2_sampling_plan_tests COMMAND ecai_cpp2_sampling_plan_tests)

# Package 3a: family-aware Stage 7 bank allocation and accounting (separate
# CTest target; depends on lhs, seed_manager, distributions, category_data,
# sampling_plan, and qer_direct for the direct-QER transform/tables).
add_executable(ecai_cpp2_family_bank_allocation_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/qer_direct.cpp
    tests/test_family_bank_allocation.cpp
)
target_include_directories(ecai_cpp2_family_bank_allocation_tests PRIVATE include)

add_test(NAME ecai_cpp2_family_bank_allocation_tests COMMAND ecai_cpp2_family_bank_allocation_tests)

# Stage 8: Jones QER calculation, a pure function of supplied sampled
# infector values (separate CTest target; depends only on sampling_plan.h
# for the InfectorParameters type, not on any sampling_plan.cpp/lhs.cpp/
# seed_manager.cpp/distributions.cpp/category_data.cpp linkage).
add_executable(ecai_cpp2_model_tests
    src/model.cpp
    tests/test_model.cpp
)
target_include_directories(ecai_cpp2_model_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_model_tests PRIVATE
    ECAI_CPP2_MODEL_SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/src/model.cpp"
)

add_test(NAME ecai_cpp2_model_tests COMMAND ecai_cpp2_model_tests)

# Stage 10: one-replicate simulation (plan -> banks -> stable assignments ->
# QER -> ECAi -> N index-ordered outputs -> frozen percentile). Separate
# CTest target; links the Stage 3-9 sources it orchestrates.
add_executable(ecai_cpp2_simulation_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/qer_direct.cpp
    src/model.cpp
    src/simulation.cpp
    src/statistics.cpp
    tests/test_simulation.cpp
)
target_include_directories(ecai_cpp2_simulation_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_simulation_tests PRIVATE
    ECAI_CPP2_SIMULATION_HEADER_PATH="${CMAKE_CURRENT_SOURCE_DIR}/include/simulation.h"
    ECAI_CPP2_SIMULATION_SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/src/simulation.cpp"
)

add_test(NAME ecai_cpp2_simulation_tests COMMAND ecai_cpp2_simulation_tests)

# Package 3b: direct-QER family integration into the Stage 10 ECAi/Probability
# simulation loops (separate CTest target; links the same Stage 3-9 sources
# as ecai_cpp2_simulation_tests, since it drives simulation.cpp's family-aware
# entry point alongside the unmodified full-Jones one).
add_executable(ecai_cpp2_qer_family_simulation_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/qer_direct.cpp
    src/model.cpp
    src/simulation.cpp
    src/statistics.cpp
    tests/test_qer_family_simulation.cpp
)
target_include_directories(ecai_cpp2_qer_family_simulation_tests PRIVATE include)

add_test(NAME ecai_cpp2_qer_family_simulation_tests COMMAND ecai_cpp2_qer_family_simulation_tests)

# Stage 11: across-replicate statistics (mean / sample SD / SE / Student-t CI
# / raw min-max / ceiling-to-5 comparator) and the multi-replicate
# orchestration entry point. Separate CTest target; links the Stage 3-10
# sources the orchestration drives.
add_executable(ecai_cpp2_statistics_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/qer_direct.cpp
    src/model.cpp
    src/simulation.cpp
    src/statistics.cpp
    tests/test_statistics.cpp
)
target_include_directories(ecai_cpp2_statistics_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_statistics_tests PRIVATE
    ECAI_CPP2_STATISTICS_HEADER_PATH="${CMAKE_CURRENT_SOURCE_DIR}/include/statistics.h"
    ECAI_CPP2_STATISTICS_SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/src/statistics.cpp"
    ECAI_CPP2_SIMULATION_HEADER_PATH="${CMAKE_CURRENT_SOURCE_DIR}/include/simulation.h"
)

add_test(NAME ecai_cpp2_statistics_tests COMMAND ecai_cpp2_statistics_tests)

# Stage 12: strict CLI parsing, category/override resolution, and the
# executable command matrix. Separate CTest target; links the Stage 3-11
# sources the resolved run drives, plus cli.
add_executable(ecai_cpp2_cli_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/model.cpp
    src/simulation.cpp
    src/statistics.cpp
    src/cli.cpp
    src/qer_direct.cpp
    tests/test_cli.cpp
)
target_include_directories(ecai_cpp2_cli_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_cli_tests PRIVATE
    ECAI_CPP2_EXECUTABLE="$<TARGET_FILE:ecai_cpp2>"
)
add_dependencies(ecai_cpp2_cli_tests ecai_cpp2)

add_test(NAME ecai_cpp2_cli_tests COMMAND ecai_cpp2_cli_tests)

# Work package 2: --qer-distribution CLI selector parsing and the
# ResolvedQerSelection descriptor resolve_run builds from it (separate CTest
# target, independent of tests/test_cli.cpp; links the same Stage 3-11 plus
# cli/qer_direct sources ecai_cpp2_cli_tests links, since resolve_run drives
# the same category/sampling-plan resolution path).
add_executable(ecai_cpp2_qer_cli_scenario_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/model.cpp
    src/simulation.cpp
    src/statistics.cpp
    src/cli.cpp
    src/qer_direct.cpp
    tests/test_qer_cli_scenario.cpp
)
target_include_directories(ecai_cpp2_qer_cli_scenario_tests PRIVATE include)

add_test(NAME ecai_cpp2_qer_cli_scenario_tests COMMAND ecai_cpp2_qer_cli_scenario_tests)

# Package 3b: executable-level regression test proving --qer-distribution
# actually propagates into the executed run (not merely into the CLI-layer
# ResolvedQerSelection descriptor tests/test_qer_cli_scenario.cpp already
# covers). Invokes the built ecai_cpp2 executable directly, like
# ecai_cpp2_cli_tests's executable command matrix.
add_executable(ecai_cpp2_qer_cli_executable_tests
    tests/test_qer_cli_executable.cpp
)
target_include_directories(ecai_cpp2_qer_cli_executable_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_qer_cli_executable_tests PRIVATE
    ECAI_CPP2_EXECUTABLE="$<TARGET_FILE:ecai_cpp2>"
)
add_dependencies(ecai_cpp2_qer_cli_executable_tests ecai_cpp2)

add_test(NAME ecai_cpp2_qer_cli_executable_tests COMMAND ecai_cpp2_qer_cli_executable_tests)

# Stage 13: metadata-grounded CSV composition and race-safe output writer.
add_executable(ecai_cpp2_output_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/category_data.cpp
    src/sampling_plan.cpp
    src/model.cpp
    src/simulation.cpp
    src/statistics.cpp
    src/cli.cpp
    src/qer_direct.cpp
    src/output.cpp
    tests/test_output.cpp
)
target_include_directories(ecai_cpp2_output_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_output_tests PRIVATE
    ECAI_CPP2_EXECUTABLE="$<TARGET_FILE:ecai_cpp2>"
    ECAI_CPP2_PROTECTED_LEGACY_TREE="${CMAKE_CURRENT_SOURCE_DIR}/../C++_legacy"
)
add_dependencies(ecai_cpp2_output_tests ecai_cpp2)
add_test(NAME ecai_cpp2_output_tests COMMAND ecai_cpp2_output_tests)

# Work package 1: typed direct-QER family/profile selection and the
# immutable Jones-fitted/Mikszewski parameter tables (separate CTest target;
# depends only on distributions.h/.cpp for the shared log10-normal
# transform, not on category_data/lhs/seed_manager/sampling_plan/model).
add_executable(ecai_cpp2_qer_direct_tests
    src/distributions.cpp
    src/qer_direct.cpp
    tests/test_qer_direct.cpp
)
target_include_directories(ecai_cpp2_qer_direct_tests PRIVATE include)

add_test(NAME ecai_cpp2_qer_direct_tests COMMAND ecai_cpp2_qer_direct_tests)

# Stage 14: bounded end-to-end scientific validation against reviewed fixtures.
add_executable(ecai_cpp2_end_to_end_tests
    src/model.cpp
    tests/test_end_to_end.cpp
)
target_include_directories(ecai_cpp2_end_to_end_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_end_to_end_tests PRIVATE
    ECAI_CPP2_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    ECAI_CPP2_EXECUTABLE="$<TARGET_FILE:ecai_cpp2>"
)
add_dependencies(ecai_cpp2_end_to_end_tests ecai_cpp2)
add_test(NAME ecai_cpp2_end_to_end_tests COMMAND ecai_cpp2_end_to_end_tests)
