cmake_minimum_required(VERSION 3.16)
project(ecai_cpp2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
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

# CLI-only C++2 implementation. Must not reference or link against the
# existing C++ tree (see IMPLEMENTATION_PLAN.md section 1). From Stage 12 the
# entry point resolves a run and drives the Stage 3-11 pipeline, so it links
# every implementation source.
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
    src/output.cpp
    src/main.cpp
)
target_include_directories(ecai_cpp2 PRIVATE include)
target_compile_definitions(ecai_cpp2 PRIVATE
    ECAI_CPP2_PROTECTED_LEGACY_TREE="${CMAKE_CURRENT_SOURCE_DIR}/../C++"
)

enable_testing()

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
# target; depends on lhs, seed_manager, and distributions).
add_executable(ecai_cpp2_sampling_plan_tests
    src/lhs.cpp
    src/seed_manager.cpp
    src/distributions.cpp
    src/sampling_plan.cpp
    tests/test_sampling_plan.cpp
)
target_include_directories(ecai_cpp2_sampling_plan_tests PRIVATE include)

add_test(NAME ecai_cpp2_sampling_plan_tests COMMAND ecai_cpp2_sampling_plan_tests)

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
    src/model.cpp
    src/simulation.cpp
    tests/test_simulation.cpp
)
target_include_directories(ecai_cpp2_simulation_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_simulation_tests PRIVATE
    ECAI_CPP2_SIMULATION_HEADER_PATH="${CMAKE_CURRENT_SOURCE_DIR}/include/simulation.h"
    ECAI_CPP2_SIMULATION_SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/src/simulation.cpp"
)

add_test(NAME ecai_cpp2_simulation_tests COMMAND ecai_cpp2_simulation_tests)

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
    tests/test_cli.cpp
)
target_include_directories(ecai_cpp2_cli_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_cli_tests PRIVATE
    ECAI_CPP2_EXECUTABLE="$<TARGET_FILE:ecai_cpp2>"
)
add_dependencies(ecai_cpp2_cli_tests ecai_cpp2)

add_test(NAME ecai_cpp2_cli_tests COMMAND ecai_cpp2_cli_tests)

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
    src/output.cpp
    tests/test_output.cpp
)
target_include_directories(ecai_cpp2_output_tests PRIVATE include)
target_compile_definitions(ecai_cpp2_output_tests PRIVATE
    ECAI_CPP2_EXECUTABLE="$<TARGET_FILE:ecai_cpp2>"
    ECAI_CPP2_PROTECTED_LEGACY_TREE="${CMAKE_CURRENT_SOURCE_DIR}/../C++"
)
add_dependencies(ecai_cpp2_output_tests ecai_cpp2)
add_test(NAME ecai_cpp2_output_tests COMMAND ecai_cpp2_output_tests)

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
