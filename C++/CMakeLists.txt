cmake_minimum_required(VERSION 3.16)
project(ASHRAE241_WR LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Warnings and optimization
if(MSVC)
    add_compile_options(/W4 /O2)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -O2)
endif()

# Core library (random_manager + model)
add_library(ashrae241 STATIC
    src/random_manager.cpp
    src/model.cpp
)
target_include_directories(ashrae241 PUBLIC src)

# Simulation executables
add_executable(ecai              src/ecai.cpp)
add_executable(probability_ecai  src/probability_ecai.cpp)
add_executable(probability_scan  src/probability_scan.cpp)
add_executable(single_probability src/single_probability.cpp)

foreach(tgt ecai probability_ecai probability_scan single_probability)
    target_link_libraries(${tgt} PRIVATE ashrae241)
endforeach()

# Analysis executables
add_executable(percentiles  analysis/percentiles.cpp)
add_executable(boxplot      analysis/boxplot.cpp)
target_link_libraries(percentiles PRIVATE ashrae241)
target_link_libraries(boxplot    PRIVATE ashrae241)

# Tests
enable_testing()
add_executable(test_random_manager tests/test_random_manager.cpp)
target_link_libraries(test_random_manager PRIVATE ashrae241)
add_test(NAME random_manager COMMAND test_random_manager)

add_executable(test_model tests/test_model.cpp)
target_link_libraries(test_model PRIVATE ashrae241)
add_test(NAME model COMMAND test_model)

add_executable(test_validation tests/test_validation.cpp)
target_link_libraries(test_validation PRIVATE ashrae241)
add_test(NAME validation COMMAND test_validation)