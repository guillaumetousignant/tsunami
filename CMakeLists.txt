cmake_minimum_required(VERSION 3.12)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()

project(tsunami VERSION 1.0.0 LANGUAGES C CXX)

add_executable(tsunami 
    src/main.cpp
    src/entities/MeshGeometryUnstructured_t.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(tsunami PUBLIC optimized OpenMP::OpenMP_CXX)
endif()

if(DEFINED ENV{ANOTHER_PATH_TRACER_DIR})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ANOTHER_PATH_TRACER_DIR})
endif()

find_package(APTracer REQUIRED)

target_include_directories(tsunami 
    PRIVATE
        include/
        $ENV{ANOTHER_PATH_TRACER_DIR}/include/ # Should not be needed but what do I know
)

target_link_libraries(tsunami PRIVATE APTracer)
target_compile_features(tsunami PRIVATE cxx_std_11)

option(BUILD_TESTS "Build all tests." OFF)

if (BUILD_TESTS)
    enable_testing()
    add_executable(unit_tests 
        tests/main_test.cpp
        tests/catch.hpp)
    target_link_libraries(unit_tests PRIVATE APTracer)
    target_compile_features(unit_tests PRIVATE cxx_std_11)

    add_test(unit_tests unit_tests)
endif()

install(TARGETS tsunami DESTINATION bin)