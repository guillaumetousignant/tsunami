cmake_minimum_required(VERSION 3.12)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()

project(tsunami VERSION 1.0.0 LANGUAGES C CXX)

add_executable(tsunami 
    src/main.cpp
    src/entities/MeshGeometryUnstructured_t.cpp
    src/shapes/MeshUnstructured_t.cpp
    src/shapes/TriangleUnstructured_t.cpp
    src/cameras/CamIndex_t.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(tsunami PUBLIC OpenMP::OpenMP_CXX) # Always uses openmp, because library is build with it
endif()

find_package(OpenGL REQUIRED COMPONENTS OpenGL) # Shouldn't be needed, should be inherited from APTracer...
find_package(APTracer REQUIRED)

target_include_directories(tsunami 
    PRIVATE
        include
)

target_link_libraries(tsunami PRIVATE APTracer::APTracer)

target_compile_features(tsunami PRIVATE cxx_std_11)
set_property(TARGET tsunami PROPERTY INTERPROCEDURAL_OPTIMIZATION True) # Always on, because library is built with it

option(BUILD_TESTS "Build all tests." OFF)

if (BUILD_TESTS)
    enable_testing()
    add_executable(unit_tests 
        tests/main_test.cpp
        tests/catch.hpp)
    target_link_libraries(unit_tests PRIVATE APTracer::APTracer)
    target_compile_features(unit_tests PRIVATE cxx_std_11)

    add_test(unit_tests unit_tests)
endif()

install(TARGETS tsunami DESTINATION bin)