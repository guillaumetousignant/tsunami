cmake_minimum_required(VERSION 3.12)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()

project(tsunami VERSION 1.0.0 LANGUAGES C CXX)

add_executable(tsunami 
    src/main.cpp
    src/entities/MeshGeometryUnstructured_t.cpp
    src/shapes/MeshUnstructured_t.cpp
    src/shapes/TriangleUnstructured_t.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(tsunami PUBLIC OpenMP::OpenMP_CXX) # Always uses openmp, because library is build with it
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

# Shouldn't be needed START
if(WIN32)
    find_library(PNG REQUIRED
        NAMES png libpng libpng_a libpng_debug libpng_a_debug
        HINTS external_windows/lib) # CHECK make image libraries optional
    find_library(JPEG REQUIRED
        NAMES jpeg libjpeg jpeg_a libjpeg_a
        HINTS external_windows/lib) # CHECK make image libraries optional
    find_library(TIFF REQUIRED
        NAMES tiff libtiff
        HINTS external_windows/lib) # CHECK make image libraries optional
    find_library(GLUT REQUIRED
        NAMES freeglut libfreeglut glut libglut
        HINTS external_windows/lib) # CHECK make image libraries optional

    target_include_directories(tsunami PUBLIC 
        $<INSTALL_INTERFACE:external_windows/include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external_windows/include>)

    SET(PNG_LIBRARY ${PNG})
    SET(JPEG_LIBRARY ${JPEG})
    SET(TIFF_LIBRARY ${TIFF})
    SET(GLUT_LIBRARY ${GLUT})
else()
    find_package(PNG REQUIRED) # CHECK make image libraries optional
    find_package(JPEG REQUIRED)
    find_package(TIFF REQUIRED)
    find_package(GLUT REQUIRED)

    target_include_directories(tsunami PUBLIC 
        ${PNG_INCLUDE_DIR}
        ${JPEG_INCLUDE_DIR}
        ${TIFF_INCLUDE_DIR}
        ${GLUT_INCLUDE_DIRS})
endif()

find_package(OpenGL REQUIRED COMPONENTS OpenGL) # CHECK make optional
add_dependencies(tsunami OpenGL::GL) # Should use OpenGL::OpenGL on newer linux hosts, but doesn't work on windows and android
target_link_libraries(tsunami PUBLIC 
    ${PNG_LIBRARY}
    ${JPEG_LIBRARY}
    ${TIFF_LIBRARY}
    OpenGL::GL
    ${GLUT_LIBRARY}) # Should use OpenGL::OpenGL on newer linux hosts, but doesn't work on windows and android
# Shouldn't be needed END

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