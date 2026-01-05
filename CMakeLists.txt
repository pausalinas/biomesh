cmake_minimum_required(VERSION 3.14)
project(BioMesh VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Find OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - parallel mesh generation enabled")
else()
    message(STATUS "OpenMP not found - mesh generation will run serially")
endif()

# Include directories
include_directories(include)

# Create library for the core functionality
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_library(biomesh_lib ${SOURCES})
target_include_directories(biomesh_lib PUBLIC include)

# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(biomesh_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Create executable for example usage
add_executable(biomesh_example examples/main.cpp)
target_link_libraries(biomesh_example biomesh_lib)

# Create executable for voxel demonstration
add_executable(voxel_demo examples/voxel_demo.cpp)
target_link_libraries(voxel_demo biomesh_lib)

# Create executable for empty voxel to GiD export
add_executable(empty_voxel_to_gid examples/empty_voxel_to_gid.cpp)
target_link_libraries(empty_voxel_to_gid biomesh_lib)

# Create executable for occupied voxel to GiD export
add_executable(occupied_voxel_to_gid examples/occupied_voxel_to_gid.cpp)
target_link_libraries(occupied_voxel_to_gid biomesh_lib)

# Create executable for filter demonstration
add_executable(filter_demo examples/filter_demo.cpp)
target_link_libraries(filter_demo biomesh_lib)

# Create executable for filter workflow example
add_executable(filter_workflow_example examples/filter_workflow_example.cpp)
target_link_libraries(filter_workflow_example biomesh_lib)

# Enable testing
enable_testing()

# Find GoogleTest
find_package(GTest QUIET)
if(GTest_FOUND)
    # Create test executable
    file(GLOB_RECURSE TEST_SOURCES "tests/*.cpp")
    add_executable(biomesh_tests ${TEST_SOURCES})
    target_link_libraries(biomesh_tests biomesh_lib GTest::gtest GTest::gtest_main)
    
    # Add test
    add_test(NAME AllTests COMMAND biomesh_tests)
else()
    message(WARNING "GoogleTest not found. Tests will not be built.")
endif()