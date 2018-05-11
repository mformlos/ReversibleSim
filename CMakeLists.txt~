cmake_minimum_required(VERSION 2.8.11)
project(HYBRIDSIM)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

list( APPEND CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -g -march=native  -ffast-math  -Wall -O3 -Wfatal-errors ")

include_directories ("${PROJECT_BINARY_DIR}")

include_directories ("${PROJECT_SOURCE_DIR}/lib")
add_subdirectory (lib)

include_directories ("${PROJECT_SOURCE_DIR}/src")
add_subdirectory (src)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Particles")
add_subdirectory(lib/Particles)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")
add_subdirectory(lib/Hydrodynamics)

include_directories ("${PROJECT_SOURCE_DIR}/lib/MolecularDynamics")
add_subdirectory(lib/MolecularDynamics)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")
add_subdirectory(lib/Analysis)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Exceptions")
add_subdirectory(lib/Exceptions)

#find_package(Eigen3 REQUIRED)


