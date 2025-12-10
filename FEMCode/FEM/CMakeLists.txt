cmake_minimum_required(VERSION 3.16)
project(FEM LANGUAGES CXX)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Include paths for your own headers ---
# Everything in FEM/include will be available as #include "FEMobject.h"
include_directories(${CMAKE_SOURCE_DIR}/include)

# --- EIGEN (Homebrew detection for Intel + M1/M2) ---
if(EXISTS /usr/local/include/eigen3)
    message(STATUS "Using Eigen in /usr/local/include/eigen3")
    include_directories(/usr/local/include/eigen3)
elseif(EXISTS /opt/homebrew/include/eigen3)
    message(STATUS "Using Eigen in /opt/homebrew/include/eigen3")
    include_directories(/opt/homebrew/include/eigen3)
else()
    message(FATAL_ERROR "Eigen not found. Install with: brew install eigen")
endif()

# --- TTL (header-only but provides a CMake target) ---
add_subdirectory(external/TTL)

# --- Source files ---
file(GLOB FEM_SOURCES
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)

# --- Executable ---
add_executable(FEM ${FEM_SOURCES})

# --- TTL include and link ---
target_include_directories(FEM PRIVATE
    ${CMAKE_SOURCE_DIR}/external/TTL/include
)
target_link_libraries(FEM PRIVATE ttl)
