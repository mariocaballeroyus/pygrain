cmake_minimum_required(VERSION 3.15)
project(pygrain)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Force all build artifacts into build/ folder
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

# Build core pygrain as static library 
add_library(pygrain_core STATIC 
    cpp/geometry.cpp
    cpp/packing.cpp
    cpp/factory.cpp
    cpp/operations.cpp
)
# Include directories for core library
target_include_directories(pygrain_core PRIVATE cpp)

# Build tests for core library
add_executable(tests_pygrain cpp/pygrain_tests.cpp)
# Include directories and library linking for tests
target_include_directories(tests_pygrain PRIVATE cpp external/catch)
target_link_libraries(tests_pygrain PRIVATE pygrain_core)

# Build Python bindings
add_subdirectory(external/pybind11)
pybind11_add_module(_pygrain3d cpp/pygrain_bindings.cpp)
# Include directories and library linking for bindings
target_include_directories(_pygrain3d PRIVATE cpp)
target_link_libraries(_pygrain3d PRIVATE pygrain_core)

# Install the extension module into the pygrain3d package
install(TARGETS _pygrain3d DESTINATION pygrain3d)