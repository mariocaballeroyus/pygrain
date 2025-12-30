cmake_minimum_required(VERSION 3.10)
project(pygrain)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB PYGRAIN_SOURCES cpp/src/*.cpp)
file(GLOB PYGRAIN_TESTS cpp/test/*.cpp)
file(GLOB PYGRAIN_BINDS cpp/bind/*.cpp)

add_library(pygrain_lib STATIC ${PYGRAIN_SOURCES})
target_include_directories(pygrain_lib PRIVATE cpp/inc)

add_executable(tests_pygrain ${PYGRAIN_TESTS})
target_include_directories(tests_pygrain PRIVATE cpp/inc external/catch)
target_link_libraries(tests_pygrain PRIVATE pygrain_lib)

add_executable(main_pygrain main.cpp)
target_include_directories(main_pygrain PRIVATE cpp/inc)
target_link_libraries(main_pygrain PRIVATE pygrain_lib)

add_subdirectory(external/pybind11)
pybind11_add_module(_pygrain ${PYGRAIN_BINDS})
target_include_directories(_pygrain PRIVATE cpp/inc)
target_link_libraries(_pygrain PRIVATE pygrain_lib)