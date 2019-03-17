# cmake_minimum_required(VERSION <specify CMake version here>)
project(retina)

set(CMAKE_CXX_STANDARD 17)
include_directories(${../eigen-eigen-323c052e1731}/src )
add_executable(retina main.cpp)
