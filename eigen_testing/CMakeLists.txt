cmake_minimum_required(VERSION 3.13)
project(eigen_testing)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eigen3 )
add_executable(eigen_testing main.cpp)