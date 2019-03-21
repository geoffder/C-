cmake_minimum_required(VERSION 3.13)
project(retina)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eigen3 )
add_executable(retina main.cpp Cell.cpp Cell.h NetworkModel.cpp NetworkModel.h Stim.cpp Stim.h utils.h eigen_types.h everything.h)
