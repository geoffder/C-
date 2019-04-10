cmake_minimum_required(VERSION 3.13)
project(retina)
set(CMAKE_CXX_STANDARD 17)


include_directories(C:/eigen3)
#include_directories(C:/eigen3-hdf5-master)
#include_directories(C:/HDF_Group/HDF5/1.10.5/include)

add_executable(retina
        main.cpp Cell.cpp Cell.h NetworkModel.cpp NetworkModel.h Stim.cpp Stim.h utils.cpp utils.h type_defs.h)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -march=native")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -fopenmp -march=native")
#SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -march=native)

# linking libraries statically to make sure that .exe has all the MinGW stuff it needs
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
# linking libgcc and libstdc++ alone still gave an error about missing libgomp-1.dll when running executable
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ -static")