cmake_minimum_required(VERSION 3.13)
project(retina)
set(CMAKE_CXX_STANDARD 17)


include_directories(C:/eigen3)
#include_directories(C:/boost_1_69_0)
#
#SET(Boost_DEBUG on)
#set(Boost_USE_STATIC_LIBS ON)
#SET(Boost_ARCHITECTURE "-x64")
#SET(BOOST_ROOT "C:/boost_1_69_0")
#SET(BOOST_LIBRARYDIR "C:/boost_1_69_0/stage/lib")
#SET(BOOST_INCLUDEDIR "C:/boost_1_69_0")
#FIND_PACKAGE(Boost COMPONENTS filesystem system REQUIRED)
# seems to fail to find these packages because the built libraries are for VS not MinGW

#include_directories("${Boost_INCLUDE_DIRS}")
#set(BOOST_LIBRARYDIR ${Boost_LIBRARY_DIR})

add_executable(retina
        main.cpp Cell.cpp Cell.h NetworkModel.cpp NetworkModel.h Stim.cpp Stim.h utils.cpp utils.h type_defs.h)

#target_link_libraries(retina ${Boost_LIBRARIES})
#MESSAGE("boost libs: ${Boost_LIBRARIES}")
#target_link_libraries(retina "lboost_filesystem")

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -march=native")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -fopenmp -march=native")
#SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -march=native)

# linking libraries statically to make sure that .exe has all the MinGW stuff it needs
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
# linking libgcc and libstdc++ alone still gave an error about missing libgomp-1.dll when running executable
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ -static")