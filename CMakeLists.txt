cmake_minimum_required(VERSION 3.16)
project(BetweennesCentrality)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O2 -g -std=c++17 -fopenmp")

set(SOURCE_EXE src/alex_main.cpp)

add_library(grAlgoLib STATIC
        include/matrix.h
        include/matrix_la.h
        include/matrix_utils.h
        include/betweenness_centrality.h
        include/mmio.h
        src/mmio.cpp)

target_include_directories(grAlgoLib PUBLIC include)

add_executable(grAlgo ${SOURCE_EXE})

target_link_libraries(grAlgo PUBLIC grAlgoLib)

add_executable(my src/main.cpp)

target_link_libraries(my PUBLIC grAlgoLib)