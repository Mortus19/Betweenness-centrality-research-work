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
        include/Graph.h
        include/StatisticTools.h
        src/StatisticTools.cpp
        src/mmio.cpp)

target_include_directories(grAlgoLib PUBLIC include)

add_executable(grAlgo ${SOURCE_EXE})
target_link_libraries(grAlgo PUBLIC grAlgoLib)

add_executable(TBC-base src/TBC-base.cpp)
target_link_libraries(TBC-base PUBLIC grAlgoLib)

add_executable(TBC-optimized-only-vertex src/TBC-optimized-only-vertex.cpp)
target_link_libraries(TBC-optimized-only-vertex PUBLIC grAlgoLib)

add_executable(TBC-optimized-all src/TBC-optimized-all.cpp)
target_link_libraries(TBC-optimized-all PUBLIC grAlgoLib)

add_executable(TBC-uniform-source-sampling src/TBC-uniform-source-sampling.cpp)
target_link_libraries(TBC-uniform-source-sampling PUBLIC grAlgoLib)

add_executable(TBC-probability-source-sampling src/TBC-probability-source-sampling.cpp)
target_link_libraries(TBC-probability-source-sampling PUBLIC grAlgoLib)

add_executable(TBC-uniform-path-sampling src/TBC-uniform-path-sampling.cpp)
target_link_libraries(TBC-uniform-path-sampling PUBLIC grAlgoLib)

add_executable(TBC-uniform-path-sampling-optimazed src/TBC-uniform-path-sampling-optimized.cpp)
target_link_libraries(TBC-uniform-path-sampling-optimazed PUBLIC grAlgoLib)

add_executable(TBC-uniform-path-sampling-optimized-vertex src/TBC-uniform-path-sampling-optimized-vertex.cpp)
target_link_libraries(TBC-uniform-path-sampling-optimized-vertex PUBLIC grAlgoLib)
