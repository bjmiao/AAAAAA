cmake_minimum_required(VERSION 3.9)
project(GraphAndSet)

set(CMAKE_CXX_STANDARD 11)

add_executable(GraphAndSet liuzhijun.cpp AdjMatrixGraph.h AdjListGraph.h SimpleSCCOnVectorGraph.h SimpleSCCOnMatrixGraph.h DisjointSet.h)