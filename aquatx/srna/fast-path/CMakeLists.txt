cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 20)
set(EXECUTABLE_OUTPUT_PATH ../)

project(collapser)
set(SOURCES
    collapser.cpp
    collapser.h)

add_executable(collapser ${SOURCES})