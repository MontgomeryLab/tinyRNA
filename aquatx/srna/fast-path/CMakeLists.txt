cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 20)
include(GNUInstallDirs)

project(collapser)
set(SOURCES
        collapser.cpp
        collapser.h
)

set(EXECUTABLE_OUTPUT_PATH ../)
add_executable(collapser ${SOURCES})