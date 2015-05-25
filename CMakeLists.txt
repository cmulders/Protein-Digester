cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

MESSAGE( STATUS "Building in ${BUILDENV} Environment" )

project(CPPDigester)

add_definitions("-std=c++11")

add_subdirectory(vendor)
include_directories(vendor)

add_subdirectory(src)
