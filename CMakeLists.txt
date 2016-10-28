cmake_minimum_required(VERSION 3.2)

find_package(MPI REQUIRED)

set(SOURCE_FILES main.cpp)
add_executable(cp-is-fun ${SOURCE_FILES})