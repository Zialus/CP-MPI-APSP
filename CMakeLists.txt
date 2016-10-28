cmake_minimum_required(VERSION 2.8)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


set(SOURCE_FILES main.cpp)
add_executable(cp-is-fun ${SOURCE_FILES})
target_link_libraries(cp-is-fun ${MPI_C_LIBRARIES})
