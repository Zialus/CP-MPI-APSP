cmake_minimum_required(VERSION 2.8)
project(cp-mpi-apsp)

#set(CMAKE_VERBOSE_MAKEFILE true)

# This project requires MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Files with the source code
set(SOURCE_FILES main.cpp)
add_executable(floyd ${SOURCE_FILES})

target_link_libraries(floyd ${MPI_LIBRARIES})

# Tests:
enable_testing()
add_test(6 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/input6.txt)
add_test(60 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/input60.txt)
add_test(300 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/input300.txt)
add_test(900 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/input900.txt)
# add_test(1800 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/input1800.txt)
