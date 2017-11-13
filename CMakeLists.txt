cmake_minimum_required(VERSION 3.0)
project(cp-mpi-apsp)

# This project requires MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

message(STATUS "EXTRA_FLAG has been set as: " $ENV{EXTRA_FLAG})
message(STATUS "Max number of procs in this machine is: " ${MPIEXEC_MAX_NUMPROCS})

# Set executable name and list the files with source code
set(EXEC floyd)
set(SOURCE_FILES main.cpp)
add_executable(${EXEC} ${SOURCE_FILES})

# Link MPI libraries to the executable
target_link_libraries(${EXEC} ${MPI_LIBRARIES})

# Tests:
enable_testing()
add_test(1617_12 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1617/input12 ../tests/1617/output12)
add_test(1617_60 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1617/input60 ../tests/1617/output60)
add_test(1617_300 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1617/input300 ../tests/1617/output300)

add_test(1718_5 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1718/input5 ../tests/1718/output5)
add_test(1718_6 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1718/input6 ../tests/1718/output6)
add_test(1718_300 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1718/input300 ../tests/1718/output300)
add_test(1718_600 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1718/input600 ../tests/1718/output600)
add_test(1718_900 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1718/input900 ../tests/1718/output900)
add_test(1718_1200 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/1718/input1200 ../tests/1718/output1200)


# add_test(1800 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/something/input1800)
# add_test(2700 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} ./${EXEC} ../tests/something/input2700)
