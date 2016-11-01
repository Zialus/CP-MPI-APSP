cmake_minimum_required(VERSION 2.8)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


set(SOURCE_FILES main.cpp)
add_executable(floyd ${SOURCE_FILES})
target_link_libraries(floyd ${MPI_LIBRARIES})

# Add a test:
enable_testing()
add_test(6 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/Input6.txt)
add_test(60 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/Input60.txt)
add_test(300 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/Input300.txt)
add_test(900 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 floyd ../Tests/Input900.txt)
