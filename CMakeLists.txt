##### GENERAL CONFIG ZONE #####

cmake_minimum_required(VERSION 3.9)
project(cp-mpi-apsp)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(MPI REQUIRED)

# Some configs might need "--oversubscribe" flag
message(STATUS "EXTRA_FLAG has been set as: " $ENV{EXTRA_FLAG})
message(STATUS "Max number of procs in this machine is: " ${MPIEXEC_MAX_NUMPROCS})

# Compiler flags to enable more warnings
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif ()

# Compiler flags to enable all warnings
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
endif ()

# Compiler flags to silence/or add some warnings
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w44365")
endif ()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

# Silence warning C4996 on Microsoft Compilers
if (MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
endif ()

##### EXECUTABLE ZONE #####
set(SOURCE_FILES src/main.cpp src/main.hpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/exec")

add_executable(floyd ${SOURCE_FILES})
target_link_libraries(floyd MPI::MPI_CXX)


##### TESTS ZONE #####
enable_testing()
add_test(NAME 1617_12 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1617/input12 ${CMAKE_SOURCE_DIR}/tests/1617/output12)
add_test(NAME 1617_60 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1617/input60 ${CMAKE_SOURCE_DIR}/tests/1617/output60)
add_test(NAME 1617_300 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1617/input300 ${CMAKE_SOURCE_DIR}/tests/1617/output300)

add_test(NAME 1718_5 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1718/input5 ${CMAKE_SOURCE_DIR}/tests/1718/output5)
add_test(NAME 1718_6 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1718/input6 ${CMAKE_SOURCE_DIR}/tests/1718/output6)
add_test(NAME 1718_300 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1718/input300 ${CMAKE_SOURCE_DIR}/tests/1718/output300)
add_test(NAME 1718_600 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1718/input600 ${CMAKE_SOURCE_DIR}/tests/1718/output600)
add_test(NAME 1718_900 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1718/input900 ${CMAKE_SOURCE_DIR}/tests/1718/output900)
add_test(NAME 1718_1200 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/1718/input1200 ${CMAKE_SOURCE_DIR}/tests/1718/output1200)

#add_test(1800 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/something/input1800)
#add_test(2700 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $ENV{EXTRA_FLAG} $<TARGET_FILE:floyd> ${CMAKE_SOURCE_DIR}/tests/something/input2700)
