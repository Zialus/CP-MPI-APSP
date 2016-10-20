set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")


set(SOURCE_FILES main.cpp)
add_executable(cp-is-fun ${SOURCE_FILES})