cmake_minimum_required(VERSION 3.10)
project(ClassMag VERSION 0.0.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -march=native -mno-avx -O3 -fPIC")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mno-avx -Og")
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_subdirectory(${CMAKE_SOURCE_DIR}/source)
add_subdirectory(${CMAKE_SOURCE_DIR}/tests)
add_subdirectory(${CMAKE_SOURCE_DIR}/executables)