add_library(
    ClassMag_ParallelTempering
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/print_pt_return.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/print_pt_return.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PT_Return.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ptProcs.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ptProcs.hpp
)

target_include_directories(
    ClassMag_ParallelTempering
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source/
)

target_link_libraries(
    ClassMag_ParallelTempering
    PUBLIC
    ClassMag_MonteCarlo
    ClassMag_Parallelism
    ClassMag_FileIO
)