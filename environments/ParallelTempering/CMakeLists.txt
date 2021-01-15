add_library(
    ClassMag_ParallelTempering
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/ptProcs.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ptProcs.hpp
)

target_include_directories(
    ClassMag_ParallelTempering
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source/
)