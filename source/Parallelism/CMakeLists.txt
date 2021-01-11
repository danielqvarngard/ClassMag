add_library(
    ClassMag_Parallelism
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/messengers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/messengers.hpp
)

target_include_directories(
    ClassMag_Parallelism
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source
)