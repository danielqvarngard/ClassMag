add_library(
    ClassMag_MonteCarlo
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/mcFunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/mcFunctions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/VectorModelManager.hpp
)

target_include_directories(
    ClassMag_MonteCarlo
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source
)