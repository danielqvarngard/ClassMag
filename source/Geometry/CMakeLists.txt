add_library(
    ClassMag_Geometry
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/predefLattices.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/reciprocal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/reciprocal.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/disorderedLattices.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/predefLattices.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/lattice.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/euclidean.hpp
)

target_include_directories(
    ClassMag_Geometry
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source/
)