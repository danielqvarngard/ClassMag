add_library(
    ClassMag_Geometry
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/predefLattices.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/predefLattices.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/lattice.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/euclidean.hpp
)