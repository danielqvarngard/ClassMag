add_library(
    ClassMag_Base
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/orderParameter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/couplingLookup.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rkky.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/rkky.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/nearestNeighbor.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/couplingLookup.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/orderParameter.hpp
)

target_link_libraries(
    ClassMag_Base
    PUBLIC
    ClassMag_Geometry
)

target_include_directories(
    ClassMag_Base
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source
    )