if (TRUE)
add_library(
    ClassMag_FileIO
    ${CMAKE_CURRENT_SOURCE_DIR}/StreamManager.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/StreamManager.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/filesystem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/filesystem.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/readLattice.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/readMCO.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/readSpins.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/readCouplings.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/stringExtractions.hpp
)
endif()

target_include_directories(
    ClassMag_FileIO
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source/
)

if (FALSE)
add_library(
    ClassMag_FileIO
    ${CMAKE_CURRENT_SOURCE_DIR}/filesystem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/filesystem.hpp
)
endif()