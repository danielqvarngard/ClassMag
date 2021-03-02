if (TRUE)
add_library(
    ClassMag_FileIO
    ${CMAKE_CURRENT_SOURCE_DIR}/StreamManager.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/StreamManager.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/filesystem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/filesystem.hpp
)
endif()

if (FALSE)
add_library(
    ClassMag_FileIO
    ${CMAKE_CURRENT_SOURCE_DIR}/filesystem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/filesystem.hpp
)
endif()