add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/FileIO)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Geometry)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Base)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/MonteCarlo)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Parallelism)

if(TRUE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Environments)
endif()