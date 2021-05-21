add_executable(
    tsaiDisorder.mpix
    tsaiDisorder.cpp
)

target_include_directories(
    tsaiDisorder.mpix
    PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    tsaiDisorder.mpix
    PRIVATE
    ClassMag_MonteCarlo
)

add_executable(
    isingpt.mpix
    isingpt.cpp
)

target_include_directories(
    isingpt.mpix
    PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    isingpt.mpix
    PRIVATE
    ClassMag_ParallelTempering
)