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