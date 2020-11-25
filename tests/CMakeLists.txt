add_executable(
    testGeometry.exe
    testGeometry.cpp
    ${CMAKE_SOURCE_DIR}/source/Geometry/include/qartesian.hpp
)

target_include_directories(
    testGeometry.exe
    PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    testGeometry.exe
    PRIVATE
    ClassMag_Geometry
)