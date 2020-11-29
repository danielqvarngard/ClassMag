add_executable(
    testGeometry.exe
    testGeometry.cpp
    ${CMAKE_SOURCE_DIR}/source/Geometry/predefLattices.cpp
    ${CMAKE_SOURCE_DIR}/source/Geometry/include/predefLattices.hpp
    ${CMAKE_SOURCE_DIR}/source/Geometry/include/euclidean.hpp
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