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

add_executable(
    testMC.exe
    testMC.cpp
)

target_include_directories(
    testMC.exe
    PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    testMC.exe
    PRIVATE
    ClassMag_MonteCarlo
    ClassMag_Geometry
)