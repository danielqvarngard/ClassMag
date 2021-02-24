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
    ClassMag_Base
    ClassMag_MonteCarlo
)

add_executable(
    testDisorder.exe
    testDisorder.cpp
    ${CMAKE_SOURCE_DIR}/source/Geometry/predefLattices.cpp
    ${CMAKE_SOURCE_DIR}/source/Geometry/include/predefLattices.hpp
    ${CMAKE_SOURCE_DIR}/source/Geometry/include/euclidean.hpp
)

target_include_directories(
    testDisorder.exe
    PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    testDisorder.exe
    PRIVATE
    ClassMag_Geometry
    ClassMag_Base
    ClassMag_MonteCarlo
)

if (TRUE)
add_executable(
    testMC.exe
    testMC.cpp
)
endif()

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


add_executable(
    heatcapTestMC.exe
    heatcapTestMC.cpp
)

target_include_directories(
    heatcapTestMC.exe
    PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    heatcapTestMC.exe
    PRIVATE
    ClassMag_MonteCarlo
    ClassMag_Geometry
)