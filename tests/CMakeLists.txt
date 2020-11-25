add_executable(
    testGeometry.exe
    testGeometry.cpp
)

target_link_libraries(
    testGeometry.exe
    PRIVATE
    ClassMag_Geometry
)