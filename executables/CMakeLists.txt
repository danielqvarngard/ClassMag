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

if(0)
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

add_executable(
    defaultPT.mpix
    defaultPT.cpp
)

target_include_directories(
    defaultPT.mpix
    PRIVATE
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    defaultPT.mpix
    PRIVATE
    ClassMag_ParallelTempering
)
endif()

if(MPI_CXX_FOUND)
    add_executable(
        approximantPaper.mpix
        ${CMAKE_CURRENT_SOURCE_DIR}/approximantPaper.cpp
    )

    target_link_libraries(
        approximantPaper.mpix
        PUBLIC
        MPI::MPI_CXX
        ClassMag_MonteCarlo
        ClassMag_FileIO
        ClassMag_Parallelism
        ClassMag_ParallelTempering
    )

    add_executable(
        approximantPaperXY.mpix
        ${CMAKE_CURRENT_SOURCE_DIR}/approximantPaperXY.cpp
    )

    target_link_libraries(
        approximantPaperXY.mpix
        PUBLIC
        MPI::MPI_CXX
        ClassMag_MonteCarlo
        ClassMag_FileIO
        ClassMag_Parallelism
        ClassMag_ParallelTempering
    )

    add_executable(
        mpitest.mpix
        ${CMAKE_CURRENT_SOURCE_DIR}/mpitest.cpp
    )
endif()

add_executable(
    freezing_temperature.exe
    ${CMAKE_CURRENT_SOURCE_DIR}/freezing_temperature.cpp
)

target_link_libraries(
    freezing_temperature.exe
    PUBLIC
    ClassMag_Base
    ClassMag_MonteCarlo
    ClassMag_FileIO
)

target_include_directories(
    freezing_temperature.exe
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source
)


add_executable(
    cw_temperature.exe
    ${CMAKE_CURRENT_SOURCE_DIR}/cw_temperature.cpp
)

target_link_libraries(
    cw_temperature.exe
    PUBLIC
    ClassMag_Base
    ClassMag_MonteCarlo
    ClassMag_FileIO
    ClassMag_Geometry
)

target_include_directories(
    cw_temperature.exe
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source
)