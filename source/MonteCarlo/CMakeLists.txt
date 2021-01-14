add_library(
    ClassMag_MonteCarlo
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/postProcessing.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mcFunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/extendedEnsembles.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/extendedEnsembles.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/postProcessing.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/mcFunctions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/VectorModelManager.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MonteCarloBase.hpp
)

target_include_directories(
    ClassMag_MonteCarlo
    PUBLIC
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(
    ClassMag_MonteCarlo
    PUBLIC
    ClassMag_Base
)