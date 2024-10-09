cmake_minimum_required(VERSION 3.1...3.29)
project(
    JJCorrFitter
    VERSION 0.3.0
    LANGUAGES CXX
)

# find necessary packages
find_package(ROOT 6.26 CONFIG REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)

# include necessary directories/features
include(GNUInstallDirs)
include(CTest)
include(cmake/Doxygen.cmake)

# adding the configurated file
add_subdirectory(configuredFiles)

# defining options which can be enabled during building
option(JJCORRFITTER_ENABLE_DOXYGEN "Enable doxygen" OFF)
option(JJCORRFITTER_ENABLE_TESTS "Enable tests" ON)

# creates Doxygen documentation with modern and fancy HTML look
if(JJCORRFITTER_ENABLE_DOXYGEN)
    JJCorrFitter_enable_doxygen("awesome-sidebar")
endif()

# enables testing
if(JJCORRFITTER_ENABLE_TESTS)
    add_subdirectory(tests)
    enable_testing()
endif()

# set the directory for binaries and libs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install")

add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor -Wold-style-cast)

# combine all .cxx and .hxx files
set(SRC_FILES
    ${PROJECT_SOURCE_DIR}/src/ChiSquaredTest.cxx
    ${PROJECT_SOURCE_DIR}/src/CorrelationFunction1D.cxx
    ${PROJECT_SOURCE_DIR}/src/Fitter.cxx
    ${PROJECT_SOURCE_DIR}/src/InteractionTermSchrodinger.cxx
    ${PROJECT_SOURCE_DIR}/src/InteractionTermPhaseShift.cxx
    ${PROJECT_SOURCE_DIR}/src/InteractionTermTPI.cxx
    ${PROJECT_SOURCE_DIR}/src/LikelihoodImpl.cxx
    ${PROJECT_SOURCE_DIR}/src/LogLikelihoodTest.cxx
    ${PROJECT_SOURCE_DIR}/src/SourceFunction1D.cxx
)
set(HEADER_FILES
    ${PROJECT_SOURCE_DIR}/include/ChiSquaredTest.hxx
    ${PROJECT_SOURCE_DIR}/include/Config.hxx
    ${PROJECT_SOURCE_DIR}/include/CorrelationFunction1D.hxx
    ${PROJECT_SOURCE_DIR}/include/CorrelationFunctionImpl.hxx
    ${PROJECT_SOURCE_DIR}/include/Fitter.hxx
    ${PROJECT_SOURCE_DIR}/include/InteractionTermImpl.hxx
    ${PROJECT_SOURCE_DIR}/include/InteractionTermSchrodinger.hxx
    ${PROJECT_SOURCE_DIR}/include/InteractionTermPhaseShift.hxx
    ${PROJECT_SOURCE_DIR}/include/InteractionTermTPI.hxx
    ${PROJECT_SOURCE_DIR}/include/LikelihoodImpl.hxx
    ${PROJECT_SOURCE_DIR}/include/LogLikelihoodTest.hxx
    ${PROJECT_SOURCE_DIR}/include/SourceFunction1D.hxx
    ${PROJECT_SOURCE_DIR}/include/SourceFunctionImpl.hxx
)

# ROOT creates his own dictionary for files where ROOT objects are used
ROOT_GENERATE_DICTIONARY(JJCorrFitterDict ${HEADER_FILES} MODULE ${PROJECT_NAME} LINKDEF include/JJCorrFitterLinkDef.h)

# creates our fitter library
add_library(JJCorrFitter SHARED ${SRC_FILES} JJCorrFitterDict.cxx ${HEADER_FILES})
target_compile_features(JJCorrFitter PRIVATE cxx_std_17)


#sets different include destinations depending on if building or installing 
target_include_directories(
    JJCorrFitter
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>" 
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
) 
target_include_directories(
    JJCorrFitter
    PUBLIC
    "$<BUILD_INTERFACE:${Boost_INCLUDE_DIR}>" 
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
) 
target_include_directories(
    JJCorrFitter
    PUBLIC
    "$<BUILD_INTERFACE:/home/jedkol/Downloads/CorAL/include>" 
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
) 
#target_include_directories(JJCorrFitter PUBLIC /home/jedkol/Downloads/CorAL/include)

target_link_libraries(JJCorrFitter 
    ${ROOT_LIBRARIES} 
    GSL::gsl 
    GSL::gslcblas
    ${Boost_LIBRARIES}
)
# I need to hardcode the CorAL unfortunately
target_link_libraries(JJCorrFitter 
    /home/jedkol/Downloads/CorAL/lib/libcoral.a 
    /home/jedkol/Downloads/CorAL/lib/libcoralutils.a 
)

#Only prepares the JJCorrFitterTargets 
install(TARGETS 
    JJCorrFitter 
    EXPORT JJCorrFitterTargets 
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

#Creates and installs the JJCorrFitterTargets at the correct location 
install(EXPORT 
    JJCorrFitterTargets
    FILE JJCorrFitterTargets.cmake
    NAMESPACE JJCorrFitter::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/JJCorrFitter
)

#Install the header files to same destination as that in target_include_directories 
install(FILES ${HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}) 

#Root wants this file in the same folder as the lib 
install(FILES build/libJJCorrFitter_rdict.pcm DESTINATION ${CMAKE_INSTALL_LIBDIR})