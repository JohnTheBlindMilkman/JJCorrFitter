cmake_minimum_required(VERSION 3.1...3.29)
project(
    JJCorrFitter
    VERSION 0.1.0
    LANGUAGES CXX
)

# find necessary packages
find_package(ROOT 6.26 CONFIG REQUIRED)

# include necessary directories/features
include(GNUInstallDirs)
include(CTest)
include(cmake/Doxygen.cmake)

# adding the configurated file
add_subdirectory(configuredFiles)

# defining options which can be enabled during building
option(JJCORRFITTER_ENABLE_DOXYGEN "Enable doxygen" OFF)
option(JJCORRFITTER_ENABLE_TESTS "Enable tests" OFF)

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
set( sources
    ${PROJECT_SOURCE_DIR}/src/
)
set( includes
    ${PROJECT_SOURCE_DIR}/include/Config.hxx
)

# ROOT creates his own dictionary for files where ROOT objects are used
ROOT_GENERATE_DICTIONARY(JJCorrFitterDict ${includes} MODULE ${PROJECT_NAME} LINKDEF include/JJCorrFitterLinkDef.h)

# creates our fitter library
add_library(JJCorrFitter SHARED ${sources} JJCorrFitterDict.cxx ${includes})
target_compile_features(JJCorrFitter PRIVATE cxx_std_17)

#sets different include destinations depending on if building or installing 
target_include_directories(
    JJCorrFitter
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>" 
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
) 
target_link_libraries(JJCorrFitter ${ROOT_LIBRARIES})

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
install(FILES ${includes} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}) 

#Root wants this file in the same folder as the lib 
install(FILES build/libJJCorrFitter_rdict.pcm DESTINATION ${CMAKE_INSTALL_LIBDIR})