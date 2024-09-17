if(JJCORRFITTER_ENABLE_TESTS)

    add_executable(draw1DCF draw1DCorrFunc.cxx)
    target_link_libraries(draw1DCF PRIVATE JJCorrFitter)
    target_compile_features(draw1DCF  PRIVATE cxx_std_17)

    #add_executable(parseOutFile parseOutFile.cxx)
    #target_link_libraries(parseOutFile PRIVATE base)
    #target_compile_features(parseOutFile  PRIVATE cxx_std_20)

    add_test(NAME draw1DCF COMMAND draw1DCF)
    #add_test(NAME parseOutFile COMMAND parseOutFile "${CMAKE_SOURCE_DIR}/tests/slurm-12345_101.out")

endif()