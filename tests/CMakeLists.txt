if(JJCORRFITTER_ENABLE_TESTS)

    add_executable(draw1DCF testCF.cxx)
    target_link_libraries(draw1DCF PRIVATE JJCorrFitter)
    target_compile_features(draw1DCF  PRIVATE cxx_std_17)

    add_executable(fit1DCF testFit.cxx)
    target_link_libraries(fit1DCF PRIVATE JJCorrFitter)
    target_compile_features(fit1DCF  PRIVATE cxx_std_17)

    add_executable(tests tests.cxx)
    target_link_libraries(tests PRIVATE JJCorrFitter PRIVATE Catch2::Catch2WithMain)
    target_compile_features(tests PRIVATE cxx_std_17)

    add_test(NAME draw1DCF COMMAND draw1DCF)
    add_test(NAME fit1DCF COMMAND fit1DCF)

    catch_discover_tests(tests)

endif()