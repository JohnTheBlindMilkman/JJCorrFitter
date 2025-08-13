if(JJCORRFITTER_ENABLE_TESTS)

    add_executable(draw1DCF testCF.cxx)
    target_link_libraries(draw1DCF PRIVATE JJCorrFitter)
    target_compile_features(draw1DCF  PRIVATE cxx_std_17)

    add_executable(calcMomRes1D testMomRes.cxx)
    target_link_libraries(calcMomRes1D PRIVATE JJCorrFitter)
    target_compile_features(calcMomRes1D  PRIVATE cxx_std_17)

    add_executable(calcMomRes1DSMASH testMomResSMASH.cxx)
    target_link_libraries(calcMomRes1DSMASH PRIVATE JJCorrFitter)
    target_compile_features(calcMomRes1DSMASH  PRIVATE cxx_std_17)

    add_executable(fit1DCF testFit.cxx)
    target_link_libraries(fit1DCF PRIVATE JJCorrFitter)
    target_compile_features(fit1DCF  PRIVATE cxx_std_17)
    add_executable(fitAll1DCF fitAll.cxx)
    target_link_libraries(fitAll1DCF PRIVATE JJCorrFitter)
    target_compile_features(fitAll1DCF  PRIVATE cxx_std_17)

    add_executable(tests tests.cxx)
    target_link_libraries(tests PRIVATE JJCorrFitter PRIVATE Catch2::Catch2WithMain)
    target_compile_features(tests PRIVATE cxx_std_17)

    add_test(NAME draw1DCF COMMAND draw1DCF)
    add_test(NAME fit1DCF COMMAND fit1DCF)

    catch_discover_tests(tests)

endif()