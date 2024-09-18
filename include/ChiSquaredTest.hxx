/**
 * @file ChiSquaredTest.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Goodness-of-fit function for the minimizer whihc utilises the Chi-squared test
 * @version 0.2.0
 * @date 2024-09-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef ChiSquaredTest_hxx
    #define ChiSquaredTest_hxx

    #include "LikelihoodImpl.hxx"

    namespace JJCorrFitter
    {
        class ChiSquaredTest : public LikelihoodImpl
        {
            private:
                /* data */
            public:
                ChiSquaredTest(/* args */);
                ~ChiSquaredTest() = default;
        };

    } // namespace JJCorrFitter
    
#endif