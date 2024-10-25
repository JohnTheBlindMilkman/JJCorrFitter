#include <catch2/catch_template_test_macros.hpp>

#include <memory>

#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"
#include "SourceFunction1D.hxx"
#include "CorrelationFunction1D.hxx"

TEST_CASE("lorem ipsum","[CorrelationFunction1D]")
{
    SECTION("Wavefunction can be calculated using the Schrodinger equation solution")
    {
        JJCorrFitter::CorrelationFunction1D cf(std::make_unique<JJCorrFitter::SourceFunction1D>(),std::make_unique<JJCorrFitter::InteractionTermSchrodinger>());
        cf.SetBinning("hCF","",200,2,200);

        SECTION("CFs can have a custom binning")
        {
            
        }
    }
}