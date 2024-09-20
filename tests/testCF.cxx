#include "TFile.h"
#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "SourceFunction1D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "ChiSquaredTest.hxx"
#include "TH1D.h"

int main()
{
    TFile otp("draw1DCF.root","recreate");

    std::unique_ptr<JJCorrFitter::CorrelationFunctionImpl> func = std::make_unique<JJCorrFitter::CorrelationFunction1D>(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func->SetParameters({1,1},{2},{});
    func->Evaluate()->Write();

    otp.Close();
}