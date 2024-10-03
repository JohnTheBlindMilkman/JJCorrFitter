#include "TFile.h"
#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "SourceFunction1D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermTPI.hxx"
#include "ChiSquaredTest.hxx"
#include "TH1D.h"

int main()
{
    std::unique_ptr<TFile> otp(TFile::Open("draw1DCF.root","recreate"));

    std::unique_ptr<JJCorrFitter::CorrelationFunctionImpl> func = std::make_unique<JJCorrFitter::CorrelationFunction1D>(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermTPI>(JJCorrFitter::InteractionTermTPI::SpinState::None)
        );

    func->SetBinning("hCF","",200,0.001,0.5);
    func->SetParameters({1,0.6491},{2.6086},{});
    func->Evaluate()->Write();

    otp->Close();
}