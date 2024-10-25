#include "TFile.h"
#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "SourceFunction1D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"
#include "ChiSquaredTest.hxx"
#include "TH1D.h"

int main()
{
    std::unique_ptr<TFile> otp(TFile::Open("draw1DCF.root","recreate"));

    JJCorrFitter::CorrelationFunction1D func(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func.SetBinning("hCF","",200,2,200);
    func.SetParameters({1,1.},{3.},{});
    func.Evaluate()->Write();

    otp->Close();
}