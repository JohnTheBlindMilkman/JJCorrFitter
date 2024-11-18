#include "TFile.h"
#include "CorrelationFunction1D.hxx"
#include "CorrelationFunction3D.hxx"
#include "SourceFunction1D.hxx"
#include "SourceFunction3D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"

int main()
{
    //std::unique_ptr<TFile> otp(TFile::Open("draw1DCF.root","recreate"));

    JJCorrFitter::CorrelationFunction1D func(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermTPI>(JJCorrFitter::InteractionTermTPI::SpinState::None)
        );

    func.SetBinning("hCF","",200,10,200);
    func.SetParameters({1,1.},{2.},{});
    func.Evaluate();

    /* otp->Save();
    otp->Close(); */
}