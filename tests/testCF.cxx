#include "TFile.h"
#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "CorrelationFunction3D.hxx"
#include "SourceFunction1D.hxx"
#include "SourceFunction3D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"
#include "ChiSquaredTest.hxx"
#include "TH3D.h"

int main()
{
    std::unique_ptr<TFile> otp(TFile::Open("draw3DCF.root","recreate"));

    JJCorrFitter::CorrelationFunction3D func(
        std::make_unique<JJCorrFitter::SourceFunction3D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func.SetBinning("hCF","",200,2,200);
    func.SetParameters({1,1.},{2.,2.,2.},{});
    func.Evaluate()->Write();

    otp->Close();
}