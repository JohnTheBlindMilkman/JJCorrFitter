#include "TFile.h"
#include "TCanvas.h"
#include "CorrelationFunction1D.hxx"
#include "CorrelationFunction3D.hxx"
#include "SourceFunction1D.hxx"
#include "DoubleGaussian1D.hxx"
#include "CauchySource1D.hxx"
#include "SourceFunction3D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"

int main()
{
    std::unique_ptr<TFile> otp(TFile::Open("draw1DCF.root","recreate"));

    JJCorrFitter::CorrelationFunction1D func1(
        std::make_unique<JJCorrFitter::CauchySource1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func1.SetBinning("hCF1","",200,2,200);
    func1.SetParameters({1,1.},{3.},{});
    func1.Evaluate()->Write();

    JJCorrFitter::CorrelationFunction1D func2(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func2.SetBinning("hCF2","",200,2,200);
    func2.SetParameters({1,1.},{3.},{});
    func2.Evaluate()->Write();

    JJCorrFitter::CorrelationFunction1D func3(
        std::make_unique<JJCorrFitter::DoubleGaussian1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func3.SetBinning("hCF3","",200,2,200);
    func3.SetParameters({1,1.},{2.,6.},{});
    func3.Evaluate()->Write();

    otp->Save();
    otp->Close();
}