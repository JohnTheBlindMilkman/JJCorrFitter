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

#include <vector>
#include <memory>

int main()
{
    // load data
    std::unique_ptr<TFile> itp(TFile::Open("/home/jedkol/Downloads/HADES/mstefan/chuan.root"));
    std::unique_ptr<TH1> hist(itp->Get<TH1D>("hist"));

    std::unique_ptr<TFile> otp(TFile::Open("draw1DCF.root","recreate"));
    hist->Write();

    JJCorrFitter::CorrelationFunction1D func1(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    std::unique_ptr<TCanvas> c(new TCanvas("c","",800,800));
    c->SetMargin(0.15,0.02,0.15,0.02);

    func1.SetBinning(hist,0,200);
    func1.SetParameters({1.02078,1.},{3.2026},{});
    func1.Evaluate()->Write("hQinvRatInteg_fitGauss");

    JJCorrFitter::CorrelationFunction1D func2(
        std::make_unique<JJCorrFitter::CauchySource1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func2.SetBinning(hist,0,200);
    func2.SetParameters({0.996768,1.},{2.94521},{});
    func2.Evaluate()->Write("hQinvRatInteg_fitCauchy");

    JJCorrFitter::CorrelationFunction1D func3(
        std::make_unique<JJCorrFitter::DoubleGaussian1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func3.SetBinning(hist,0,200);
    func3.SetParameters({1.00859,1.},{2.35467,5.05293},{});
    func3.Evaluate()->Write("hQinvRatInteg_fitDoubleGauss");

    otp->Save();
    otp->Close();
}