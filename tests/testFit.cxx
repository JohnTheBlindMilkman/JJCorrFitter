#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "SourceFunction1D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"
#include "ChiSquaredTest.hxx"
#include "LogLikelihoodTest.hxx"

#include "Minuit2/Minuit2Minimizer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

int main()
{
    using ParType = JJCorrFitter::Fitter::ParType;

    TFile *itp = TFile::Open("/home/jedkol/Downloads/HADES/HADES-CrAP/output/1Dcorr_0_10_cent_forHAL.root");
    std::unique_ptr<TH1> hist(itp->Get<TH1D>("hQinvRatKt2"));
    //std::unique_ptr<TH1> signal(itp->Get<TH1D>("hQinvSignKt2"));
    //std::unique_ptr<TH1> background(itp->Get<TH1D>("hQinvBckgKt2"));

    std::unique_ptr<JJCorrFitter::CorrelationFunction1D> func = std::make_unique<JJCorrFitter::CorrelationFunction1D>(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermPhaseShift>()
        );
    func->SetBinning(hist,2,70);

    JJCorrFitter::Fitter fitter
    (
        std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("GSLSimAn","")),
        std::make_unique<JJCorrFitter::ChiSquaredTest>(std::move(hist),std::move(func))
    );

    fitter.SetParameter(ParType::Generic,"N",1.);
    fitter.SetParameter(ParType::Generic,"Lambda",0.4,0.001,0.,1.);
    fitter.SetParameter(ParType::Source,"Rinv",2.,0.001,1.,6.);

    fitter.SetPrintLevel(1);
    fitter.SetTolerance(1e-9);

    fitter.PrintInfo();

    fitter.Fit();

    std::unique_ptr<TH1> cf = fitter.GetFitFunction();
    hist = fitter.GetDataHistogram();
    hist->GetXaxis()->SetRangeUser(0,200);

    std::unique_ptr<TCanvas> c(new TCanvas("c","",800,800));
    hist->Draw();
    cf->SetLineColor(kRed);
    cf->Draw("c same");

    std::unique_ptr<TFile> otp(TFile::Open("fit1DCF.root","recreate"));
    c->Write();
    hist->Write();
    cf->Write();

    otp->Close();
}