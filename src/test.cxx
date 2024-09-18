#include "TFile.h"
#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "ChiSquaredTest.hxx"
#include "TH1D.h"

int main()
{
    TFile *itp = TFile::Open("/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_DR_forHAL.root");
    std::unique_ptr<TH1D> hist(itp->Get<TH1D>("hQinvDRKt2"));
    // hist->SetDirectory(nullptr) ?

    JJCorrFitter::Fitter fitter
    (
        std::move(hist),
        std::make_unique<JJCorrFitter::CorrelationFunction1D>(),
        std::make_unique<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad")),
        std::make_unique<JJCorrFitter::ChiSquaredTest>()
    );
    hist = nullptr; // I moved it so hist is undefined, hence the new assignment

    fitter.SetParameter(0,"Lambda",1.);
    fitter.SetParameter(1,"N",1.,0.1,0,10e4);
    fitter.SetParameter(2,"Rinv",2.,0.1,1,6);

    fitter.Fit();

    TFile *otp = TFile::Open("draw1DCF.root","recreate");

    // fitter.GetFitFunction()->Write();
    // hist->Write();

    otp->Close();
}