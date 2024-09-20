#include "TFile.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "SourceFunction1D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "ChiSquaredTest.hxx"
#include "TH1D.h"

int main()
{
    using ParType = JJCorrFitter::Fitter::ParType;

    TFile *itp = TFile::Open("/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_DR_forHAL.root");
    std::unique_ptr<TH1D> hist(itp->Get<TH1D>("hQinvDRKt2"));
    // hist->SetDirectory(nullptr) ?

    std::unique_ptr<JJCorrFitter::CorrelationFunction1D> func = std::make_unique<JJCorrFitter::CorrelationFunction1D>(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    JJCorrFitter::Fitter fitter
    (
        std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad")),
        std::make_unique<JJCorrFitter::ChiSquaredTest>(std::move(hist),std::move(func))
    );
    hist = nullptr; // I moved it so hist is undefined, hence the new assignment

    fitter.SetParameter(ParType::Generic,"N",1.,0.1,0,10e4);
    fitter.SetParameter(ParType::Generic,"Lambda",1.);
    fitter.SetParameter(ParType::Source,"Rinv",2.,0.1,1,6);

    fitter.Fit();

    //TFile otp("fit1DCF.root","recreate");

    // fitter.GetFitFunction()->Write();

    //otp.Close();
}