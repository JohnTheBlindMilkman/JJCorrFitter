#include "TFile.h"
#include "Fitter1D.hxx"

int main()
{
    TFile *itp = TFile::Open("/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_DR_forHAL.root");
    std::unique_ptr<TH1D> hist(itp->Get<TH1D>("hQinvDRKt2"));
    // hist->SetDirectory(nullptr) ?

    JJCorrFitter::Fitter1D fitter(std::move(hist));
    hist = nullptr;

    fitter.SetParameter(JJCorrFitter::Fitter1D::Parameter::Lambda,1.);
    fitter.SetParameter(JJCorrFitter::Fitter1D::Parameter::Norm,1.,0,10e4);
    fitter.SetParameter(JJCorrFitter::Fitter1D::Parameter::Radius,2.,1,6);

    fitter.Fit();

    TFile *otp = TFile::Open("draw1DCF.root","recreate");

    // fitter.GetFitFunction()->Write();
    // hist->Write();

    otp->Close();
}