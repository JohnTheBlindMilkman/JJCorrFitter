#include "Fitter.hxx"
#include "CorrelationFunction3D.hxx"
#include "SourceFunction3D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"
#include "ChiSquaredTest.hxx"

#include "Minuit2/Minuit2Minimizer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "../../HADES/HADES-CrAP/Externals/Palettes.hxx"

#include <chrono>

void prepareGraph(const std::unique_ptr<TH1> &hist, Color_t col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(6);
    hist->SetLineColor(col);

    hist->GetXaxis()->SetTitleOffset();
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
    hist->GetZaxis()->SetTitleSize(0.06);
    hist->GetZaxis()->SetLabelSize(0.06);
    hist->GetZaxis()->SetNdivisions(506);
}

int main()
{
    using ParType = JJCorrFitter::Fitter::ParType;

    // load data
    std::unique_ptr<TFile> itp(TFile::Open("/home/jedkol/Downloads/JJCorrFitter/macros/input/draw3DCF_sim_errs_pts.root"));
    std::unique_ptr<TH1> hist(itp->Get<TH3D>("hQoslRatInteg_gauss"));

    // create CF object
    std::unique_ptr<JJCorrFitter::CorrelationFunction3D> func = std::make_unique<JJCorrFitter::CorrelationFunction3D>
    (
        std::make_unique<JJCorrFitter::SourceFunction3D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
    );
    func->SetBinning(hist,12,80);

    JJCorrFitter::CorrelationFunction3D func2
    (
        std::make_unique<JJCorrFitter::SourceFunction3D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
    );
    func2.SetBinning("fitFullInteg",hist->GetTitle(),hist->GetNbinsX() * 4,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

    // create fitter object
    JJCorrFitter::Fitter fitter
    (
        std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad")),
        // std::make_unique<JJCorrFitter::LogLikelihoodTest>(std::move(hist), std::move(signal), std::move(background), std::move(func))
        std::make_unique<JJCorrFitter::ChiSquaredTest>(std::move(hist), std::move(func))
    );

    // set parameters with limits
    fitter.SetParameter(ParType::Generic,"N",1.,0.001,0.95,1.05);
    fitter.SetParameter(ParType::Generic,"Lambda",1./* ,0.001,0.,1. */);
    fitter.SetParameter(ParType::Source,"Rout",2.5,0.001,0.5,6.);
    fitter.SetParameter(ParType::Source,"Rside",2.5,0.001,0.5,6.);
    fitter.SetParameter(ParType::Source,"Rlong",2.5,0.001,0.5,6.);

    // set fit tolerance and ROOT::Minimizer print level
    fitter.SetPrintLevel(0);
    fitter.SetTolerance(1e-9);

    // print setup summary
    fitter.PrintInfo();

    auto start = std::chrono::high_resolution_clock::now();
    // perform fit
    fitter.Fit();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Fitting time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "\n";

    // get CF
    std::unique_ptr<TH1> cf = fitter.GetFitFunction();
    cf->SetLineColor(JJColor::fWutSecondaryColors[3]);
    hist = fitter.GetDataHistogram();
    hist->GetXaxis()->SetRangeUser(0,200);
    hist->SetTitle(";k^{*}_{out} [MeV/c];k^{*}_{side} [MeV/c];k^{*}_{long} [MeV/c];C(k^{*})");
    prepareGraph(hist,JJColor::fWutMainColors[1]);

    auto params = fitter.GetFitParameterValues();
    auto errors = fitter.GetFitParameterErrors();

    func2.SetParameters({params[0],params[1]},{params[2],params[3],params[4]},{});
    auto cfFull = func2.EvaluateAtPlanes(4);
    cfFull->SetLineColor(JJColor::fWutSecondaryColors[3]);
    cfFull->SetLineStyle(kDashed);

    std::unique_ptr<TFile> otp(TFile::Open("fit3DCF.root","recreate"));
    hist->Write();
    cf->Write();
    cfFull->Write();

    otp->Save();
    otp->Close();
}