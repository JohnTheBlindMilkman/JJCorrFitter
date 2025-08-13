#include "Fitter.hxx"
#include "CorrelationFunction1D.hxx"
#include "SourceFunction1D.hxx"
#include "DoubleGaussian1D.hxx"
#include "CauchySource1D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "InteractionTermPhaseShift.hxx"
#include "InteractionTermTPI.hxx"
#include "ChiSquaredTest.hxx"
#include "LogLikelihoodTest.hxx"

#include "Minuit2/Minuit2Minimizer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "../../HADES/HADES-CrAP/Externals/Palettes.hxx"

void prepareGraph(const std::unique_ptr<TH1> &hist, Color_t col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(6);
    hist->SetLineColor(col);

    //hist->GetXaxis()->SetLabelSize();
    //hist->GetXaxis()->SetLabelOffset();
    //hist->GetXaxis()->SetTitleSize();
    //hist->GetXaxis()->SetTitleOffset();

    hist->GetXaxis()->SetTitleOffset();
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
}

int main()
{
    using ParType = JJCorrFitter::Fitter::ParType;

    // load data
    std::unique_ptr<TFile> itp(TFile::Open("/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_Purity_MomRes_forHAL.root"));
    std::unique_ptr<TFile> otp(TFile::Open("fitAll1DCF.root","recreate"));

    const double ktBins[] = {0,255,382.5,482.5,575,657.5,737.5,815,895,980,1070,1172.5,1300,1490};
    const double rapBins[] = {0.,0.372,0.48,0.568,0.656,0.748,0.868};
    TH2D radii("radii",";kt;y",14,ktBins,7,rapBins);

    for (const auto &kt : {1,2,3,4,5,6,7,8,9,10,11,12,13,14})
        for (const auto &y : {1,2,3,4,5,6,7})
        {
            try
            {
                std::unique_ptr<TH1> hist(itp->Get<TH1D>(TString::Format("hQinvRatKt%dY%d",kt,y)));

                // create CF object
                std::unique_ptr<JJCorrFitter::CorrelationFunction1D> func = std::make_unique<JJCorrFitter::CorrelationFunction1D>
                (
                    std::make_unique<JJCorrFitter::SourceFunction1D>(),
                    //std::make_unique<JJCorrFitter::InteractionTermTPI>(JJCorrFitter::InteractionTermTPI::SpinState::None)
                    std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
                );
                func->SetBinning(hist,4,100);

                // create fitter object
                JJCorrFitter::Fitter fitter
                (
                    std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad")),
                    std::make_unique<JJCorrFitter::ChiSquaredTest>(std::move(hist),std::move(func))
                );

                // set parameters with limits
                fitter.SetParameter(ParType::Generic,"N",1.,0.001,0.95,1.05);
                fitter.SetParameter(ParType::Generic,"Lambda",1./* ,0.001,0.,1. */);
                fitter.SetParameter(ParType::Source,"Rinv",2.5,0.001,0.5,6.);
                //fitter.SetParameter(ParType::Source,"Rinv2",6.,0.001,0.5,10.);

                // set fit tolerance and ROOT::Minimizer print level
                fitter.SetPrintLevel(0);
                fitter.SetTolerance(1e-9);

                // print setup summary
                //fitter.PrintInfo();

                // perform fit
                fitter.Fit();

                // get CF
                std::unique_ptr<TH1> cf = fitter.GetFitFunction();
                cf->SetName(TString::Format("fitKt%dY%d",kt,y));
                hist = fitter.GetDataHistogram();
                hist->GetXaxis()->SetRangeUser(0,200);
                hist->SetTitle(";k^{*} [MeV/c];C(k^{*})");
                hist->SetName(TString::Format("dataKt%dY%d",kt,y));
                prepareGraph(hist,JJColor::fWutMainColors[1]);

                // make text box with fit params
                std::unique_ptr<TPaveText> tpt(new TPaveText(0.6,0.85,0.98,0.98,"nbr"));
                tpt->SetTextAngle(12);
                tpt->SetFillColor(0);

                std::vector<std::string> parNames{"N","#lambda","R_{inv}"};
                auto params = fitter.GetFitParameterValues();
                auto errors = fitter.GetFitParameterErrors();
                for (std::size_t i = 0 ; i < parNames.size(); ++i)
                {
                    tpt->AddText(TString::Format("%s = %lf +/- %lf",parNames.at(i).c_str(),params.at(i),errors.at(i)));
                }

                radii.SetBinContent(kt,y,params.at(2));
                radii.SetBinError(kt,y,errors.at(2));

                std::unique_ptr<TCanvas> c(new TCanvas(TString::Format("canvKt%dY%d",kt,y),"",800,800));
                c->SetMargin(0.15,0.02,0.15,0.02);
                JJColor::CreatePrimaryWutGradient();

                hist->Draw();
                cf->SetLineColor(JJColor::fWutSecondaryColors[3]);
                cf->Draw("c same");
                tpt->Draw("same");

                c->Write();
                hist->Write();
                cf->Write();
            }
            catch (...)
            {

            }
        }
    radii.Write();
    otp->Save();
    otp->Close();
}