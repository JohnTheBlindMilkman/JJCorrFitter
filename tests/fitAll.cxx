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

    constexpr std::array<int,12> ktBinNums{1,2,3,4,5,6,7,8,9,10,11,12};
    constexpr std::array<int,9> rapBinNums{1,2,3,4,5,6,7,8,9};

    // load data
    std::unique_ptr<TFile> itp(TFile::Open("/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_Purity_MomRes_forHAL.root"));
    std::unique_ptr<TFile> otp(TFile::Open("fitAll1DCF.root","recreate"));

    const double ktBins[ktBinNums.size() + 1] = {300,450,600,750,900,1050,1200,1350,1500,1650,1800,1950,2100};
    const double rapBins[rapBinNums.size() + 1] = {0.09,0.19,0.29,0.39,0.49,0.59,0.69,0.79,0.89,0.99};
    TH2D chisqr = TH2D("chisqrKtY",";k_{T};y^{1,2};#chi^{2}",ktBinNums.size(),ktBins,rapBinNums.size(),rapBins);

    std::vector<TH1D> radii;
    for (const auto &i : rapBinNums)
    {
        radii.emplace_back(TH1D(TString::Format("radiiKtY_%d",i),";k_{T};R_{inv}",ktBinNums.size(),ktBins));
    }

    std::unique_ptr<TH1> hist(itp->Get<TH1D>("hQinvRatInteg"));
    // create CF object
    std::unique_ptr<JJCorrFitter::CorrelationFunction1D> func = std::make_unique<JJCorrFitter::CorrelationFunction1D>
    (
        std::make_unique<JJCorrFitter::CauchySource1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
    );
    func->SetBinning(hist,12,80);

    JJCorrFitter::CorrelationFunction1D func2
    (
        std::make_unique<JJCorrFitter::CauchySource1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
    );
    func2.SetBinning(hist,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

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
    fitter.PrintInfo();

    TH1D radiiKt = TH1D("radiiKt",";k_{T};R_{inv}",ktBinNums.size(),ktBins);
    TH1D chisqrKt = TH1D("chisqrKt",";k_{T};#chi^{2}",ktBinNums.size(),ktBins);

    for (const auto &kt : ktBinNums)
    {
        std::cout << "kt: " << kt << "\n";
        hist = std::unique_ptr<TH1>(itp->Get<TH1D>(TString::Format("hQinvRatKt%d",kt)));
        fitter.SetDataHistogram(std::move(hist));            

        // perform fit
        fitter.Fit();

        // get CF
        std::unique_ptr<TH1> cf = fitter.GetFitFunction();
        cf->SetName(TString::Format("fitKt%d",kt));
        hist = fitter.GetDataHistogram();
        hist->GetXaxis()->SetRangeUser(0,200);
        hist->SetTitle(";k^{*} [MeV/c];C(k^{*})");
        hist->SetName(TString::Format("dataKt%d",kt));
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

        func2.SetParameters({params[0],params[1]},{params[2]},{});
        auto cfFull = func2.Evaluate();
        cfFull->SetName(TString::Format("fitFullKt%d",kt));
        cfFull->SetLineColor(JJColor::fWutSecondaryColors[3]);
        cfFull->SetLineStyle(kDashed);

        if (fitter.GetStatus() == 0 || fitter.GetStatus() == 1)
        {
            radiiKt.SetBinContent(kt,params.at(2));
            radiiKt.SetBinError(kt,errors.at(2));
            chisqrKt.SetBinContent(kt,fitter.GetMinValue());
        }
        else
        {
            radiiKt.SetBinContent(kt,0.);
            radiiKt.SetBinError(kt,0.);
            chisqrKt.SetBinContent(kt,0.);
        }

        std::unique_ptr<TCanvas> c(new TCanvas(TString::Format("canvKt%d",kt),"",800,800));
        c->SetMargin(0.15,0.02,0.15,0.02);
        JJColor::CreatePrimaryWutGradient();

        hist->Draw();
        cf->SetLineColor(JJColor::fWutSecondaryColors[3]);
        cf->Draw("c same");
        cfFull->Draw("c same");
        tpt->Draw("same");

        c->Write();
        hist->Write();
        cf->Write();
        cfFull->Write();
    }
    radiiKt.Write();
    chisqrKt.Write();

    TH1D radiiY = TH1D("radiiY",";y^{1,2};R_{inv}",rapBinNums.size(),rapBins);
    TH1D chisqrY = TH1D("chisqrY",";y^{1,2};#chi^{2}",rapBinNums.size(),rapBins);

    for (const auto &y : rapBinNums)
    {
        std::cout << "y: " << y << "\n";
        hist = std::unique_ptr<TH1>(itp->Get<TH1D>(TString::Format("hQinvRatY%d",y)));
        fitter.SetDataHistogram(std::move(hist));            

        // perform fit
        fitter.Fit();

        // get CF
        std::unique_ptr<TH1> cf = fitter.GetFitFunction();
        cf->SetName(TString::Format("fitY%d",y));
        hist = fitter.GetDataHistogram();
        hist->GetXaxis()->SetRangeUser(0,200);
        hist->SetTitle(";k^{*} [MeV/c];C(k^{*})");
        hist->SetName(TString::Format("dataY%d",y));
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

        func2.SetParameters({params[0],params[1]},{params[2]},{});
        auto cfFull = func2.Evaluate();
        cfFull->SetName(TString::Format("fitFullY%d",y));
        cfFull->SetLineColor(JJColor::fWutSecondaryColors[3]);
        cfFull->SetLineStyle(kDashed);

        if (fitter.GetStatus() == 0 || fitter.GetStatus() == 1)
        {
            radiiY.SetBinContent(y,params.at(2));
            radiiY.SetBinError(y,errors.at(2));
            chisqrY.SetBinContent(y,fitter.GetMinValue());
        }
        else
        {
            radiiY.SetBinContent(y,0.);
            radiiY.SetBinError(y,0.);
            chisqrY.SetBinContent(y,0.);
        }

        std::unique_ptr<TCanvas> c(new TCanvas(TString::Format("canvY%d",y),"",800,800));
        c->SetMargin(0.15,0.02,0.15,0.02);
        JJColor::CreatePrimaryWutGradient();

        hist->Draw();
        cf->SetLineColor(JJColor::fWutSecondaryColors[3]);
        cf->Draw("c same");
        cfFull->Draw("c same");
        tpt->Draw("same");

        c->Write();
        hist->Write();
        cf->Write();
        cfFull->Write();
    }
    radiiY.Write();
    chisqrY.Write();

    for (const auto &y : rapBinNums)
        for (const auto &kt : ktBinNums)
        {
            std::cout << "kt: " << kt << "\ty: " << y << "\n";
            hist = std::unique_ptr<TH1>(itp->Get<TH1D>(TString::Format("hQinvRatKt%dY%d",kt,y)));
            fitter.SetDataHistogram(std::move(hist));            

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

            func2.SetParameters({params[0],params[1]},{params[2]},{});
            auto cfFull = func2.Evaluate();
            cfFull->SetName(TString::Format("fitFullKt%dY%d",kt,y));
            cfFull->SetLineColor(JJColor::fWutSecondaryColors[3]);
            cfFull->SetLineStyle(kDashed);

            if (fitter.GetStatus() == 0 || fitter.GetStatus() == 1)
            {
                radii[y - 1].SetBinContent(kt,params.at(2));
                radii[y - 1].SetBinError(kt,errors.at(2));
                chisqr.SetBinContent(kt,y,fitter.GetMinValue());
            }
            else
            {
                radii[y - 1].SetBinContent(kt,0.);
                radii[y - 1].SetBinError(kt,0.);
                chisqr.SetBinContent(kt,y,0.);
            }

            std::unique_ptr<TCanvas> c(new TCanvas(TString::Format("canvKt%dY%d",kt,y),"",800,800));
            c->SetMargin(0.15,0.02,0.15,0.02);
            JJColor::CreatePrimaryWutGradient();

            hist->Draw();
            cf->SetLineColor(JJColor::fWutSecondaryColors[3]);
            cf->Draw("c same");
            cfFull->Draw("c same");
            tpt->Draw("same");

            c->Write();
            hist->Write();
            cf->Write();
            cfFull->Write();
        }

    for (const auto &hist : radii)
        hist.Write();
    
    chisqr.Write();
    otp->Save();
    otp->Close();
}