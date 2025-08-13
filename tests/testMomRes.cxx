#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "CorrelationFunction1D.hxx"
#include "SourceFunction1D.hxx"
#include "InteractionTermSchrodinger.hxx"
#include "../../HADES/HADES-CrAP/macros/MacroUtils.hxx"

std::unique_ptr<TH1D> FromKStarToQInv(const std::unique_ptr<TH1> &hist)
{
    const int nBins = hist->GetNbinsX();
    const double newMin = hist->GetXaxis()->GetXmin() ;
    const double newMax = hist->GetXaxis()->GetXmax() * 2.; // "*2" to convert from k* to qinv

    std::unique_ptr<TH1D> hOtp(new TH1D(hist->GetName(),hist->GetTitle(),nBins,newMin,newMax));
    hOtp->GetXaxis()->SetTitle("q_{inv} [MeV/c]");

    for (int i = 1; i <= nBins; ++i)
    {
        hOtp->SetBinContent(i,hist->GetBinContent(i));
        hOtp->SetBinError(i,hist->GetBinError(i));
    }

    return std::move(hOtp);
}

double Integral(const std::unique_ptr<TH1D> &cfIdeal, const std::unique_ptr<TH1D> &sigmaVals, double measured, double low, double high)
{
    double sumNum = 0, sumDen = 0;
    int nBins = 0;

    int lowBin = cfIdeal->FindBin(low);
    int highBin = cfIdeal->FindBin(high);

    for (int bin = lowBin; bin <= highBin; ++bin)
    {
        sumNum += cfIdeal->GetBinContent(bin) * cfIdeal->GetBinCenter(bin) * cfIdeal->GetBinCenter(bin) * std::exp(-std::pow(measured - cfIdeal->GetBinCenter(bin),2) / (2 * sigmaVals->GetBinContent(bin) * sigmaVals->GetBinContent(bin)));
        sumDen += cfIdeal->GetBinCenter(bin) * cfIdeal->GetBinCenter(bin) * std::exp(-std::pow(measured - cfIdeal->GetBinCenter(bin),2) / (2 * sigmaVals->GetBinContent(bin) * sigmaVals->GetBinContent(bin)));
        ++nBins;
    }
    return (sumDen > 0) ? sumNum / sumDen : 0.; // no need to do "(sumNum / nBins) / (sumDen / nBins);" nBins cancels out
}

double IntegralErr(const std::unique_ptr<TH1D> &cfIdeal, const std::unique_ptr<TH1D> &sigmaVals, double measured, double low, double high)
{
    double sumNum = 0, sumDen = 0;
    int nBins = 0;

    int lowBin = cfIdeal->FindBin(low);
    int highBin = cfIdeal->FindBin(high);

    // TODO: propagate uncertainty; right now I just took the equation and replaced values with errors
    for (int bin = lowBin; bin <= highBin; ++bin)
    {
        sumNum += cfIdeal->GetBinError(bin) * cfIdeal->GetBinCenter(bin) * cfIdeal->GetBinCenter(bin) * std::exp(-std::pow(measured - cfIdeal->GetBinCenter(bin),2) / (2 * sigmaVals->GetBinError(bin) * sigmaVals->GetBinError(bin)));
        sumDen += cfIdeal->GetBinCenter(bin) * cfIdeal->GetBinCenter(bin) * std::exp(-std::pow(measured - cfIdeal->GetBinCenter(bin),2) / (2 * sigmaVals->GetBinError(bin) * sigmaVals->GetBinError(bin)));
        ++nBins;
    }
    return (sumDen > 0) ? sumNum / sumDen : 0.; // no need to do "(sumNum / nBins) / (sumDen / nBins);" nBins cancels out
}

std::unique_ptr<TH1D> CalculateSmearedFunction(const std::unique_ptr<TH1D> &cfIdeal, const std::unique_ptr<TH1D> &sigmaVals)
{
    double kStarIdeal = 0;
    int nBins = cfIdeal->GetNbinsX();
    std::unique_ptr<TH1D> hOut(static_cast<TH1D*>(cfIdeal->Clone()));

    for (int bin = 1; bin <= nBins; ++bin)
    {
        kStarIdeal = cfIdeal->GetBinCenter(bin);
        hOut->SetBinContent(
            bin,
            Integral(
                cfIdeal,
                sigmaVals,
                kStarIdeal,
                kStarIdeal - sigmaVals->GetBinContent(bin),
                kStarIdeal + sigmaVals->GetBinContent(bin))
            );
        hOut->SetBinError(
            bin,
            IntegralErr(
                cfIdeal,
                sigmaVals,
                kStarIdeal,
                kStarIdeal - sigmaVals->GetBinContent(bin),
                kStarIdeal + sigmaVals->GetBinContent(bin))
            );
        //hOut->SetBinError(bin,0);
    }

    return std::move(hOut);
}

int main()
{
    std::unique_ptr<TFile> inp(TFile::Open("/home/jedkol/lxpool/hades-crap/output/1DMomResSigma_0_40_cent.root"));
    std::unique_ptr<TH1D> hist(inp->Get<TH1D>("hQinvSigma"));

    JJCorrFitter::CorrelationFunction1D func1(
        std::make_unique<JJCorrFitter::SourceFunction1D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    func1.SetBinning("hCFIdeal","",hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax() / 2.);

    std::unique_ptr<TFile> otp(TFile::Open("momRes1D.root","recreate"));

    for (const auto& radius : {1.})
    {
        func1.SetParameters({1,0.902764},{2.82633},{});
        std::unique_ptr<TH1> hCFIdealKstar = func1.Evaluate();

        std::unique_ptr<TH1D> hCFIdealQinv = FromKStarToQInv(hCFIdealKstar);
        hCFIdealQinv->SetName(TString::Format("hCFIdeal_%d",static_cast<int>(radius)));

        std::unique_ptr<TH1D> hCFSmearedQinv = CalculateSmearedFunction(hCFIdealQinv,hist);
        hCFSmearedQinv->SetName(TString::Format("hCFSmeared_%d",static_cast<int>(radius)));

        std::unique_ptr<TH1D> hRatio(static_cast<TH1D*>(hCFIdealQinv->Clone(TString::Format("hRatio_%d",static_cast<int>(radius)))));
        hRatio->Divide(hCFSmearedQinv.get());
        hRatio->Sumw2();
        //JJUtils::Generic::SetErrors(hRatio.get(),hCFIdealQinv.get(),hCFSmearedQinv.get());

        otp->cd();
        hCFSmearedQinv->Write();
        hCFIdealQinv->Write();
        hRatio->Write();
    }
}