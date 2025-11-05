#include "TFile.h"
#include "TCanvas.h"
#include "CorrelationFunction3D.hxx"
#include "SourceFunction3D.hxx"
#include "InteractionTermSchrodinger.hxx"

#include <memory>
#include <chrono>
#include <random>

void SimulateErrors(const std::unique_ptr<TH1> &data)
{
    static constexpr auto linFunc = [](double x){return -0.01 * x + 1;};
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution<double> dist(100,50);

    const int binsX = data->GetNbinsX();
    const int binsY = data->GetNbinsY();
    const int binsZ = data->GetNbinsZ();

    for (int binX = 1; binX <= binsX; ++binX)
        for (int binY = 1; binY <= binsY; ++binY)
            for (int binZ = 1; binZ <= binsZ; ++binZ)
            {
                const double qOut = data->GetXaxis()->GetBinCenter(binX);
                const double qSide = data->GetYaxis()->GetBinCenter(binY);
                const double qLong = data->GetZaxis()->GetBinCenter(binZ);
                const double qInv = std::sqrt(qOut * qOut + qSide * qSide + qLong * qLong);

                const double error = data->GetBinError(binX,binY,binZ);
                data->SetBinError(binX,binY,binZ,(qInv <= 100) ? error * (1 + linFunc(qInv) * dist(gen)) : error);
            }
}

void SimulatePoints(const std::unique_ptr<TH1> &data)
{
    static constexpr auto linFunc = [](double x){return -0.01 * x + 1;};
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution<double> dist(0,0.01);

    const int binsX = data->GetNbinsX();
    const int binsY = data->GetNbinsY();
    const int binsZ = data->GetNbinsZ();

    for (int binX = 1; binX <= binsX; ++binX)
        for (int binY = 1; binY <= binsY; ++binY)
            for (int binZ = 1; binZ <= binsZ; ++binZ)
            {
                const double qOut = data->GetXaxis()->GetBinCenter(binX);
                const double qSide = data->GetYaxis()->GetBinCenter(binY);
                const double qLong = data->GetZaxis()->GetBinCenter(binZ);
                const double qInv = std::sqrt(qOut * qOut + qSide * qSide + qLong * qLong);

                const double cont = data->GetBinContent(binX,binY,binZ);
                data->SetBinContent(binX,binY,binZ,(qInv <= 100) ? cont * (1 + linFunc(qInv) * dist(gen)) : cont);
            }
}

int main()
{
    // load data
    // std::unique_ptr<TFile> itp(TFile::Open("/home/jedkol/lxpool/hades-crap/output/3Dcorr_0_10_cent_forHAL.root"));
    // std::unique_ptr<TH1> hist(itp->Get<TH3D>("hQoslRatInteg"));

    std::unique_ptr<TFile> otp(TFile::Open("draw3DCF.root","recreate"));
    // hist->Write();

    JJCorrFitter::CorrelationFunction3D func1(
        std::make_unique<JJCorrFitter::SourceFunction3D>(),
        std::make_unique<JJCorrFitter::InteractionTermSchrodinger>()
        );

    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>("c","",800,800);
    c->SetMargin(0.15,0.02,0.15,0.02);

    // func1.SetBinning(hist,0,200);
    func1.SetBinning("fit_gauss","",50,0,200);
    func1.SetParameters({1.,1.},{2,3,4},{});
    auto start = std::chrono::high_resolution_clock::now();
    auto data = func1.EvaluateAtPlanes(1);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Evaluation time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "\n";

    SimulateErrors(data);
    SimulatePoints(data);
    data->Write("hQoslRatInteg_gauss");
}