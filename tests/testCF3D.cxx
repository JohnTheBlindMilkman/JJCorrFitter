#include "TFile.h"
#include "TCanvas.h"
#include "CorrelationFunction3D.hxx"
#include "SourceFunction3D.hxx"
#include "InteractionTermSchrodinger.hxx"

#include <memory>
#include <chrono>

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

    std::unique_ptr<TCanvas> c(new TCanvas("c","",800,800));
    c->SetMargin(0.15,0.02,0.15,0.02);

    // func1.SetBinning(hist,0,200);
    func1.SetBinning("fit_gauss","",50,0,200);
    func1.SetParameters({1.,1.},{2,3,4},{});
    auto start = std::chrono::high_resolution_clock::now();
    func1.EvaluateAtEdges()->Write("hQoslRatInteg_gauss");
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Evaluation time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "\n";
}