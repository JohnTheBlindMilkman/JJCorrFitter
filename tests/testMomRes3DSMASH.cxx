#include "/home/jedkol/lxpool/virgoSmash/macros/JJSmashUtils.hxx"
#include "/home/jedkol/Downloads/indicators/single_include/indicators/indicators.hpp"
#include "../../HADES/HADES-CrAP/macros/MacroUtils.hxx"

#include <fstream>
#include <random>

#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TH3D.h"
#include "ROOT/TSeq.hxx"
#include "Math/LorentzVector.h"

#include "InteractionTermSchrodinger.hxx"

using LorenzVec = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >;

std::tuple<double,double,double,double,double,double> CalcKinematics(const LorenzVec &mom1, const LorenzVec &pos1, const LorenzVec &mom2, const LorenzVec &pos2)
{
    constexpr double FmToGev = 0.197327;
    constexpr double GevToFm = 1. / FmToGev;
    constexpr double GeVToMeV = 1000.;
    
    // Calculate pair variables
    double tPx = mom1.Px() + mom2.Px();
    double tPy = mom1.Py() + mom2.Py();
    double tPz = mom1.Pz() + mom2.Pz();
    double tE  = mom1.E() + mom2.E();
    double tPt = tPx * tPx + tPy * tPy;
    double tMt = tE * tE - tPz * tPz;
    double tM  = std::sqrt(tMt - tPt);
    tMt = std::sqrt(tMt);
    tPt = std::sqrt(tPt);

    // Boost to LCMS
    double tBeta = tPz / tE;
    double tGamma = tE / tMt;	
    double mKStarLong = tGamma * (mom1.Pz() - tBeta * mom1.E());
    double tE1L = tGamma * (mom1.E()  - tBeta * mom1.Pz());   
    
    // Rotate in transverse plane
    double mKStarOut  = ( mom1.Px() * tPx + mom1.Py() * tPy) / tPt;
    double mKStarSide = (-mom1.Px() * tPy + mom1.Py() * tPx) / tPt;

    // Boost to pair cms
    mKStarOut = tMt / tM * (mKStarOut - tPt / tMt * tE1L);

    double tDX = pos1.X() * FmToGev - pos2.X() * FmToGev;
    double tDY = pos1.Y() * FmToGev - pos2.Y() * FmToGev;
    double mRLong = pos1.Z() * FmToGev - pos2.Z() * FmToGev;
    double mDTime = pos1.T() * FmToGev - pos2.T() * FmToGev;

    double mROut = (tDX * tPx + tDY * tPy) / tPt;
    double mRSide = (-tDX * tPy + tDY * tPx) / tPt;

    double mRSidePairCMS = mRSide;

    double mRLongPairCMS = tGamma * (mRLong - tBeta * mDTime);
    double mDTimePairLCMS = tGamma * (mDTime - tBeta * mRLong);

    tBeta = tPt / tMt;
    tGamma = tMt / tM;

    double mROutPairCMS = tGamma * (mROut - tBeta * mDTimePairLCMS);

    double mKStar = sqrt(mKStarSide * mKStarSide + mKStarOut * mKStarOut + mKStarLong * mKStarLong);
    double mRStar = std::sqrt(mROutPairCMS * mROutPairCMS + mRSidePairCMS * mRSidePairCMS + mRLongPairCMS * mRLongPairCMS);
    double mCosTheta = (mKStarOut * mROutPairCMS + mKStarSide * mRSidePairCMS + mKStarLong * mRLongPairCMS) / (mKStar * mRStar);

    return std::make_tuple(mRStar * GevToFm, mKStar * GeVToMeV, mCosTheta, mKStarOut * GeVToMeV, mKStarSide * GeVToMeV, mKStarLong * GeVToMeV);
}

bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) 
{
   if (value.GetSetupStatus() < 0) 
   {
      std::cerr << "Error " << value.GetSetupStatus() << "setting up reader for " << value.GetBranchName() << '\n';
      return false;
   }

   return true;
}

int main()
{
    const std::string inputfileBase = "/home/jedkol/lustre/hades/user/kjedrzej/SmashResults/AuAu_1p23AGeV_0_10cent_new/particle_list_";
    constexpr std::size_t numberOfFiles = 10;
    
    TChain *chain = new TChain("tree");
    for (const auto &i : ROOT::TSeqUL(numberOfFiles))
    {
        chain->Add(TString::Format("%s%ld.root",inputfileBase.c_str(),i));
    }

    TTreeReader reader(chain);
    TTreeReaderValue<int> nTracks(reader,"tracks");
    TTreeReaderValue<std::vector<LorenzVec> > fourPosIdeal(reader,"posIdeal");
    TTreeReaderValue<std::vector<LorenzVec> > fourMomIdeal(reader,"momIdeal");
    TTreeReaderValue<std::vector<LorenzVec> > fourMomSmeared(reader,"momSmear");
    TTreeReaderValue<std::vector<int> > procTypeOrigin(reader,"typeOrigin");
    TTreeReaderValue<std::vector<double> > formationTime(reader,"formTime");

    indicators::BlockProgressBar bar{indicators::option::BarWidth{80},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::ShowPercentage{true},
                    indicators::option::ShowRemainingTime{true},
                    indicators::option::PrefixText{"Calculating CFs "},
                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};

    double kStarIdeal, kStarSmeared, rStarIdeal, rStarSmeared, cosThetaIdeal, cosThetaSmeared, weight;
    std::string tmp, inputFileName;
    int counter = 1;

    JJCorrFitter::InteractionTermSchrodinger wavefunction;

    constexpr double momBegin = 0.5;
    constexpr double momStep = 1;
    std::size_t momCounter = -1;
    std::vector<double> momBins(500,0.);
    std::generate(momBins.begin(),momBins.end(),[&momBegin,&momStep,&momCounter]{return momBegin + momStep * (++momCounter);});
    wavefunction.SetMomentumBins(std::move(momBins));

    constexpr double rBegin = 0.5;
    constexpr double rStep = 0.5;
    std::size_t rCounter = -1;
    std::vector<double> rBins(200);
    std::generate(rBins.begin(),rBins.end(),[&rBegin,&rStep,&rCounter]{return rBegin + rStep * (++rCounter);});
    wavefunction.SetDistanceBins(std::move(rBins));

    constexpr std::size_t elems = 200;
    constexpr double stepCt = 2. / elems;
    std::size_t ctCounter = 0;
    std::vector<double> ctBins(elems,0.);
    std::generate(ctBins.begin(),ctBins.end(),[&stepCt,&ctCounter]{return -1 + (++ctCounter) * stepCt;});
    ctBins.push_back(1 + stepCt);
    wavefunction.SetCosThetaBins(std::move(ctBins));

    wavefunction.PopulateGrid();

    TH3D *hIdealNum = new TH3D("hIdealNum",";q^{ideal}_{out} [MeV/c];q^{ideal}_{side} [MeV/c];q^{ideal}_{long} [MeV/c];C(q^{ideal}_{out},q^{ideal}_{side},q^{ideal}_{long})",125,0,500,125,0,500,125,0,500);
    TH3D *hIdealDen = new TH3D("hIdealDen",";q^{ideal}_{out} [MeV/c];q^{ideal}_{side} [MeV/c];q^{ideal}_{long} [MeV/c];C(q^{ideal}_{out},q^{ideal}_{side},q^{ideal}_{long})",125,0,500,125,0,500,125,0,500);
    TH3D *hSmearedNum = new TH3D("hSmearedNum",";q^{smeared}_{out} [MeV/c];q^{smeared}_{side} [MeV/c];q^{smeared}_{long} [MeV/c];C(q^{smeared}_{out},q^{smeared}_{side},q^{smeared}_{long})",125,0,500,125,0,500,125,0,500);
    TH3D *hSmearedDen = new TH3D("hSmearedDen",";q^{smeared}_{out} [MeV/c];q^{smeared}_{side} [MeV/c];q^{smeared}_{long} [MeV/c];C(q^{smeared}_{out},q^{smeared}_{side},q^{smeared}_{long})",125,0,500,125,0,500,125,0,500);

    bool firstEntry = true;
    indicators::show_console_cursor(false);

    while (reader.Next())
    {
        if (firstEntry) 
        {
            // Check that branches exist and their types match our expectation.
            if (!CheckValue(nTracks)) break;
            if (!CheckValue(fourPosIdeal)) break;
            if (!CheckValue(fourMomIdeal)) break;
            if (!CheckValue(fourMomSmeared)) break;
            if (!CheckValue(procTypeOrigin)) break;
            if (!CheckValue(formationTime)) break;
            firstEntry = false;
        }

        for (std::size_t it = 0; it < *nTracks; it++)
        {
            for (std::size_t jt = it + 1; jt < *nTracks; jt++)
            {
                const auto &[rStarIdeal,kStarIdeal,cosThetaIdeal,kStarOutIdeal,kStarSideIdeal,kStarLongIdeal] = CalcKinematics(fourMomIdeal->at(it),fourPosIdeal->at(it),fourMomIdeal->at(jt),fourPosIdeal->at(jt));
                const auto &[rStarSmeared,kStarSmeared,cosThetaSmeared,kStarOutSmeared,kStarSideSmeared,kStarLongSmeared] = CalcKinematics(fourMomSmeared->at(it),fourPosIdeal->at(it),fourMomSmeared->at(jt),fourPosIdeal->at(jt));

                wavefunction.SetMomentum(kStarIdeal);
                weight = (kStarIdeal >= 499) ? 1. : wavefunction.GetValue(rStarIdeal,cosThetaIdeal);

                hIdealNum->Fill(2 * kStarOutIdeal, 2 * kStarSideIdeal, 2 * kStarLongIdeal, weight);
                hIdealDen->Fill(2 * kStarOutIdeal, 2 * kStarSideIdeal, 2 * kStarLongIdeal);

                hSmearedNum->Fill(2 * kStarOutSmeared, 2 * kStarSideSmeared, 2 * kStarLongSmeared, weight);
                hSmearedDen->Fill(2 * kStarOutSmeared, 2 * kStarSideSmeared, 2 * kStarLongSmeared);
            }
        }

        if (--counter == 0)
        {
            bar.set_progress(100 * static_cast<float>(reader.GetCurrentEntry()) / reader.GetEntries());
            counter = 100;
        }
    }

    bar.mark_as_completed();
    indicators::show_console_cursor(true);

    TH3D *hRatioIdeal = static_cast<TH3D*>(hIdealNum->Clone("hRatioIdeal"));
    hRatioIdeal->Divide(hIdealDen);
    hRatioIdeal->SetTitle(";q^{ideal}_{out} [MeV/c];q^{ideal}_{side} [MeV/c];q^{ideal}_{long} [MeV/c];C(q^{ideal}_{out},q^{ideal}_{side},q^{ideal}_{long})");
    JJUtils::Generic::SetErrorsDivide(hRatioIdeal,hIdealNum,hIdealDen);

    TH3D *hRatioSmeared = static_cast<TH3D*>(hSmearedNum->Clone("hRatioSmeared"));
    hRatioSmeared->Divide(hSmearedDen);
    hRatioSmeared->SetTitle(";q^{smeared}_{out} [MeV/c];q^{smeared}_{side} [MeV/c];q^{smeared}_{long} [MeV/c];C(q^{smeared}_{out},q^{smeared}_{side},q^{smeared}_{long})");
    JJUtils::Generic::SetErrorsDivide(hRatioSmeared,hSmearedNum,hSmearedDen);

    TH3D *hRatio = static_cast<TH3D*>(hRatioIdeal->Clone("hRatio"));
    hRatio->Divide(hRatioSmeared);
    //hRatio->Sumw2();
    hRatio->SetTitle(";q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];C(q^{ideal}_{out},q^{ideal}_{side},q^{ideal}_{long})/C(q^{smeared}_{out},q^{smeared}_{side},q^{smeared}_{long})");
    //JJUtils::Generic::SetErrors(hRatio,hRatioIdeal,hRatioSmeared); // IMO those are not correlated

    std::unique_ptr<TFile> otpFile(TFile::Open("momRes3DSmashInvP_Acceptance.root","recreate"));
    hIdealNum->Write();
    hIdealDen->Write();
    hSmearedNum->Write();
    hSmearedDen->Write();
    hRatioIdeal->Write();
    hRatioSmeared->Write();
    hRatio->Write();
}