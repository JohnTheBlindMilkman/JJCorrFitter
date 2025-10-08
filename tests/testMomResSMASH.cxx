#include "/home/jedkol/lxpool/virgoSmash/macros/JJSmashUtils.hxx"
#include "/home/jedkol/Downloads/indicators/single_include/indicators/indicators.hpp"
#include "../../HADES/HADES-CrAP/macros/MacroUtils.hxx"

#include <fstream>
#include <random>

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"

#include "InteractionTermSchrodinger.hxx"

// I'm so sorry for this

TFile *inpFile = TFile::Open("/home/jedkol/lxpool/hades-crap/output/momentum_resolution.root");
TF1 *fMomMu = inpFile->Get<TF1>("fMomMuFit");
TF1 *fMomSig = inpFile->Get<TF1>("fMomSigFit");
TF1 *fPhiMu = inpFile->Get<TF1>("fPhiMuFit");
TF1 *fPhiSig = inpFile->Get<TF1>("fPhiSigFit");
TF1 *fThetaMu = inpFile->Get<TF1>("fThetaMuFit");
TF1 *fThetaSig = inpFile->Get<TF1>("fThetaSigFit");

std::random_device rd;
std::mt19937 mersenneTwister(rd());

std::tuple<double,double,double> CalcKinematics(const JJUtils::ParticleOscar97 &part1, const JJUtils::ParticleOscar97 &part2)
{
    constexpr double FmToGev = 0.197327;
    constexpr double GevToFm = 1. / FmToGev;
    constexpr double GeVToMeV = 1000.;
    
    // Calculate pair variables
    double tPx = part1.px + part2.px;
    double tPy = part1.py + part2.py;
    double tPz = part1.pz + part2.pz;
    double tE  = part1.p0 + part2.p0;
    double tPt = tPx * tPx + tPy * tPy;
    double tMt = tE * tE - tPz * tPz;
    double tM  = std::sqrt(tMt - tPt);
    tMt = std::sqrt(tMt);
    tPt = std::sqrt(tPt);

    // Boost to LCMS
    double tBeta = tPz / tE;
    double tGamma = tE / tMt;	
    double mKStarLong = tGamma * (part1.pz - tBeta * part1.p0);
    double tE1L = tGamma * (part1.p0  - tBeta * part1.pz);   
    
    // Rotate in transverse plane
    double mKStarOut  = ( part1.px*tPx + part1.py*tPy)/tPt;
    double mKStarSide = (-part1.px*tPy + part1.py*tPx)/tPt;

    // Boost to pair cms
    mKStarOut = tMt/tM * (mKStarOut - tPt/tMt * tE1L);

    double tDX = part1.x * FmToGev - part2.x * FmToGev;
    double tDY = part1.y * FmToGev - part2.y * FmToGev;
    double mRLong = part1.z * FmToGev - part2.z * FmToGev;
    double mDTime = part1.t * FmToGev - part2.t * FmToGev;

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

    return std::make_tuple(mRStar * GevToFm, mKStar * GeVToMeV, mCosTheta);
}

void smear(JJUtils::ParticleOscar13 &part)
{
    constexpr float GeVtoMeV{1000.f};
    constexpr float MeVtoGeV{0.001f};
    constexpr float DegToRad{3.14159f / 180.f};

    double mom = part.fourMom.P();
    double phi = part.fourMom.Phi();
    double theta = part.fourMom.Theta();

    std::normal_distribution<double> rngMom(fMomMu->Eval(mom * GeVtoMeV),fMomSig->Eval(mom * GeVtoMeV));
    std::normal_distribution<double> rngPhi(fPhiMu->Eval(mom * GeVtoMeV),fPhiSig->Eval(mom * GeVtoMeV));
    std::normal_distribution<double> rngTheta(fThetaMu->Eval(mom * GeVtoMeV),fThetaSig->Eval(mom * GeVtoMeV));

    mom = 1. / mom;
    mom += rngMom(mersenneTwister) * MeVtoGeV;
    mom = 1. / mom;
    phi += rngPhi(mersenneTwister) * DegToRad;
    theta += rngTheta(mersenneTwister) * DegToRad;

    part.px = mom * std::cos(phi) * std::sin(theta);
    part.py = mom * std::sin(phi) * std::sin(theta);
    part.pz = mom * std::cos(theta);
    part.p0 = std::sqrt(part.px * part.px + part.py * part.py + part.pz * part.pz + part.mass * part.mass);
}

int main()
{
    const std::string inputfileBase{"/home/jedkol/lustre/hades/user/kjedrzej/SmashResults/AuAu_1p23AGeV_0_10cent_new/particle_list_"};
    constexpr std::size_t numberOfEvents = 10000;
    constexpr std::size_t numberOfFiles = 10;
    
    indicators::BlockProgressBar bar{indicators::option::BarWidth{80},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::ShowPercentage{true},
                    indicators::option::ShowRemainingTime{true},
                    indicators::option::PrefixText{"Smearing momenta "},
                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                    indicators::option::MaxProgress{numberOfEvents*numberOfFiles}};

    double kStarIdeal, kStarSmeared, rStarIdeal, rStarSmeared, cosThetaIdeal, cosThetaSmeared, weight;
    std::string tmp, inputFileName;
    int numberOfTracks, counter = 1;
    std::vector<JJUtils::ParticleOscar13> idealArray,smearedArray;

    JJCorrFitter::InteractionTermSchrodinger wavefunction;
    std::ifstream inputfile;

    TH1D *hIdealNum = new TH1D("hIdealNum",";q^{ideal}_{inv} [MeV/c];C(q^{ideal}_{inv})",750,0,3000);
    TH1D *hIdealDen = new TH1D("hIdealDen",";q^{ideal}_{inv} [MeV/c];C(q^{ideal}_{inv})",750,0,3000);
    TH1D *hSmearedNum = new TH1D("hSmearedNum",";q^{smeared}_{inv} [MeV/c];C(q^{smeared}_{inv})",750,0,3000);
    TH1D *hSmearedDen = new TH1D("hSmearedDen",";q^{smeared}_{inv} [MeV/c];C(q^{smeared}_{inv})",750,0,3000);
    TH1D *hRStarDiff = new TH1D("hRStarDiff",";#left| k^{*}_{ideal} - k^{*}_{smeared} #right| [fm];",1000,0,1);

    indicators::show_console_cursor(false);

    for (std::size_t file = 0; file < numberOfFiles; ++file)
    {
        inputFileName = inputfileBase + std::to_string(file) + ".oscar";
        inputfile.open(inputFileName);
        if (! inputfile.is_open())
        {
            std::string msg = "inputFile " + inputFileName + " could not be opened";
            continue;
        }

        for (int i : {1,2,3}) // omit first 3 lines of the inputfile
        {
            std::getline(inputfile,tmp);
        }

        for (std::size_t event = 0; event < numberOfEvents; ++event)
        {
            inputfile >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> numberOfTracks; // get the number of tracks for each event

            idealArray.clear();
            idealArray.resize(0);
            smearedArray.clear();
            smearedArray.resize(0);

            for (int track = 0; track < numberOfTracks; ++track)
            {
                JJUtils::ParticleOscar13 part;
                inputfile >> part;

                if (part.pdg == 2212)
                {
                    JJUtils::PropagatePositions(part);
                    idealArray.push_back(part);
                    smear(part);
                    smearedArray.push_back(part);
                }
            }

            for (std::size_t it = 0; it < idealArray.size(); it++)
            {
                for (std::size_t jt = it + 1; jt < idealArray.size(); jt++)
                {
                    std::tie(rStarIdeal,kStarIdeal,cosThetaIdeal) = CalcKinematics(JJUtils::ParticleOscar97(idealArray[it]),JJUtils::ParticleOscar97(idealArray[jt]));
                    std::tie(rStarSmeared,kStarSmeared,cosThetaSmeared) = CalcKinematics(JJUtils::ParticleOscar97(smearedArray[it]),JJUtils::ParticleOscar97(smearedArray[jt]));
                    
                    hRStarDiff->Fill(std::abs(rStarIdeal - rStarSmeared));

                    wavefunction.SetMomentum(kStarIdeal);
                    weight = wavefunction.GetValue(rStarIdeal,cosThetaIdeal);

                    hIdealNum->Fill(2 * kStarIdeal, weight);
                    hIdealDen->Fill(2 * kStarIdeal);

                    hSmearedNum->Fill(2 * kStarSmeared, weight);
                    hSmearedDen->Fill(2 * kStarSmeared);
                }
            }

            if (--counter == 0)
            {
                bar.set_progress(file * numberOfEvents + event);
                counter = 100;
            }

            std::getline(inputfile,tmp); //omit last line
        }

        inputfile.close();
    }

    bar.mark_as_completed();
    indicators::show_console_cursor(true);

    TH1D *hRatioIdeal = static_cast<TH1D*>(hIdealNum->Clone("hRatioIdeal"));
    hRatioIdeal->Divide(hIdealDen);
    hRatioIdeal->SetTitle(";q_{inv} [MeV/c];C(q^{ideal}_{inv})");
    JJUtils::Generic::SetErrorsDivide(hRatioIdeal,hIdealNum,hIdealDen);

    TH1D *hRatioSmeared = static_cast<TH1D*>(hSmearedNum->Clone("hRatioSmeared"));
    hRatioSmeared->Divide(hSmearedDen);
    hRatioSmeared->SetTitle(";q_{inv} [MeV/c];C(q^{smeared}_{inv})");
    JJUtils::Generic::SetErrorsDivide(hRatioSmeared,hSmearedNum,hSmearedDen);

    TH1D *hRatio = static_cast<TH1D*>(hRatioIdeal->Clone("hRatio"));
    hRatio->Divide(hRatioSmeared);
    //hRatio->Sumw2();
    hRatio->SetTitle(";q_{inv} [MeV/c];C(q^{ideal}_{inv})/C(q^{smeared}_{inv})");
    //JJUtils::Generic::SetErrors(hRatio,hRatioIdeal,hRatioSmeared); // IMO those are not correlated

    std::unique_ptr<TFile> otpFile(TFile::Open("momRes1DSmashInvP_Acceptance.root","recreate"));
    hIdealNum->Write();
    hIdealDen->Write();
    hSmearedNum->Write();
    hSmearedDen->Write();
    hRatioIdeal->Write();
    hRatioSmeared->Write();
    hRatio->Write();
    hRStarDiff->Write();
}