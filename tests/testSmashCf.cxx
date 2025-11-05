#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TChain.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "ROOT/TSeq.hxx"
#include "Math/LorentzVector.h"

#include "/home/jedkol/Downloads/HADES/JJFemtoMixer/JJFemtoMixer.hxx"
#include "/home/jedkol/Downloads/indicators/single_include/indicators/indicators.hpp"

#include "InteractionTermSchrodinger.hxx"

struct Track
{
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > m_idealMom, m_smearedMom, m_idealPos;
    int m_origin;
    double m_formationTime;
    Track(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &idealMom, 
    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &smearedMom, 
    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &idealPos, 
    int origin, double formationTime) : 
        m_idealMom(idealMom), m_smearedMom(smearedMom), m_idealPos(idealPos), m_origin(origin), m_formationTime(formationTime) {}
};

struct Event
{
    std::string m_ID;
    std::size_t m_nTracks;
    std::vector<std::shared_ptr<Track> > m_tracks;
    Event(std::string ID, std::size_t nTracks) : m_ID(ID), m_nTracks(nTracks) {}
    [[nodiscard]] const std::string& GetID() const noexcept {return m_ID;}
};

struct Pair
{
    std::shared_ptr<Track> m_track1, m_track2;
    double m_kStar, m_rStar, m_kT, m_rapidity, m_cosTheta;

    Pair(const std::shared_ptr<Track> &trck1, const std::shared_ptr<Track> &trck2) : m_track1(trck1), m_track2(trck2) 
    {
        CalcKinematics(m_track1, m_track2);
    }
    void CalcKinematics(const std::shared_ptr<Track> &track1, const std::shared_ptr<Track> &track2)
    {
        constexpr double FmToGev = 0.197327;
        constexpr double GevToFm = 1. / FmToGev;
        constexpr double GeVToMeV = 1000.;
        
        // Calculate pair variables
        double tPx = track1->m_idealMom.Px() + track2->m_idealMom.Px();
        double tPy = track1->m_idealMom.Py() + track2->m_idealMom.Py();
        double tPz = track1->m_idealMom.Pz() + track2->m_idealMom.Pz();
        double tE  = track1->m_idealMom.E() + track2->m_idealMom.E();
        double tPt = tPx * tPx + tPy * tPy;
        double tMt = tE * tE - tPz * tPz;
        double tM  = std::sqrt(tMt - tPt);
        tMt = std::sqrt(tMt);
        m_kT = std::sqrt(tPt);
        m_rapidity = (track1->m_idealMom.Rapidity() + track2->m_idealMom.Rapidity()) / 2.;

        // Boost to LCMS
        double tBeta = tPz / tE;
        double tGamma = tE / tMt;	
        double mKStarLong = tGamma * (track1->m_idealMom.Pz() - tBeta * track1->m_idealMom.E());
        double tE1L = tGamma * (track1->m_idealMom.E()  - tBeta * track1->m_idealMom.Pz());   
        
        // Rotate in transverse plane
        double mKStarOut  = ( track1->m_idealMom.Px() * tPx + track1->m_idealMom.Py() * tPy) / m_kT;
        double mKStarSide = (-track1->m_idealMom.Px() * tPy + track1->m_idealMom.Py() * tPx) / m_kT;

        // Boost to pair cms
        mKStarOut = tMt / tM * (mKStarOut - m_kT / tMt * tE1L);

        double tDX = track1->m_idealPos.X() * FmToGev - track2->m_idealPos.X() * FmToGev;
        double tDY = track1->m_idealPos.Y() * FmToGev - track2->m_idealPos.Y() * FmToGev;
        double mRLong = track1->m_idealPos.Z() * FmToGev - track2->m_idealPos.Z() * FmToGev;
        double mDTime = track1->m_idealPos.T() * FmToGev - track2->m_idealPos.T() * FmToGev;

        double mROut = (tDX * tPx + tDY * tPy) / m_kT;
        double mRSide = (-tDX * tPy + tDY * tPx) / m_kT;

        double mRSidePairCMS = mRSide;

        double mRLongPairCMS = tGamma * (mRLong - tBeta * mDTime);
        double mDTimePairLCMS = tGamma * (mDTime - tBeta * mRLong);

        tBeta = m_kT / tMt;
        tGamma = tMt / tM;

        double mROutPairCMS = tGamma * (mROut - tBeta * mDTimePairLCMS);

        m_kStar = std::sqrt(mKStarSide * mKStarSide + mKStarOut * mKStarOut + mKStarLong * mKStarLong);
        m_rStar = std::sqrt(mROutPairCMS * mROutPairCMS + mRSidePairCMS * mRSidePairCMS + mRLongPairCMS * mRLongPairCMS);
        m_cosTheta = (mKStarOut * mROutPairCMS + mKStarSide * mRSidePairCMS + mKStarLong * mRLongPairCMS) / (m_kStar * m_rStar);

        m_kStar *= GeVToMeV;
        m_rStar *= GevToFm;
    }
};

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
    const std::string inputfileBase = "/home/jedkol/lustre/hades/user/kjedrzej/SmashResults/AuAu_1p23AGeV_0_80cent_hardEoS/particle_list_";
    constexpr std::size_t numberOfFiles = 1;

    JJCorrFitter::InteractionTermSchrodinger wf;

    constexpr double momBegin = 0.5;
    constexpr double momStep = 1;
    std::size_t momCounter = -1;
    std::vector<double> momBins(500,0.);
    std::generate(momBins.begin(),momBins.end(),[&momBegin,&momStep,&momCounter]{return momBegin + momStep * (++momCounter);});
    wf.SetMomentumBins(std::move(momBins));

    constexpr double rBegin = 0.5;
    constexpr double rStep = 0.5;
    std::size_t rCounter = -1;
    std::vector<double> rBins(200);
    std::generate(rBins.begin(),rBins.end(),[&rBegin,&rStep,&rCounter]{return rBegin + rStep * (++rCounter);});
    wf.SetDistanceBins(std::move(rBins));

    constexpr std::size_t elems = 200;
    constexpr double stepCt = 2. / elems;
    std::size_t ctCounter = 0;
    std::vector<double> ctBins(elems,0.);
    std::generate(ctBins.begin(),ctBins.end(),[&stepCt,&ctCounter]{return -1 + (++ctCounter) * stepCt;});
    ctBins.push_back(1 + stepCt);
    wf.SetCosThetaBins(std::move(ctBins));

    wf.PopulateGrid();
    
    TChain *chain = new TChain("tree");
    for (const auto &i : ROOT::TSeqUL(numberOfFiles))
    {
        chain->Add(TString::Format("%s%ld.root",inputfileBase.c_str(),i));
    }

    TTreeReader reader(chain);
    TTreeReaderValue<int> nTracks(reader,"tracks");
    TTreeReaderValue<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > fourPosIdeal(reader,"posIdeal");
    TTreeReaderValue<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > fourMomIdeal(reader,"momIdeal");
    TTreeReaderValue<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > fourMomSmeared(reader,"momSmear");
    TTreeReaderValue<std::vector<int> > procTypeOrigin(reader,"typeOrigin");
    TTreeReaderValue<std::vector<double> > formationTime(reader,"formTime");

    indicators::BlockProgressBar bar{indicators::option::BarWidth{80},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::ShowPercentage{true},
                    indicators::option::ShowRemainingTime{true},
                    indicators::option::PrefixText{"Calculating CFs "},
                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                    indicators::option::MaxProgress{reader.GetEntries()}};

    std::shared_ptr<Event> event;
    std::shared_ptr<Track> track;

    Mixing::JJFemtoMixer<Event,Track,Pair> mixer;
    mixer.SetMaxBufferSize(0);
    // mixer.SetEventHashingFunction();
    // mixer.SetPairHashingFunction();
    // mixer.SetPairCuttingFunction();

    std::map<std::string,TH1D> histSign, histBckg;
    TH2D hWeight("hWeight",";cos(#theta);|#Psi|^{2}",100,-1,1,100,0,2);

    bool firstEntry = true;
    std::size_t counter = 1;

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

        event = std::make_shared<Event>(std::to_string(reader.GetCurrentEntry()),*nTracks);

        for (const auto &i : ROOT::TSeqUL(*nTracks))
        {
            track = std::make_shared<Track>(fourMomIdeal->at(i),fourMomSmeared->at(i),fourPosIdeal->at(i),procTypeOrigin->at(i),formationTime->at(i));
            event->m_tracks.push_back(track);
        }

        auto signal = mixer.AddEvent(event,event->m_tracks);

        for (const auto &[key,pairs] : signal)
        {
            for (const auto &pair : pairs)
            {
                if (histSign.find(key) == histSign.end())
                {
                    histSign.emplace(key, TH1D(TString::Format("hkStarSign_%s",key.data()),"Signal of Protons 0-10%% centrality;k^{*} [MeV/c];CF(k^{*})",750,0,1500));
                    histBckg.emplace(key, TH1D(TString::Format("hkStarBckg_%s",key.data()),"Backgound of Protons 0-10%% centrality;k^{*} [MeV/c];CF(k^{*})",750,0,1500));
                }

                wf.SetMomentum(pair->m_kStar);
                double weight = (pair->m_kStar >= 499) ? 1. : wf.GetValue(pair->m_rStar, pair->m_cosTheta);
                histSign.at(key).Fill(pair->m_kStar, weight);
                histBckg.at(key).Fill(pair->m_kStar);
                hWeight.Fill(pair->m_cosTheta, weight);
            }
        }

        if (--counter == 0)
        {
            bar.set_progress(reader.GetCurrentEntry());
            counter = 100;
        }
    }

    bar.mark_as_completed();
    indicators::show_console_cursor(true);
    
    // mixer.PrintStatus();
    
    TFile otpFile("result_hardEoS.root","recreate");

    hWeight.Write();
    for (const auto &[key,hist] : histSign)
        hist.Write();

    for (const auto &[key,hist] : histBckg)
        hist.Write();
}