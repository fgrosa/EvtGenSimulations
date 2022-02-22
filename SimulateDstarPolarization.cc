R__LOAD_LIBRARY(EvtGen)
R__ADD_INCLUDE_PATH($EVTGEN_ROOT/include)

#include <array>
#include <string>
#include <vector>
#include <map>
#include <deque>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TSpline.h>
#include <TNtuple.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Pythia8/Pythia.h"
#include "EvtGen/EvtGen.hh"
#include "Pythia8Plugins/EvtGen.h"

using namespace Pythia8;

namespace
{
    enum decayer
    {
        kPythia8 = 0,
        kEvtGen
    };

    enum tunes
    {
        kMonash = 0,
        kCRMode0,
        kCRMode2,
        kCRMode3
    };

    enum processes
    {
        kSoftQCD = 0,
        kHardQCD
    };
}

//__________________________________________________________________________________________________
void SimulateDstarPolarization(int nEvents=10000, int decayer=kEvtGen, int tune=kMonash, int process=kHardQCD, float energy=13000, int seed=42, std::string outFileNameRoot="AnalysisResults_Dstar_polarization_pythia_HardQCD_Monash_EvtGen.root");
template<typename T, typename T2>
bool IsFromBeauty(T& mothers, T2& pythia, int& motherIdx);

//__________________________________________________________________________________________________
void SimulateDstarPolarization(int nEvents, int decayer, int tune, int process, float energy, int seed, std::string outFileNameRoot)
{
    //__________________________________________________________
    // create and configure pythia generator

    Pythia pythia;
    if(process == kSoftQCD)
    {
        pythia.readString("SoftQCD:all = on");
    }
    else if(process == kHardQCD)
    {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    }

    // set tune
    if(tune == kMonash)
    {
        pythia.readString(Form("Tune:pp = 14"));
    }
    else if(tune == kCRMode0)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode2)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode3)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }

    EvtGenDecays *evtgen = nullptr;
    if(decayer == kPythia8) 
    {
        // keep only interesting decays
        pythia.readString("421:onMode = off");
        pythia.readString("413:onMode = off");
        pythia.readString("421:onIfMatch = 211 321");
        pythia.readString("413:onIfMatch = 211 421");
    }
    else if(decayer == kEvtGen)
    {
        // keep only interesting decays
        evtgen = new EvtGenDecays(&pythia, "./DECAY_2010.DEC", "./evt.pdl");
        evtgen->readDecayFile("DstarDecay.dec");
    }
    else {
        std::cerr << "Invalid decayer selected! Exit" << std::endl;
        return;
    }

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    //__________________________________________________________
    // define outputs
    TNtuple* tupleDstars = new TNtuple("tupleDstars", "tupleDstars", "pTDstar:yDstar:cosThetaStar:cosThetaStarB:pTD0:yD0:pTPi:yPi:fromB");
    auto hPxCorr = new TH2F("hPxCorr", ";#it{p}_{x}(D*^{+}) lab frame;#it{p}_{x}(B) D*^{+} rest frame", 200, -10., 10., 200, -10., 10.); 
    auto hPyCorr = new TH2F("hPyCorr", ";#it{p}_{x}(D*^{y}) lab frame;#it{p}_{y}(B) D*^{+} rest frame", 200, -10., 10., 200, -10., 10.); 
    auto hPzCorr = new TH2F("hPzCorr", ";#it{p}_{x}(D*^{z}) lab frame;#it{p}_{z}(B) D*^{+} rest frame", 200, -10., 10., 200, -10., 10.); 
    auto hPhiCorr = new TH2F("hPhiCorr", ";#varphi(D*^{+}) lab frame;#varphi(B) D*^{+} rest frame", 100, 0., TMath::Pi(), 100, 0., TMath::Pi()); 
    auto hEtaCorr = new TH2F("hEtaCorr", ";#eta(D*^{+}) lab frame;#eta(B) D*^{+} rest frame", 100, -2.5, 2.5, 100, -2.5, 2.5); 

    //__________________________________________________________
    // perform the simulation
    for (auto iEvent{0}; iEvent<nEvents; ++iEvent)
    {
        if(!pythia.next())
            continue;
        if(evtgen)
            evtgen->decay();
        for(auto iPart{2}; iPart<pythia.event.size(); ++iPart)
        {
            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);
            if (absPdg != 413)
                continue;

            float pT = pythia.event[iPart].pT();
            float px = pythia.event[iPart].px();
            float py = pythia.event[iPart].py();
            float pz = pythia.event[iPart].pz();
            float phi = pythia.event[iPart].phi();
            float eta = pythia.event[iPart].eta();
            float p = std::sqrt(px*px + py*py + pz*pz);
            float y = pythia.event[iPart].y();
            if (std::abs(y) > 1.)
                continue;

            auto mothers = pythia.event[iPart].motherList();
            auto dauList = pythia.event[iPart].daughterList();

            int motherIdx = -1;
            bool isFromB = IsFromBeauty(mothers, pythia, motherIdx);

            std::array<float, 2> ptDau{}, yDau{};
            ROOT::Math::PxPyPzMVector fourVecD0, fourVecPi;
            for (auto& dau : dauList)
            {
                auto absPdgDau = std::abs(pythia.event[dau].id());
                if (absPdgDau == 421) {
                    ptDau[0] = pythia.event[dau].pT();
                    yDau[0] = pythia.event[dau].y();
                    fourVecD0 = ROOT::Math::PxPyPzMVector(pythia.event[dau].px(), pythia.event[dau].py(), pythia.event[dau].pz(), TDatabasePDG::Instance()->GetParticle(421)->Mass());
                }
                else if (absPdgDau == 211) {
                    ptDau[1] = pythia.event[dau].pT();
                    yDau[1] = pythia.event[dau].y();
                    fourVecPi = ROOT::Math::PxPyPzMVector(pythia.event[dau].px(), pythia.event[dau].py(), pythia.event[dau].pz(), TDatabasePDG::Instance()->GetParticle(211)->Mass());
                }
            }

            float sumPt = 0., sumY = 0.;
            if(std::accumulate(ptDau.begin(), ptDau.end(), sumPt) < 1.e-10 || std::abs(std::accumulate(yDau.begin(), yDau.end(), sumY)) < 1.e-10)
                continue;

            ROOT::Math::PxPyPzMVector fourVecDstar = fourVecPi + fourVecD0;
            ROOT::Math::Boost boostv12{fourVecDstar.BoostToCM()};
            ROOT::Math::PxPyPzMVector fourVecD0CM = boostv12(fourVecD0);
            ROOT::Math::XYZVector threeVecD0CM = fourVecD0CM.Vect();
            ROOT::Math::XYZVector helicityVec = ROOT::Math::XYZVector(px / p, py / p, pz / p);
            float cosThetaStar = std::abs(helicityVec.Dot(threeVecD0CM) / std::sqrt(threeVecD0CM.Mag2()));
            float pxB = 0., pyB = 0., pzB = 0.;

            if(isFromB && motherIdx >= 0) {
                auto pdgMotherB = std::abs(pythia.event[motherIdx].id());
                if(isFromB && ( (pdgMotherB > 500 && pdgMotherB < 600) || (pdgMotherB > 5000 && pdgMotherB < 6000) ) ) {
                    auto fourVecB = ROOT::Math::PxPyPzMVector(pythia.event[motherIdx].px(), pythia.event[motherIdx].py(), pythia.event[motherIdx].pz(), TDatabasePDG::Instance()->GetParticle(pdgMotherB)->Mass());
                    ROOT::Math::PxPyPzMVector fourVecBCM = boostv12(fourVecB);
                    pxB = fourVecBCM.Px();
                    pyB = fourVecBCM.Py();
                    pzB = fourVecBCM.Pz();
                    hPxCorr->Fill(px, pxB);
                    hPyCorr->Fill(py, pyB);
                    hPzCorr->Fill(pz, pzB);
                    hPhiCorr->Fill(phi, fourVecBCM.Phi());
                    hEtaCorr->Fill(eta, fourVecBCM.Eta());
                }
            }

            float pB = std::sqrt(pxB*pxB + pyB*pyB + pzB*pzB);
            if (pB == 0.)
                pB = 1;
            ROOT::Math::XYZVector helicityVecB = ROOT::Math::XYZVector(pxB / pB, pyB / pB, pzB / pB);
            float cosThetaStarB = std::abs(helicityVecB.Dot(threeVecD0CM) / std::sqrt(threeVecD0CM.Mag2()));

            float array4tuple[9] = {pT, y, cosThetaStar, cosThetaStarB, ptDau[0], yDau[0], ptDau[1], yDau[1], (float)isFromB};
            tupleDstars->Fill(array4tuple);
        }
    }

    // save root output file
    TFile outFile(outFileNameRoot.data(), "recreate");
    tupleDstars->Write();
    hPxCorr->Write();
    hPyCorr->Write();
    hPzCorr->Write();
    hPhiCorr->Write();
    hEtaCorr->Write();
    outFile.Close();
}

//__________________________________________________________________________________________________
template<typename T, typename T2>
bool IsFromBeauty(T& mothers, T2& pythia, int& motherIdx)
{
    for(auto& mom : mothers)
    {
        int absPdgMom = std::abs(pythia.event[mom].id());
        if(absPdgMom == 5 || absPdgMom/100 == 5 || absPdgMom/1000 == 5 ||
           (absPdgMom-10000)/100 == 5 || (absPdgMom-20000)/100 == 5 || (absPdgMom-30000)/100 == 5 ||
           (absPdgMom-100000)/100 == 5 || (absPdgMom-200000)/100 == 5 || (absPdgMom-300000)/100 == 5)
        {  // remove beauty feed-down
            motherIdx = mom;
            return true;
        }
    }
    return false;
}
