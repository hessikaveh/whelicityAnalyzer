// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "AnalysisDataFormats/TopObjects/interface/TtEvent.h"
#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "ttamwtsolver.h"
#include "MSETools.h"

using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
typedef math::XYZPoint Point;
typedef math::XYZTLorentzVector myLorentzVector;
struct MyLepton
{
    string type;
    int flavour;
    float mass;
    //pat::Particle pat_particle;
    // pat::Electron pat_particle;
    pat::GenericParticle pat_particle;
    //pat::Lepton<reco::Particle> pat_particle;

};
struct by_pt{
    bool operator()(MyLepton const &l1, MyLepton const &l2){
        return l1.pat_particle.pt() > l2.pat_particle.pt();
    }
};

//
// class declaration
//



class MiniAnalyzer : public edm::EDAnalyzer
{
public:
    explicit MiniAnalyzer (const edm::ParameterSet&);
    ~MiniAnalyzer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const& iRun, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    bool isMediumMuon(const reco::Muon & recoMu);
    double SF( double x);
 void printCutFlowResult(vid::CutFlowResult &cutflow);

    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronToken_;
    edm::EDGetTokenT<pat::TauCollection> tauToken_;
    edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
    edm::EDGetTokenT<pat::JetCollection>  jetToken_;
    edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;
    edm::EDGetTokenT<LHEEventProduct> lheEvtInfo_;
    edm::EDGetTokenT<double> rho_;

    //Boolean variables
    bool isPosMu = false;
    bool isNegMu = false;
    bool isPosEl = false;
    bool isNegEl = false;
    bool isDiMuon = false;
    bool isDiElectron = false;
    bool isElMu = false;
    bool isMuEl = false;
    bool isDiLeptonic = false;
const float m_el = 0.000511;
       const float m_mu = 0.10566;

    // function declaration
    bool IsSoftMuon(const pat::Muon&, const reco::Vertex&);

    // histogram declaration
    edm::Service<TFileService> fs;
    const reco::GenParticle* mother;

    TH1F* h_NPV;
    //After Lepton Selection
    TH1F* h_ALS_etaLMuMu;
    TH1F* h_ALS_ptLMuMu;
    TH1F* h_ALS_etaLElEl;
    TH1F* h_ALS_ptLElEl;
    TH1F* h_ALS_etaLElMu;
    TH1F* h_ALS_ptLElMu;
    TH1F* h_ALS_etaLDiLep;
    TH1F* h_ALS_ptLDiLep;
    //After MET Requirements
    TH1F* h_AMS_ptLepMuMu;
    TH1F* h_AMS_ptLepElMu;
    TH1F* h_AMS_ptLepElEl;
    TH1F* h_AMS_ptLepDiLep;
    TH1F* h_AMS_mLepMuMu;
    TH1F* h_AMS_mLepElEl;
    TH1F* h_AMS_mLepElMu;
    TH1F* h_AMS_mLepDiLep;
    TH1F* h_AMS_METMuMu;
    TH1F* h_AMS_METElEl;
    TH1F* h_AMS_METElMu;
    TH1F* h_AMS_METDiLep;
    //After jet Selectin
    TH1F* h_AJS_etaLepMuMu;
    TH1F* h_AJS_ptLepMuMu;
    TH1F* h_AJS_etaLepElEl;
    TH1F* h_AJS_ptLepElEl;
    TH1F* h_AJS_etaLepElMu;
    TH1F* h_AJS_ptLepElMu;
    TH1F* h_AJS_etaLepDiLep;
    TH1F* h_AJS_ptLepDiLep;
    //After Btag requirement
    TH1F* h_ABS_etaLepMuMu;
    TH1F* h_ABS_ptLepMuMu;
    TH1F* h_ABS_etaLepElEl;
    TH1F* h_ABS_ptLepElEl;
    TH1F* h_ABS_etaLepElMu;
    TH1F* h_ABS_ptLepElMu;
    TH1F* h_ABS_etaLepDiLep;
    TH1F* h_ABS_ptLepDiLep;
    TH1F* h_ABS_mLepMuMu;
    TH1F* h_ABS_mLepElEl;
    TH1F* h_ABS_mLepElMu;
    TH1F* h_ABS_mLepDiLep;

    TH1F* h_ABS_mLepNoVetoMuMu;
    TH1F* h_ABS_mLepNoVetoElEl;
    TH1F* h_ABS_mLepNoVetoElMu;
    TH1F* h_ABS_mLepNoVetoDiLep;
    TH1F* h_AoneBS_mLepNoVetoMuMu;
    TH1F* h_AoneBS_mLepNoVetoElEl;
    TH1F* h_AoneBS_mLepNoVetoElMu;
    TH1F* h_ALS_mLepNoVetoMuMu;
    TH1F* h_ALS_mLepNoVetoElEl;
    TH1F* h_ALS_mLepNoVetoElMu;
    TH1F* h_AJS_mLepNoVetoMuMu;
    TH1F* h_AJS_mLepNoVetoElEl;
    TH1F* h_AJS_mLepNoVetoElMu;
    TH1F* h_AMS_mLepNoVetoMuMu;
    TH1F* h_AMS_mLepNoVetoElEl;
    TH1F* h_AMS_mLepNoVetoElMu;
    TH1F* h_ATS_mLepNoVetoMuMu;
    TH1F* h_ATS_mLepNoVetoElEl;
    TH1F* h_ATS_mLepNoVetoElMu;
    TH1F* h_ATS_mLepNoVetoDiLep;

    TH1F* h_ABS_METMuMu;
    TH1F* h_ABS_METElEl;
    TH1F* h_ABS_METElMu;
    TH1F* h_ABS_METDiLep;
    TH1F* h_ABS_NJetsMuMu;
    TH1F* h_ABS_NJetsElMu;
    TH1F* h_ABS_NJetsElEl;
    TH1F* h_ABS_NJetsDiLep;
    TH1F* h_ABS_NBJetsMuMu;
    TH1F* h_ABS_NBJetsElEl;
    TH1F* h_ABS_NBJetsElMu;
    TH1F* h_ABS_NBJetsDiLep;
    TH1F* h_ABS_ptLeadingJetMuMu;
    TH1F* h_ABS_ptLeadingJetElEl;
    TH1F* h_ABS_ptLeadingJetElMu;
    TH1F* h_ABS_ptLeadingJetDiLep;
    TH1F* h_ABS_etaLeadingJetMuMu;
    TH1F* h_ABS_etaLeadingJetElEl;
    TH1F* h_ABS_etaLeadingJetElMu;
    TH1F* h_ABS_etaLeadingJetDiLep;

    //After tt reconstruction
    TH1F* h_NBJetsMuMu;
    TH1F* h_NBJetsElMu;
    TH1F* h_NBJetsElEl;
    TH1F* h_NBJetsDiLep;
    TH1F* h_cosMuMu;
    TH1F* h_cosElEl;
    TH1F* h_cosElMu;
    TH1F* h_cosDiLep;
    TH1F* h_cosGenMuMu;
    TH1F* h_cosGenElEl;
    TH1F* h_cosGenElMu;
    TH1F* h_cosGen;
    TH1F* h_mTTbarMuMu;
    TH1F* h_mTTbarElEl;
    TH1F* h_mTTbarElMu;
    TH1F* h_TTbarM;
    TH1F* h_GenTTbarM;
    TH1F* h_ptTMuMu;
    TH1F* h_ptTElEl;
    TH1F* h_ptTElMu;
    TH1F* h_ptTDiLep;
    TH1F* h_yTMuMu;
    TH1F* h_yTElEl;
    TH1F* h_yTElMu;
    TH1F* h_yTDiLep;
    TH1F* h_ptWMuMu;
    TH1F* h_ptWElEl;
    TH1F* h_ptWElMu;
    TH1F* h_ptWDiLep;
    TH1F* h_yWMuMu;
    TH1F* h_yWElEl;
    TH1F* h_yWElMu;
    TH1F* h_yWDiLep;
    TH1F* h_PtMu;
    TH1F* h_etaMu;
    TH2F* h_truthRecoMuMu;
    TH2F* h_truthRecoElEl;
    TH2F* h_truthRecoElMu;

    TH1F* h_ATS_etaLepMuMu;
    TH1F* h_ATS_ptLepMuMu;
    TH1F* h_ATS_etaLepElEl;
    TH1F* h_ATS_ptLepElEl;
    TH1F* h_ATS_etaLepElMu;
    TH1F* h_ATS_ptLepElMu;

    TH1F* h_ATS_mLepMuMu;
    TH1F* h_ATS_mLepElEl;
    TH1F* h_ATS_mLepElMu;


    TH1F* h_ATS_METMuMu;
    TH1F* h_ATS_METElEl;
    TH1F* h_ATS_METElMu;
    TH1F* h_ATS_NJetsMuMu;
    TH1F* h_ATS_NJetsElMu;
    TH1F* h_ATS_NJetsElEl;
    TH1F* h_ATS_NBJetsMuMu;
    TH1F* h_ATS_NBJetsElEl;
    TH1F* h_ATS_NBJetsElMu;


    TH2F* h_truthRecoCos;
    TH2F* h_truthRecoMTT;
    TH1F* h_RecoMinusTruthCos;
    TH1F* h_RecoMinusTruthMTT;
    TH1F* h_Nevents;
    TH1F* h_Nevents_weighted;
    TH1F* h_Nevents_AVS;
    TH1F* h_Nevents_ALS;
    TH1F* h_NeventsMuMu_ALS;
    TH1F* h_NeventsElMu_ALS;
    TH1F* h_NeventsElEl_ALS;
    TH1F* h_Nevents_AT;
    TH1F* h_Nevents_AJS;
    TH1F* h_NeventsMuMu_AJS;
    TH1F* h_NeventsElEl_AJS;
    TH1F* h_NeventsElMu_AJS;
    TH1F* h_Nevents_ABS;
    TH1F* h_NeventsMuMu_ABS;
    TH1F* h_NeventsElEl_ABS;
    TH1F* h_NeventsElMu_ABS;
    TH1F* h_Nevents_AMS;
    TH1F* h_NeventsMuMu_AMS;
    TH1F* h_NeventsElEl_AMS;
    TH1F* h_NeventsElMu_AMS;
    TH1F* h_Weight;
    TH1F* h_Nevents_DiMu;
    TH1F* h_Nevents_DiEl;
    TH1F* h_Nevents_ElMu;
    TH1F* h_Nevents_top;
    TH1F* h_NeventsMuMu_top;
    TH1F* h_NeventsElMu_top;
    TH1F* h_NeventsElEl_top;
    TH1F* h_Nevents_ATS;
    TH1F* h_NeventsMuMu_ATS;
    TH1F* h_NeventsElEl_ATS;
    TH1F* h_NeventsElMu_ATS;


    TH1F* h_NinMuMu_ALS;
    TH1F* h_NoutMuMu_ALS;
    TH1F* h_NinElMu_ALS;
    TH1F* h_NoutElMu_ALS;
    TH1F* h_NinElEl_ALS;
    TH1F* h_NoutElEl_ALS;
    TH1F* h_NinMuMu_AJS;
    TH1F* h_NoutMuMu_AJS;
    TH1F* h_NinElMu_AJS;
    TH1F* h_NoutElMu_AJS;
    TH1F* h_NinElEl_AJS;
    TH1F* h_NoutElEl_AJS;
    TH1F* h_NinMuMu_AMS;
    TH1F* h_NoutMuMu_AMS;
    TH1F* h_NinElMu_AMS;
    TH1F* h_NoutElMu_AMS;
    TH1F* h_NinElEl_AMS;
    TH1F* h_NoutElEl_AMS;
    TH1F* h_NinMuMu_ABS;
    TH1F* h_NoutMuMu_ABS;
    TH1F* h_NinElMu_ABS;
    TH1F* h_NoutElMu_ABS;
    TH1F* h_NinElEl_ABS;
    TH1F* h_NoutElEl_ABS;
    TH1F* h_NinNoMET_ABS;
    TH1F* h_NoutNoMET_ABS;
    TH1F* h_NinMuMu_NoMET_ABS;
    TH1F* h_NoutMuMu_NoMET_ABS;
    TH1F* h_NinElMu_NoMET_ABS;
    TH1F* h_NoutElMu_NoMET_ABS;
    TH1F* h_NinElEl_NoMET_ABS;
    TH1F* h_NoutElEl_NoMET_ABS;
    TH1F* h_NinMuMu_AoneBS;
    TH1F* h_NoutMuMu_AoneBS;
    TH1F* h_NinElMu_AoneBS;
    TH1F* h_NoutElMu_AoneBS;
    TH1F* h_NinElEl_AoneBS;
    TH1F* h_NoutElEl_AoneBS;
    TH1F* h_NinMuMu_NoMET_AoneBS;
    TH1F* h_NoutMuMu_NoMET_AoneBS;
    TH1F* h_NinElMu_NoMET_AoneBS;
    TH1F* h_NoutElMu_NoMET_AoneBS;
    TH1F* h_NinElEl_NoMET_AoneBS;
    TH1F* h_NoutElEl_NoMET_AoneBS;

    TH1F* h_NinMuMu_ATS;
    TH1F* h_NoutMuMu_ATS;
    TH1F* h_NinElMu_ATS;
    TH1F* h_NoutElMu_ATS;
    TH1F* h_NinElEl_ATS;
    TH1F* h_NoutElEl_ATS;
    TH1F* h_NinMuMu_NoMET_ATS;
    TH1F* h_NoutMuMu_NoMET_ATS;
    TH1F* h_NinElMu_NoMET_ATS;
    TH1F* h_NoutElMu_NoMET_ATS;
    TH1F* h_NinElEl_NoMET_ATS;
    TH1F* h_NoutElEl_NoMET_ATS;

    double coriso= 9999; // initialise to dummy value
    double coriso2= 9999; // initialise to dummy value





    //    TtFullLepKinSolver* amwtsolver;
    TtAMWTSolver* amwtSolver;
    EffectiveAreas effectiveAreas_;
    edm::EDGetTokenT<TtGenEvent> ttgenEvt_;
    bool isData;
    string ptRes;
    string phiRes;
    string sfRes;
    string ptResData;
    string phiResData;
    string sfResData;
    int NEvent=0;
    string outfileName;
    //TTree defenition for storing important stuff
    TFile* f_outFile;
    TTree* t_outTree;
    vector<Double_t> RecoCos;
    vector<Double_t> TruthCos;
    vector<Double_t> RecoCosMuMu;
    vector<Double_t> TruthCosMuMu;
    vector<Double_t> RecoMTT;
    vector<Double_t> TruthMTT;
    vector<Double_t> TruthMTTMuMu;
    vector<Double_t> RecoMTTMuMu;
    int nGoodVtxs = 0;

    edm::EDGetTokenT<edm::TriggerResults> triggerResluts_;

    //    int n_afterVertex = 0;
    //    int n_afterHLT = 0;
    //    int n_afterDiLepton = 0;
    //    int n_afterDiMu = 0;
    //    int n_afterDiEl = 0;
    //    int n_afterElMu = 0;
    //    int n_DiMuHLT = 0;
    //    int n_DiElHLT = 0;
    //    int n_ElMuHLT = 0;
    //    int n_afterMet = 0;
    //    int n_after2Jets = 0;
    //    int n_after2BJets = 0;
    //    int n_afterTop = 0;
    edm::View<pat::PackedGenParticle> genColl;
    edm::EDGetTokenT<LHERunInfoProduct> lheInfo_;

    TLorentzVector PosLep;
    TLorentzVector NegLep;
    TLorentzVector DiLep;
    bool b_Mu1=0,b_Mu2=0,b_Mu3=0,b_Mu4=0,b_MuMu1=0,b_MuMu2=0,b_MuMu3=0,b_MuMu4=0,b_El1=0,b_El2=0,b_El3=0,b_ElEl1=0,b_ElEl2=0,b_ElMu1=0,b_ElMu2=0,b_ElMu3=0,b_ElMu4=0,b_ElMu5=0,b_ElMu6=0;
    bool isPythia;
    edm::EDGetTokenT<edm::TriggerResults> triggerFilters_;
    bool Flag_goodVertices=0, Flag_globalTightHalo2016Filter=0,Flag_HBHENoiseFilter=0,Flag_HBHENoiseIsoFilter=0,Flag_BadPFMuonFilter=0,Flag_BadChargedCandidateFilter=0,Flag_EcalDeadCellTriggerPrimitiveFilter=0,Flag_eeBadScFilter=0;
    bool DiMu;
    bool DiEl;
    bool ElMu;
    string egammaSF;
    TFile* f_egammaSF;
    float SF_posEl = 1;
    float SF_negEl = 1;
    float SF_posMu = 1;
    float SF_negMu = 1;
    bool isRunGH = false;
    string muonISOSF;
    string muonIDSF;
    TFile* f_muonISOSF;
    TFile* f_muonIDSF;
    string egammaTkSF;
    TFile* f_egammaTkSF;
    string btagSf;
    BTagCalibrationReader reader;
    string muonTkSF;
    TFile* f_muonTkSF;
    bool isSingleMuon;
    bool isSingleElectron;
    bool isRunH;
    // ID decisions objects
    edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_;

    string ee_sf;
    TFile* f_ee_sf;
    string em_sf;
    TFile* f_em_sf;
    string mm_sf;
    TFile* f_mm_sf;

    vector<double> t_cos;

    vector<myLorentzVector> t_Leptons;
    vector<double> t_Leptons_pt;
    vector<double> t_Leptons_eta;
    vector<double> t_Leptons_phi;
    vector<double> t_Leptons_e;
    vector<int> t_Leptons_charge;
    vector<string> t_Leptons_type;
    vector<double> t_Leptons_etaSC;
    vector<double> t_Leptons_d0;
    vector<double> t_Leptons_dz;

    vector<double> t_gen_Leptons_pt;
    vector<double> t_gen_Leptons_eta;
    vector<double> t_gen_Leptons_phi;
    vector<double> t_gen_Leptons_e;
    vector<int> t_gen_Leptons_id;
    vector<int> t_gen_Leptons_charge;

    vector<double> t_gen_daughters_pt;
    vector<double> t_gen_daughters_eta;
    vector<double> t_gen_daughters_phi;
    vector<double> t_gen_daughters_e;
    vector<int> t_gen_daughters_id;

    vector<double> t_gen_mothers_pt;
    vector<double> t_gen_mothers_eta;
    vector<double> t_gen_mothers_phi;
    vector<double> t_gen_mothers_e;
    vector<int> t_gen_mothers_id;


    vector<myLorentzVector> t_Jets;
    vector<double> t_Jets_pt;
    vector<double> t_Jets_eta;
    vector<double> t_Jets_phi;
    vector<double> t_Jets_e;
    vector<int> t_Jets_hadflav;

    vector<bool> t_Jets_LBdiscriminator;
    vector<bool> t_Jets_MBdiscriminator;
    vector<bool> t_Jets_TBdiscriminator;
    vector<double> t_Jets_dis;

    vector<double> t_Jets_pt_resolution;
    vector<double> t_Jets_phi_resolution;

    vector<double> t_Jets_sf_nominal;
    vector<double> t_Jets_sf_up;
    vector<double> t_Jets_sf_down;


    vector<double> t_gen_Jets_pt;
    vector<double> t_gen_Jets_eta;
    vector<double> t_gen_Jets_phi;
    vector<double> t_gen_Jets_e;
    vector<double> t_gen_Jets_hadflav;
    vector<double> t_gen_Jets_prtnflav;

    vector<double> t_Met_pt;
    vector<double> t_Met_px;
    vector<double> t_Met_py;
    vector<double> t_Met_sum;

    vector<myLorentzVector> t_bJets;
    vector<double> t_bJets_pt;
    vector<double> t_bJets_eta;
    vector<double> t_bJets_phi;
    vector<double> t_bJets_e;
    int t_Leptons_size;
    int t_Muons_size;
    int t_Electrons_size;
    vector<myLorentzVector> t_top;
    vector<double> t_top_pt;
    vector<double> t_top_eta;
    vector<double> t_top_phi;
    vector<double> t_top_e;

    vector<myLorentzVector> t_antitop;
vector<double> t_antitop_pt;
    vector<double> t_antitop_eta;
    vector<double> t_antitop_phi;
    vector<double> t_antitop_e;

    vector<myLorentzVector> t_W;
vector<double> t_W_pt;
    vector<double> t_W_eta;
    vector<double> t_W_phi;
    vector<double> t_W_e;
    vector<myLorentzVector> t_antiW;
    vector<double> t_antiW_pt;
    vector<double> t_antiW_eta;
    vector<double> t_antiW_phi;
    vector<double> t_antiW_e;
    vector<myLorentzVector> t_nu;
    vector<double> t_nu_pt;
    vector<double> t_nu_eta;
    vector<double> t_nu_phi;
    vector<double> t_nu_e;
    vector<myLorentzVector> t_antinu;
    vector<double> t_antinu_pt;
    vector<double> t_antinu_eta;
    vector<double> t_antinu_phi;
    vector<double> t_antinu_e;
    
    int t_Run;
    int t_Event;
    int t_bunch;
    int t_lumi;
    vector<double> w_bjets;
    double w_bjet;
    double w_posEl;
    double w_negEl;
    double w_posMu;
    double w_negMu;
    double w_ee;
    double w_em;
    double w_mm;
    double w_top;
    double w_mc;
    float t_rho;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSrc_;
    int t_num_PU_vertices;
    int t_PU_BunchCrossing;
    int t_num_PU_gen_vertices;
    int t_num_PV;
    double t_PU_weight;
    double t_PU_weight_secondWay;
    edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
    JME::JetResolution* ptResol;
    JME::JetResolution* phiResol;
    JME::JetResolutionScaleFactor* jerSF;
    edm::LumiReWeighting LumiWeights_;
    string s_pileup_data;
    string s_pileup_mc;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleIdFullInfoMapToken_;
    edm::EDGetTokenT<bool> BadChCandFilterToken_;
    edm::EDGetTokenT<bool> BadPFMuonFilterToken_;


     bool verboseIdFlag_;

};

