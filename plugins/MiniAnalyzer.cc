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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h"
#include "ttamwtsolver.h"
#include "MSETools.h"




#include "DileptonAnalyticalSolver.h"

using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer  {
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


    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
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
    TH2F* h_truthRecoCos;
    TH2F* h_truthRecoMTT;
    TH1F* h_truthMinusRecoCos;
    TH1F* h_truthMinusRecoMTT;

    double coriso= 999; // initialise to dummy value
    double coriso2= 999; // initialise to dummy value
    double DiLepMass;




    //    TtFullLepKinSolver* amwtsolver;
    TtAMWTSolver* amwtSolver;
    EffectiveAreas effectiveAreas_;
    edm::EDGetTokenT<TtGenEvent> ttgenEvt_;
    bool isData;
    string ptRes;
    string phiRes;
    string sfRes;

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

    int n_afterVertex = 0;
    int n_afterHLT = 0;
    int n_afterDiLepton = 0;
    int n_afterDiMu = 0;
    int n_afterDiEl = 0;
    int n_afterElMu = 0;
    int n_DiMuHLT = 0;
    int n_DiElHLT = 0;
    int n_ElMuHLT = 0;
    int n_afterMet = 0;
    int n_after2Jets = 0;
    int n_after2BJets = 0;
    int n_afterTop = 0;
    edm::View<pat::PackedGenParticle> genColl;
    edm::EDGetTokenT<LHERunInfoProduct> lheInfo_;

    TLorentzVector PosLep;
    TLorentzVector NegLep;
    TLorentzVector DiLep;
    std::vector<double> *SlimmedMuon_Pt;
    bool b_Mu1=0,b_Mu2=0,b_ElEl1=0,b_ElEl2=0,b_MuEl1=0,b_MuEl2=0,b_ElMu1=0,b_ElMu2=0,b_ElMu3=0,b_ElMu4=0;



};


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    prunedGenToken_(mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken_(mayConsume<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    genEvtInfo_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfo"))),
    lheEvtInfo_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvtInfo"))),
    rho_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
    effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() ),
    ttgenEvt_( mayConsume<TtGenEvent>(iConfig.getParameter<edm::InputTag>("ttgen"))),
    isData((iConfig.getParameter<bool>("isData"))),
    ptRes( (iConfig.getParameter<string>("ptRes"))),
    phiRes( (iConfig.getParameter<string>("phiRes")) ),
    sfRes( (iConfig.getParameter<string>("sfRes")) ),
    outfileName((iConfig.getParameter<string>("outFileName"))),
    triggerResluts_(mayConsume<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
    lheInfo_(mayConsume<LHERunInfoProduct,InRun>(iConfig.getParameter<edm::InputTag>("externalLHEProducer")))

{
    // initializing the solver
    amwtSolver = new TtAMWTSolver(isData,172.44,172.50,0.1,80.4,4.8,ptRes,phiRes,sfRes);

    //    std::vector<double> nupars= {30.7137,56.2880,23.0744,59.1015,24.9145};
    //    amwtsolver = new TtFullLepKinSolver(172.5,172.5,10,nupars,80.4,4.8);

    // Initializing  output root file;
    f_outFile = TFile::Open(outfileName.c_str(),"RECREATE");
    // cout << outfileName.c_str() << endl;


}

MiniAnalyzer::~MiniAnalyzer()
{
    delete f_outFile;
}

void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;


    /////////////FILLINGINPUTS/////////////////////
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);

    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);

    Handle<edm::View<reco::GenParticle> > pruned;
    Handle<edm::View<pat::PackedGenParticle> > packed;
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    edm::Handle<LHEEventProduct> lhEvtInfo;
    edm::Handle<TtGenEvent> genEvent;
    if(!isData){
        iEvent.getByToken(prunedGenToken_,pruned);
        iEvent.getByToken(packedGenToken_,packed);
        iEvent.getByToken( genEvtInfo_, genEvtInfo );
        iEvent.getByToken( lheEvtInfo_, lhEvtInfo);
        iEvent.getByToken(ttgenEvt_, genEvent);
        genColl = *packed;

    }
    
    edm::Handle<double> rhoHandle_;
    iEvent.getByToken( rho_,rhoHandle_);
    float rho = *rhoHandle_;

    amwtSolver->SetRho(rho);

    isPosMu = false;
    isNegMu = false;
    isPosEl = false;
    isNegEl = false;
    isDiMuon = false;
    isDiElectron = false;
    isElMu = false;
    isMuEl = false;
    isDiLeptonic = false;

    /////////////////////HLT
    ///
    edm::Handle<edm::TriggerResults> trigResults; //our trigger result object

    iEvent.getByToken(triggerResluts_,trigResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
    

    std::string Mu1,Mu2,ElEl1,ElEl2,ElMu1,ElMu2,ElMu3,ElMu4,MuEl1,MuEl2;
    if(isData){
        // Double Muon  36.811 fb^-1 Mu1 or Mu2
        Mu1="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
        Mu2="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";

        // Double Electron 36.615 fb^-1
        ElEl1="HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v";

        // Double Electron 36.811 fb^-1
        ElEl2="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

        // MuEl 36.811 fb^-1 MuEl1 OR MuEl2
        MuEl1="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
        MuEl2="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

        // El Mu  36.810 ELMu1 OR ElMu2 OR ElMu3
        ElMu1="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
        ElMu2="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v";
        ElMu3="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
        ElMu4="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v";
        //cout << "safe here" << endl;
    }
    else{
        Mu1="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7";
        Mu2="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6";

        ElEl1="HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v8";

        ElEl2="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9";

        MuEl1 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9";
        MuEl2 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4";

        ElMu1="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9";
        ElMu2="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3";
        ElMu3="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4";
        ElMu4="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9";
        // cout << "why here" << endl;


    }
    for (unsigned int i = 0, n = trigResults->size(); i < n; ++i) {
        std::string nameHLT,st;
        nameHLT = trigNames.triggerName(i);
        st = nameHLT.substr(0, nameHLT.size()-1);
        if(st.compare(Mu1) == 0)
        {
            b_Mu1 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(Mu2) == 0)
        {
            b_Mu2 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(ElEl2) == 0)
        {
            b_ElEl2 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(MuEl1) == 0)
        {
            b_MuEl1 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(MuEl2) == 0)
        {
            b_MuEl2 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(ElMu1) == 0)
        {
            b_ElMu1 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(ElMu2) == 0)
        {
            b_ElMu2 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(ElMu3) == 0)
        {
            b_ElMu3 = trigResults->accept(i);
            cout << st << endl;
        }
        else if(st.compare(ElMu4) == 0)
        {
            b_ElMu4 = trigResults->accept(i);
            cout << st << endl;
        }

        //        std::cout << "Trigger " << trigNames.triggerName(i) << "i: "<<i

        //                  <<": " << (trigResults->accept(i) ? "PASS" : "fail (or not run)")
        //                 << std::endl;
    }

    if(!isData){
        b_Mu1=trigResults->accept(trigNames.triggerIndex(Mu1));

        b_Mu2=trigResults->accept(trigNames.triggerIndex(Mu2));

        //    bool b_ElEl1=trigResults->accept(trigNames.triggerIndex(ElEl1));
        b_ElEl2=trigResults->accept(trigNames.triggerIndex(ElEl2));

        b_MuEl1=trigResults->accept(trigNames.triggerIndex(MuEl1));
        b_MuEl2=trigResults->accept(trigNames.triggerIndex(MuEl2));

        b_ElMu1=trigResults->accept(trigNames.triggerIndex(ElMu1));
        b_ElMu2=trigResults->accept(trigNames.triggerIndex(ElMu2));
        b_ElMu3=trigResults->accept(trigNames.triggerIndex(ElMu3));
        b_ElMu4=trigResults->accept(trigNames.triggerIndex(ElMu4));
    }

    //    cout << b_Mu1 << "      "<< b_Mu2 << endl;
    
    /////////////////////////////////////////////////
    // make sure we have a good vertex //////////////
    /////////////////////////////////////////////////
    ++NEvent;
    //    // cout << "number of Events " << NEvent << endl;
    if (vertices->empty()) return;
    VertexCollection::const_iterator PV = vertices->end();
    for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
        if ( !(vtx->isFake())
             && vtx->ndof() >= 4. && vtx->position().Rho() < 2.0
             && fabs(vtx->position().Z()) < 24.0) {
            PV = vtx;
            break;
        }
    }
    if (PV==vertices->end()) return;
    ++n_afterVertex;
    // count how many good vertices we have
    nGoodVtxs = 0;
    for (VertexCollection::const_iterator vtx = vertices->begin();vtx != vertices->end(); ++vtx) {
        if ( !(vtx->isFake()) && vtx->ndof() >= 4. && vtx->position().Rho() <= 2.0 && fabs(vtx->position().Z()) <= 24.) nGoodVtxs++;
    }
    //////////////EVENT WEIGHT/////////
    /// finding weight so we can
    /// compare different backgrounds
    double theWeight;
    if(isData)
    {
        theWeight = 1;
    }
    else
    {
        theWeight = genEvtInfo->weight();
        int whichWeight = 1;
        theWeight *= lhEvtInfo->weights()[whichWeight].wgt/lhEvtInfo->originalXWGTUP();

    }

    h_NPV->Fill(nGoodVtxs,theWeight);


    //////////////MUON////////////////
    /// finding pos and neg muon with highest pt
    ///
    pat::Muon posMu;
    pat::Muon negMu;


    for (pat::MuonCollection::const_iterator mup = muons->begin();mup != muons->end(); ++mup){
        SlimmedMuon_Pt->push_back(mup->pt());

        if( !(mup->charge() > 0)) continue;
        if( !(mup->pt() > 20.0 )) continue;
        if( !(fabs(mup->eta()) < 2.4 )) continue;
        if( !(mup->isPFMuon())) continue;
        if( !(mup->isGlobalMuon() || mup->isTrackerMuon())) continue;
        if( mup->isIsolationValid()){
            reco::MuonPFIsolation pfR04 = mup->pfIsolationR04();
            coriso = pfR04.sumChargedHadronPt + std::max(0., pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
        }
        if (!(coriso/mup->pt() <  0.15)) continue;
        if(mup->pt() > posMu.pt() ) posMu = *mup;
        //        // cout << mup->pt() <<endl;
        //        // cout << posMu.pt() << endl;
        h_PtMu->Fill(mup->pt(),theWeight);
        h_ALS_ptLMuMu->Fill(posMu.pt(),theWeight);
        h_ALS_etaLMuMu->Fill(posMu.eta(),theWeight);
        h_etaMu->Fill(mup->eta(),theWeight);
        isPosMu = true;

    }

    for (pat::MuonCollection::const_iterator mum = muons->begin();mum != muons->end(); ++mum){
        if( !(mum->charge() < 0)) continue;
        if( !(mum->pt() > 20.0 )) continue;
        if( !(fabs(mum->eta()) < 2.4 )) continue;
        if( !(mum->isPFMuon())) continue;
        if( !(mum->isGlobalMuon() || mum->isTrackerMuon())) continue;
        if( mum->isIsolationValid()){
            reco::MuonPFIsolation pfR04 = mum->pfIsolationR04();
            coriso2 = pfR04.sumChargedHadronPt + std::max(0., pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
        }
        if (!(coriso2/mum->pt() <  0.12)) continue;
        if(mum->pt() > negMu.pt() ) negMu = *mum;
        h_PtMu->Fill(mum->pt(),theWeight);
        h_ALS_etaLMuMu->Fill(negMu.eta(),theWeight);
        h_ALS_ptLMuMu->Fill(negMu.pt(),theWeight);
        h_etaMu->Fill(mum->eta(),theWeight);
        isNegMu = true;


    }


    //    // cout << posMu.charge() << " charges " << negMu.charge() << endl;
    ////////////////////ELECTRONS//////////////////////////
    /// electron identification
    ///
    pat::Electron posEl;
    pat::Electron negEl;

    for (pat::ElectronCollection::const_iterator elp = electrons->begin();elp != electrons->end();++elp){
        if( !( elp->charge() > 0 ) ) continue;
        if( !( elp->pt() > 20 ) ) continue;
        if( !( fabs(elp->eta() ) < 2.5)) continue;
        if( !( elp->gsfTrack()->dz() < 0.04 )) continue;
        if( !(elp->passConversionVeto())) continue;
        if( !(elp->gsfTrack()->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) <= 0)) continue;
        GsfElectron::PflowIsolationVariables pfIso = elp->pfIsolationVariables();
        static double relCombIsoEA = 999.;
        float abseta = fabs(elp->superCluster()->eta());
        float eA = effectiveAreas_.getEffectiveArea(abseta);
        relCombIsoEA = (( pfIso.sumChargedHadronPt
                          + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho) )
                        / elp->pt() );
        if( !(relCombIsoEA < 0.15)) continue;
        if( !(elp->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") > 0.5 )) continue;
        //        // cout << elp->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") << "MVA Id "<< endl;
        double dR =0;
        for(pat::MuonCollection::const_iterator dumMu = muons->begin();dumMu != muons->end();++dumMu){
            if(!(dumMu->isGlobalMuon())) continue;
            dR = ROOT::Math::VectorUtil::DeltaR(dumMu->p4(),elp->p4());
            if(!(dR > 0.1)) break;
        }
        if(!(dR > 0.1)) continue;
        if( elp->pt() > posEl.pt() ) posEl = *elp;
        isPosEl = true;

    }

    for (pat::ElectronCollection::const_iterator elm = electrons->begin();elm != electrons->end();++elm){
        if( !( elm->charge() < 0 ) ) continue;
        if( !( elm->pt() > 20 ) ) continue;
        if( !( fabs(elm->eta() ) < 2.5)) continue;
        if( !( elm->gsfTrack()->dz() < 0.04 )) continue;
        if( !(elm->passConversionVeto())) continue;
        if( !(elm->gsfTrack()->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) <= 0)) continue;
        GsfElectron::PflowIsolationVariables pfIso = elm->pfIsolationVariables();
        static double relCombIsoEA = 999.;
        float abseta = fabs(elm->superCluster()->eta());
        float eA = effectiveAreas_.getEffectiveArea(abseta);
        relCombIsoEA = (( pfIso.sumChargedHadronPt
                          + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho) )
                        / elm->pt() );
        if( !(relCombIsoEA < 0.15)) continue;
        if( !(elm->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") > 0.5 )) continue;
        // cout << elm->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") << "MVA Id "<< endl;
        double dR =0;
        for(pat::MuonCollection::const_iterator dumMu = muons->begin();dumMu != muons->end();++dumMu){
            if(!(dumMu->isGlobalMuon())) continue;
            dR = ROOT::Math::VectorUtil::DeltaR(dumMu->p4(),elm->p4());
            if(!(dR > 0.1)) break;
        }
        if(!(dR > 0.1)) continue;
        if( elm->pt() > negEl.pt() ) negEl = *elm;
        isNegEl = true;

    }

    //After Lepton Selection
    if(isPosMu && isNegMu) isDiMuon = true;
    if(isDiMuon)++n_afterDiMu;
    if(isPosEl && isNegEl) isDiElectron = true;
    if(isDiElectron)++n_afterDiEl;
    if( isPosEl && isNegMu)  isElMu = true;
    if( isNegEl && isPosMu ) isMuEl = true;
    if(isElMu || isMuEl) ++n_afterElMu;
    if(isDiMuon || isDiElectron || isElMu ) isDiLeptonic = true;
    if(!isDiLeptonic)return;
    ++n_afterDiLepton;
    if(isDiMuon && !(b_Mu1 || b_Mu2) ) return;
    if(isDiElectron && !(b_ElEl2)) return;
    if(isElMu && !(b_ElMu1 || b_ElMu2 || b_ElMu3 || b_ElMu4 || b_MuEl1 || b_MuEl2)) return;
    if(isMuEl && !(b_ElMu1 || b_ElMu2 || b_ElMu3 || b_ElMu4 || b_MuEl1 || b_MuEl2)) return;
    ++n_afterHLT;


    if(isDiMuon)
    {
        PosLep.SetPtEtaPhiM(posMu.pt(),posMu.eta(),posMu.phi(),posMu.mass());
        NegLep.SetPtEtaPhiM(negMu.pt(),negMu.eta(),negMu.phi(),negMu.mass());
        DiLep = PosLep + NegLep;
        if((DiLep.M() < 106 && DiLep.M() > 76)) return;

        h_ALS_etaLMuMu->Fill(posMu.eta(),theWeight);
        h_ALS_etaLMuMu->Fill(negMu.eta(),theWeight);
        h_ALS_ptLMuMu->Fill(posMu.pt(),theWeight);
        h_ALS_ptLMuMu->Fill(negMu.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posMu.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negMu.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posMu.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negMu.pt(),theWeight);

    }
    if(isDiElectron)
    {
        PosLep.SetPtEtaPhiM(posEl.pt(),posEl.eta(),posEl.phi(),posEl.mass());
        NegLep.SetPtEtaPhiM(negEl.pt(),negEl.eta(),negEl.phi(),negEl.mass());
        DiLep = PosLep + NegLep;
        if((DiLep.M() < 106 && DiLep.M() > 76)) return;
        h_ALS_etaLElEl->Fill(posEl.eta(),theWeight);
        h_ALS_etaLElEl->Fill(negEl.eta(),theWeight);
        h_ALS_ptLElEl->Fill(posEl.pt(),theWeight);
        h_ALS_ptLElEl->Fill(negEl.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posEl.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negEl.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posEl.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negEl.pt(),theWeight);

    }
    if(isElMu)
    {
        PosLep.SetPtEtaPhiM(posEl.pt(),posEl.eta(),posEl.phi(),posEl.mass());
        NegLep.SetPtEtaPhiM(negMu.pt(),negMu.eta(),negMu.phi(),negMu.mass());
        DiLep = PosLep + NegLep;
        if((DiLep.M() < 106 && DiLep.M() > 76)) return;
        h_ALS_etaLElMu->Fill(posEl.eta(),theWeight);
        h_ALS_etaLElMu->Fill(negMu.eta(),theWeight);
        h_ALS_ptLElMu->Fill(posEl.pt(),theWeight);
        h_ALS_ptLElMu->Fill(negMu.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posEl.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negMu.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posEl.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negMu.pt(),theWeight);

    }
    if(isMuEl)
    {
        PosLep.SetPtEtaPhiM(posMu.pt(),posMu.eta(),posMu.phi(),posMu.mass());
        NegLep.SetPtEtaPhiM(negEl.pt(),negEl.eta(),negEl.phi(),negEl.mass());
        DiLep = PosLep + NegLep;
        if((DiLep.M() < 106 && DiLep.M() > 76)) return;
        h_ALS_etaLElMu->Fill(posMu.eta(),theWeight);
        h_ALS_etaLElMu->Fill(negEl.eta(),theWeight);
        h_ALS_ptLElMu->Fill(posMu.pt(),theWeight);
        h_ALS_ptLElMu->Fill(negEl.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posMu.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negEl.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posMu.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negEl.pt(),theWeight);

    }



    ///////////////////JETS//////////////////////////////
    vector<pat::Jet> bjets;
    vector<pat::Jet> njets;
    //    pat::Jet j;j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    for(pat::JetCollection::const_iterator jet_it=jets->begin(); jet_it != jets->end();++jet_it){
        if( !( jet_it->pt() > 30 )) continue;
        if( !( fabs(jet_it->eta()) < 2.4 )) continue;

        if( !(jet_it->neutralHadronEnergyFraction() < 0.99 && jet_it->neutralEmEnergyFraction() < 0.99 && (jet_it->chargedMultiplicity() + jet_it->neutralMultiplicity())> 1.
              && jet_it->chargedHadronEnergyFraction() > 0. && jet_it->chargedEmEnergyFraction() < 0.99 && jet_it->chargedMultiplicity() > 1.)) continue;
        njets.push_back(*jet_it);
        if( !(jet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") >  0.460 )) continue;          //loose working point for btaging new value is 0.605
        bjets.push_back(*jet_it);

    }

    if(njets.size() < 2) return;
    ++n_after2Jets;
    //After jet selection
    if(njets.size() > 0)
    {
        if(isDiMuon)
        {


            h_AJS_ptLepMuMu->Fill(posMu.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(posMu.pt(),theWeight);
            h_AJS_ptLepMuMu->Fill(negMu.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(negMu.pt(),theWeight);

            h_AJS_etaLepMuMu->Fill(posMu.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(posMu.eta(),theWeight);
            h_AJS_etaLepMuMu->Fill(negMu.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(negMu.eta(),theWeight);
        }
        if(isDiElectron)
        {

            h_AJS_ptLepElEl->Fill(posEl.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(posEl.pt(),theWeight);
            h_AJS_ptLepElEl->Fill(negEl.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(negEl.pt(),theWeight);

            h_AJS_etaLepElEl->Fill(posEl.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(posEl.eta(),theWeight);
            h_AJS_etaLepElEl->Fill(negEl.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(negEl.eta(),theWeight);


        }
        if(isElMu)
        {

            h_AJS_ptLepElMu->Fill(posEl.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(posEl.pt(),theWeight);
            h_AJS_ptLepElMu->Fill(negMu.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(negMu.pt(),theWeight);

            h_AJS_etaLepElMu->Fill(posEl.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(posEl.eta(),theWeight);
            h_AJS_etaLepElMu->Fill(negMu.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(negMu.eta(),theWeight);

        }
        if(isMuEl)
        {

            h_AJS_ptLepElMu->Fill(posMu.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(posMu.pt(),theWeight);
            h_AJS_ptLepElMu->Fill(negEl.pt(),theWeight);
            h_AJS_ptLepDiLep->Fill(negEl.pt(),theWeight);

            h_AJS_etaLepElMu->Fill(posMu.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(posMu.eta(),theWeight);
            h_AJS_etaLepElMu->Fill(negEl.eta(),theWeight);
            h_AJS_etaLepDiLep->Fill(negEl.eta(),theWeight);


        }
    }


    ////////////////////METS PF1//////////////////////////
    const pat::MET &met = mets->front();
    if(isDiMuon && met.pt() < 30) return;
    if(isDiElectron && met.pt() < 30) return;
    if((isElMu || isMuEl) && met.pt() <0) return;
    ++n_afterMet;
    if(met.pt() > 0)
    {
        if(isDiMuon)
        {

            DiLepMass = posMu.mt() + negMu.mt();
            h_AMS_mLepMuMu->Fill(DiLep.M(),theWeight);
            h_AMS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_AMS_ptLepMuMu->Fill(posMu.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(posMu.pt(),theWeight);
            h_AMS_ptLepMuMu->Fill(negMu.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(negMu.pt(),theWeight);
            h_AMS_METMuMu->Fill(met.pt(),theWeight);
            h_AMS_METDiLep->Fill(met.pt(),theWeight);
        }
        if(isDiElectron)
        {
            DiLepMass = posEl.mt() + negEl.mt();
            h_AMS_mLepElEl->Fill(DiLep.M(),theWeight);
            h_AMS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_AMS_ptLepElEl->Fill(posEl.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(posEl.pt(),theWeight);
            h_AMS_ptLepElEl->Fill(negEl.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(negEl.pt(),theWeight);
            h_AMS_METElEl->Fill(met.pt(),theWeight);
            h_AMS_METDiLep->Fill(met.pt(),theWeight);


        }
        if(isElMu)
        {
            DiLepMass = posEl.mt() + negMu.mt();
            h_AMS_mLepElMu->Fill(DiLep.M(),theWeight);
            h_AMS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_AMS_ptLepElMu->Fill(posEl.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(posEl.pt(),theWeight);
            h_AMS_ptLepElMu->Fill(negMu.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(negMu.pt(),theWeight);
            h_AMS_METElMu->Fill(met.pt(),theWeight);
            h_AMS_METDiLep->Fill(met.pt(),theWeight);

        }
        if(isMuEl)
        {
            DiLepMass = posMu.mt() + negEl.mt();
            h_AMS_mLepElMu->Fill(DiLep.M(),theWeight);
            h_AMS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_AMS_ptLepElMu->Fill(posMu.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(posMu.pt(),theWeight);
            h_AMS_ptLepElMu->Fill(negEl.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(negEl.pt(),theWeight);
            h_AMS_METElMu->Fill(met.pt(),theWeight);
            h_AMS_METDiLep->Fill(met.pt(),theWeight);


        }
    }

    if(bjets.size() < 2) return;
    ++n_after2BJets;
    if(bjets.size() > 0)
    {
        if(isDiMuon)
        {
            h_ABS_NBJetsMuMu->Fill(bjets.size());
            h_ABS_NJetsMuMu->Fill(njets.size());
            h_ABS_NBJetsDiLep->Fill(bjets.size());
            h_ABS_NJetsDiLep->Fill(njets.size());
            h_ABS_etaLeadingJetMuMu->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetMuMu->Fill(bjets.at(0).pt(),theWeight);
            h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);

            DiLepMass = posMu.mt() + negMu.mt();
            h_ABS_mLepMuMu->Fill(DiLep.M(),theWeight);
            h_ABS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_ABS_ptLepMuMu->Fill(posMu.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(posMu.pt(),theWeight);
            h_ABS_ptLepMuMu->Fill(negMu.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(negMu.pt(),theWeight);
            h_ABS_METMuMu->Fill(met.pt(),theWeight);
            h_ABS_METDiLep->Fill(met.pt(),theWeight);
            h_ABS_etaLepMuMu->Fill(posMu.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(posMu.eta(),theWeight);
            h_ABS_etaLepMuMu->Fill(negMu.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(negMu.eta(),theWeight);
        }
        if(isDiElectron)
        {
            h_ABS_NBJetsElEl->Fill(bjets.size());
            h_ABS_NJetsElEl->Fill(njets.size());
            h_ABS_NBJetsDiLep->Fill(bjets.size());
            h_ABS_NJetsDiLep->Fill(njets.size());
            h_ABS_etaLeadingJetElEl->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetElEl->Fill(bjets.at(0).pt(),theWeight);
            h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);
            DiLepMass = posEl.mt() + negEl.mt();
            h_ABS_mLepElEl->Fill(DiLep.M(),theWeight);
            h_ABS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_ABS_ptLepElEl->Fill(posEl.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(posEl.pt(),theWeight);
            h_ABS_ptLepElEl->Fill(negEl.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(negEl.pt(),theWeight);
            h_ABS_METElEl->Fill(met.pt(),theWeight);
            h_ABS_METDiLep->Fill(met.pt(),theWeight);
            h_ABS_etaLepElEl->Fill(posEl.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(posEl.eta(),theWeight);
            h_ABS_etaLepElEl->Fill(negEl.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(negEl.eta(),theWeight);


        }
        if(isElMu)
        {
            h_ABS_NBJetsElMu->Fill(bjets.size());
            h_ABS_NJetsElMu->Fill(njets.size());
            h_ABS_NBJetsDiLep->Fill(bjets.size());
            h_ABS_NJetsDiLep->Fill(njets.size());
            h_ABS_etaLeadingJetElMu->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetElMu->Fill(bjets.at(0).pt(),theWeight);
            h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);
            DiLepMass = posEl.mt() + negMu.mt();
            h_ABS_mLepElMu->Fill(DiLep.M(),theWeight);
            h_ABS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_ABS_ptLepElMu->Fill(posEl.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(posEl.pt(),theWeight);
            h_ABS_ptLepElMu->Fill(negMu.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(negMu.pt(),theWeight);
            h_ABS_METElMu->Fill(met.pt(),theWeight);
            h_ABS_METDiLep->Fill(met.pt(),theWeight);
            h_ABS_etaLepElMu->Fill(posEl.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(posEl.eta(),theWeight);
            h_ABS_etaLepElMu->Fill(negMu.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(negMu.eta(),theWeight);

        }
        if(isMuEl)
        {
            h_ABS_NBJetsElMu->Fill(bjets.size());
            h_ABS_NJetsElMu->Fill(njets.size());
            h_ABS_NBJetsDiLep->Fill(bjets.size());
            h_ABS_NJetsDiLep->Fill(njets.size());
            h_ABS_etaLeadingJetElMu->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetElMu->Fill(bjets.at(0).pt(),theWeight);
            h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
            h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);
            DiLepMass = posMu.mass() + negEl.mass();
            h_ABS_mLepElMu->Fill(DiLep.M(),theWeight);
            h_ABS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_ABS_ptLepElMu->Fill(posMu.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(posMu.pt(),theWeight);
            h_ABS_ptLepElMu->Fill(negEl.pt(),theWeight);
            h_ABS_ptLepDiLep->Fill(negEl.pt(),theWeight);
            h_ABS_METElMu->Fill(met.pt(),theWeight);
            h_ABS_METDiLep->Fill(met.pt(),theWeight);
            h_ABS_etaLepElMu->Fill(posMu.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(posMu.eta(),theWeight);
            h_ABS_etaLepElMu->Fill(negEl.eta(),theWeight);
            h_ABS_etaLepDiLep->Fill(negEl.eta(),theWeight);


        }
    }

    //////////////////////GEN LEVEL////////////////
    if(!isData){
        TLorentzVector W1;
        TLorentzVector W2;
        TLorentzVector t1;
        TLorentzVector t2;
        TLorentzVector ttbar;
        TLorentzVector genPosLep;
        TLorentzVector genNegLep;
        const pat::PackedGenParticle matchedPosMu = mse::getMatchedGenParticle(posMu, genColl,13);

        const reco::GenParticle posWMuMu = mse::getMotherPacked(matchedPosMu);
        const reco::GenParticle topMuMu = mse::getMother(posWMuMu);
        if(matchedPosMu.pt() != 0 )  cout << matchedPosMu.pt() << "muon match mother"<< posWMuMu.pdgId()<<"topMuMu Mother pdgId"<< topMuMu.pdgId() << endl;
        const pat::PackedGenParticle matchedNegMu = mse::getMatchedGenParticle(negMu, genColl,13);
        const reco::GenParticle negWMuMu = mse::getMotherPacked(matchedNegMu);
        const reco::GenParticle antitopMuMu = mse::getMother(negWMuMu);
        if(matchedPosMu.pt()> 0 && matchedNegMu.pt() > 0 && posWMuMu.pdgId() == 24 && topMuMu.pdgId() == 6 && negWMuMu.pdgId() == -24 && antitopMuMu.pdgId() == -6  )
        {
            genPosLep.Clear();
            genNegLep.Clear();
            W1.Clear();
            W2.Clear();
            t1.Clear();
            t2.Clear();
            ttbar.Clear();
            genPosLep.SetPtEtaPhiM(matchedPosMu.pt(),matchedPosMu.eta(),matchedPosMu.phi(),matchedPosMu.mass());
            genNegLep.SetPtEtaPhiM(matchedNegMu.pt(),matchedNegMu.eta(),matchedNegMu.phi(),matchedNegMu.mass());

            W1.SetPtEtaPhiM(posWMuMu.pt(),posWMuMu.eta(),posWMuMu.phi(),posWMuMu.mass());
            W2.SetPtEtaPhiM(negWMuMu.pt(),negWMuMu.eta(),negWMuMu.phi(),negWMuMu.mass());

            t1.SetPtEtaPhiM(topMuMu.pt(),topMuMu.eta(),topMuMu.phi(),topMuMu.mass());
            t2.SetPtEtaPhiM(antitopMuMu.pt(),antitopMuMu.eta(),antitopMuMu.phi(),antitopMuMu.mass());

            ttbar = t1 + t2;
            h_GenTTbarM->Fill(ttbar.M(),theWeight);
            TruthMTT.push_back(ttbar.M());
            TruthMTTMuMu.push_back(ttbar.M());
            //            cout << W1.T() << " t component " << t1.T() << endl;
            genPosLep.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            //            cout << W1.T() << "after boost"<< endl;

            float theta1 = (W1.Angle(genPosLep.Vect()));
            h_cosGenMuMu->Fill(TMath::Cos(theta1),theWeight);
            h_cosGen->Fill(TMath::Cos(theta1),theWeight);
            TruthCos.push_back(TMath::Cos(theta1));
            //            TruthCos.push_back(theta1);
            TruthCosMuMu.push_back(TMath::Cos(theta1));
            genNegLep.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(genNegLep.Vect()));
            h_cosGenMuMu->Fill(TMath::Cos(theta2),theWeight);
            h_cosGen->Fill(TMath::Cos(theta2),theWeight);
            TruthCos.push_back(TMath::Cos(theta2));
            //            TruthCos.push_back(theta2);
            TruthCosMuMu.push_back(TMath::Cos(theta2));




        }

        const pat::PackedGenParticle matchedPosEl = mse::getMatchedGenParticle(posEl, genColl,11);

        const reco::GenParticle posWElEl = mse::getMotherPacked(matchedPosEl);
        const reco::GenParticle topElEl = mse::getMother(posWElEl);
        if(matchedPosEl.pt() != 0 )  cout << matchedPosEl.pt() << "Elon match mother"<< posWElEl.pdgId()<<"topElEl Mother pdgId"<< topElEl.pdgId() << endl;
        const pat::PackedGenParticle matchedNegEl = mse::getMatchedGenParticle(negEl, genColl,11);
        const reco::GenParticle negWElEl = mse::getMotherPacked(matchedNegEl);
        const reco::GenParticle antitopElEl = mse::getMother(negWElEl);
        if(matchedPosEl.pt()> 0 && matchedNegEl.pt() > 0 && posWElEl.pdgId() == 24 && topElEl.pdgId() == 6 && negWElEl.pdgId() == -24 && antitopElEl.pdgId() == -6  )
        {
            genPosLep.Clear();
            genNegLep.Clear();
            W1.Clear();
            W2.Clear();
            t1.Clear();
            t2.Clear();
            ttbar.Clear();
            genPosLep.SetPtEtaPhiM(matchedPosEl.pt(),matchedPosEl.eta(),matchedPosEl.phi(),matchedPosEl.mass());
            genNegLep.SetPtEtaPhiM(matchedNegEl.pt(),matchedNegEl.eta(),matchedNegEl.phi(),matchedNegEl.mass());

            W1.SetPtEtaPhiM(posWElEl.pt(),posWElEl.eta(),posWElEl.phi(),posWElEl.mass());
            W2.SetPtEtaPhiM(negWElEl.pt(),negWElEl.eta(),negWElEl.phi(),negWElEl.mass());

            t1.SetPtEtaPhiM(topElEl.pt(),topElEl.eta(),topElEl.phi(),topElEl.mass());
            t2.SetPtEtaPhiM(antitopElEl.pt(),antitopElEl.eta(),antitopElEl.phi(),antitopElEl.mass());

            ttbar = t1 + t2;
            h_GenTTbarM->Fill(ttbar.M(),theWeight);
            TruthMTT.push_back(ttbar.M());


            genPosLep.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());

            float theta1 = (W1.Angle(genPosLep.Vect()));
            h_cosGenElEl->Fill(TMath::Cos(theta1),theWeight);
            h_cosGen->Fill(TMath::Cos(theta1),theWeight);
            TruthCos.push_back(TMath::Cos(theta1));
            //            TruthCos.push_back(theta1);
            genNegLep.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(genNegLep.Vect()));
            h_cosGenElEl->Fill(TMath::Cos(theta2),theWeight);
            h_cosGen->Fill(TMath::Cos(theta2),theWeight);
            TruthCos.push_back(TMath::Cos(theta2));
            //            TruthCos.push_back(theta2);

        }

        const pat::PackedGenParticle matchedPosEMEl = mse::getMatchedGenParticle(posEl, genColl,11);

        const reco::GenParticle posWElMu = mse::getMotherPacked(matchedPosEMEl);
        const reco::GenParticle topElMu = mse::getMother(posWElMu);
        if(matchedPosEl.pt() != 0 )  cout << matchedPosEMEl.pt() << "Elon match mother"<< posWElMu.pdgId()<<"topElEl Mother pdgId"<< topElMu.pdgId() << endl;
        const pat::PackedGenParticle matchedNegEMMu = mse::getMatchedGenParticle(negMu, genColl,13);
        const reco::GenParticle negWElMu = mse::getMotherPacked(matchedNegEMMu);
        const reco::GenParticle antitopElMu = mse::getMother(negWElMu);
        if(matchedPosEMEl.pt() > 0 && matchedNegEMMu.pt() > 0 && posWElMu.pdgId() == 24 && topElMu.pdgId() == 6 && negWElMu.pdgId() == -24 && antitopElMu.pdgId() == -6  )
        {
            genPosLep.Clear();
            genNegLep.Clear();
            W1.Clear();
            W2.Clear();
            t1.Clear();
            t2.Clear();
            ttbar.Clear();
            genPosLep.SetPtEtaPhiM(matchedPosEMEl.pt(),matchedPosEMEl.eta(),matchedPosEMEl.phi(),matchedPosEMEl.mass());
            genNegLep.SetPtEtaPhiM(matchedNegEMMu.pt(),matchedNegEMMu.eta(),matchedNegEMMu.phi(),matchedNegEMMu.mass());

            W1.SetPtEtaPhiM(posWElMu.pt(),posWElMu.eta(),posWElMu.phi(),posWElMu.mass());
            W2.SetPtEtaPhiM(negWElMu.pt(),negWElMu.eta(),negWElMu.phi(),negWElMu.mass());

            t1.SetPtEtaPhiM(topElMu.pt(),topElMu.eta(),topElMu.phi(),topElMu.mass());
            t2.SetPtEtaPhiM(antitopElMu.pt(),antitopElMu.eta(),antitopElMu.phi(),antitopElMu.mass());

            ttbar = t1 + t2;
            h_GenTTbarM->Fill(ttbar.M(),theWeight);
            TruthMTT.push_back(ttbar.M());


            genPosLep.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(genPosLep.Vect()));
            h_cosGenElMu->Fill(TMath::Cos(theta1),theWeight);
            h_cosGen->Fill(TMath::Cos(theta1),theWeight);
            TruthCos.push_back(TMath::Cos(theta1));
            //            TruthCos.push_back(theta1);
            genNegLep.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(genNegLep.Vect()));
            h_cosGenElMu->Fill(TMath::Cos(theta2),theWeight);
            h_cosGen->Fill(TMath::Cos(theta2),theWeight);
            TruthCos.push_back(TMath::Cos(theta2));
            //            TruthCos.push_back(theta2);
        }

        const pat::PackedGenParticle matchedPosEMMu = mse::getMatchedGenParticle(posMu, genColl,13);

        const reco::GenParticle posWElMu2 = mse::getMotherPacked(matchedPosEMMu);
        const reco::GenParticle topElMu2 = mse::getMother(posWElMu2);
        if(matchedPosEl.pt() != 0 )  cout << matchedPosEMMu.pt() << "Elon match mother"<< posWElMu2.pdgId()<<"topElEl Mother pdgId"<< topElMu2.pdgId() << endl;
        const pat::PackedGenParticle matchedNegEMEl = mse::getMatchedGenParticle(negEl, genColl,11);
        const reco::GenParticle negWElMu2 = mse::getMotherPacked(matchedNegEMEl);
        const reco::GenParticle antitopElMu2 = mse::getMother(negWElMu2);
        if(matchedPosEMMu.pt()> 0 && matchedNegEMEl.pt() > 0 && posWElMu2.pdgId() == 24 && topElMu2.pdgId() == 6 && negWElMu2.pdgId() == -24 && antitopElMu2.pdgId() == -6  )
        {
            genPosLep.Clear();
            genNegLep.Clear();
            W1.Clear();
            W2.Clear();
            t1.Clear();
            t2.Clear();
            ttbar.Clear();
            genPosLep.SetPtEtaPhiM(matchedPosEMMu.pt(),matchedPosEMMu.eta(),matchedPosEMMu.phi(),matchedPosEMMu.mass());
            genNegLep.SetPtEtaPhiM(matchedNegEMEl.pt(),matchedNegEMEl.eta(),matchedNegEMEl.phi(),matchedNegEMEl.mass());

            W1.SetPtEtaPhiM(posWElMu2.pt(),posWElMu2.eta(),posWElMu2.phi(),posWElMu2.mass());
            W2.SetPtEtaPhiM(negWElMu2.pt(),negWElMu2.eta(),negWElMu2.phi(),negWElMu2.mass());

            t1.SetPtEtaPhiM(topElMu2.pt(),topElMu2.eta(),topElMu2.phi(),topElMu2.mass());
            t2.SetPtEtaPhiM(antitopElMu2.pt(),antitopElMu2.eta(),antitopElMu2.phi(),antitopElMu2.mass());

            ttbar = t1 + t2;
            h_GenTTbarM->Fill(ttbar.M(),theWeight);
            TruthMTT.push_back(ttbar.M());


            genPosLep.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(genPosLep.Vect()));
            h_cosGenElMu->Fill(TMath::Cos(theta1),theWeight);
            h_cosGen->Fill(TMath::Cos(theta1),theWeight);
            TruthCos.push_back(TMath::Cos(theta1));
            //            TruthCos.push_back(theta1);
            genNegLep.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(genNegLep.Vect()));
            h_cosGenElMu->Fill(TMath::Cos(theta2),theWeight);
            h_cosGen->Fill(TMath::Cos(theta2),theWeight);
            TruthCos.push_back(TMath::Cos(theta2));
            //            TruthCos.push_back(theta2);
        }



    }


    //////////////////////TOP RECO/////////////////
    //    TtDilepEvtSolution asol;

    //    asol.setGenEvt(genEvent);
    ////////DiMuon:
    if(isDiMuon && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negMu.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posMu.p4()) < 0.5) ){
        // cout << "IM HERE!MuMu" << endl;
        TLorentzVector W1;
        TLorentzVector W2;
        TLorentzVector t1;
        TLorentzVector t2;
        TLorentzVector ttbar;
        TLorentzVector lepPos;
        TLorentzVector lepNeg;
        TLorentzVector BJet;
        TLorentzVector BBJet;
        TLorentzVector dilepton;
        lepPos.SetPtEtaPhiE(posMu.pt(),posMu.eta(),posMu.phi(),posMu.energy());
        lepNeg.SetPtEtaPhiE(negMu.pt(),negMu.eta(),negMu.phi(),negMu.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());
        //        cout << posMu.genParticle()->mother()->pdgId() << endl;

                amwtSolver->SetConstraints(met.px(),met.py());

                TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);

//        amwtsolver->SetConstraints(met.px(),met.py());
//        TtFullLepKinSolver::NeutrinoSolution nuSol= amwtsolver->getNuSolution(lepPos,lepNeg,BJet,BBJet);
        // cout << nuSol.neutrino.p4() << " neutrino "<< endl;

        dilepton = lepPos + lepNeg;

        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 && dilepton.M() > 76)){
            // cout << "what 0" << endl;
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;
            // cout << "what 1" << endl;
            h_yWMuMu->Fill(W1.Rapidity(),theWeight);
            h_yWMuMu->Fill(W2.Rapidity(),theWeight);
            h_yWDiLep->Fill(W1.Rapidity(),theWeight);
            h_yWDiLep->Fill(W2.Rapidity(),theWeight);
            h_ptWMuMu->Fill(W1.Pt(),theWeight);
            h_ptWMuMu->Fill(W2.Pt(),theWeight);
            h_ptWDiLep->Fill(W1.Pt(),theWeight);
            h_ptWDiLep->Fill(W2.Pt(),theWeight);

            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            // cout << "what 2" << endl;
            h_yTMuMu->Fill(t1.Rapidity(),theWeight);
            h_yTMuMu->Fill(t2.Rapidity(),theWeight);
            h_yTDiLep->Fill(t1.Rapidity(),theWeight);
            h_yTDiLep->Fill(t2.Rapidity(),theWeight);
            h_ptTMuMu->Fill(t1.Pt(),theWeight);
            h_ptTMuMu->Fill(t2.Pt(),theWeight);
            h_ptTDiLep->Fill(t1.Pt(),theWeight);
            h_ptTDiLep->Fill(t2.Pt(),theWeight);
            h_mTTbarMuMu->Fill(ttbar.M(),theWeight);
            h_TTbarM->Fill(ttbar.M(),theWeight);

            ++n_afterTop;
            RecoMTT.push_back(ttbar.M());
            RecoMTTMuMu.push_back(ttbar.M());


            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            RecoCos.push_back(TMath::Cos(theta1));
            //            RecoCos.push_back(theta1);
            RecoCosMuMu.push_back(TMath::Cos(theta1));
            h_cosMuMu->Fill(TMath::Cos(theta1),theWeight);
            h_cosDiLep->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            RecoCos.push_back(TMath::Cos(theta2));
            //            RecoCos.push_back(theta2);
            RecoCosMuMu.push_back(TMath::Cos(theta2));
            h_cosMuMu->Fill(TMath::Cos(theta2),theWeight);
            h_cosDiLep->Fill(TMath::Cos(theta2),theWeight);

            if(!isData){
                TLorentzVector genW1;
                TLorentzVector genW2;
                TLorentzVector gent1;
                TLorentzVector gent2;
                TLorentzVector genttbar;
                TLorentzVector genPosLep;
                TLorentzVector genNegLep;
                const pat::PackedGenParticle matchedPosMu = mse::getMatchedGenParticle(posMu, genColl,13);

                const reco::GenParticle posWMuMu = mse::getMotherPacked(matchedPosMu);
                const reco::GenParticle topMuMu = mse::getMother(posWMuMu);
                if(matchedPosMu.pt() != 0 )  cout << matchedPosMu.pt() << "muon match mother"<< posWMuMu.pdgId()<<"topMuMu Mother pdgId"<< topMuMu.pdgId() << endl;
                const pat::PackedGenParticle matchedNegMu = mse::getMatchedGenParticle(negMu, genColl,13);
                const reco::GenParticle negWMuMu = mse::getMotherPacked(matchedNegMu);
                const reco::GenParticle antitopMuMu = mse::getMother(negWMuMu);
                if(matchedPosMu.pt()> 0 && matchedNegMu.pt() > 0 && posWMuMu.pdgId() == 24 && topMuMu.pdgId() == 6 && negWMuMu.pdgId() == -24 && antitopMuMu.pdgId() == -6  )
                {
                    genPosLep.Clear();
                    genNegLep.Clear();
                    genW1.Clear();
                    genW2.Clear();
                    gent1.Clear();
                    gent2.Clear();
                    genttbar.Clear();
                    genPosLep.SetPtEtaPhiM(matchedPosMu.pt(),matchedPosMu.eta(),matchedPosMu.phi(),matchedPosMu.mass());
                    genNegLep.SetPtEtaPhiM(matchedNegMu.pt(),matchedNegMu.eta(),matchedNegMu.phi(),matchedNegMu.mass());

                    genW1.SetPtEtaPhiM(posWMuMu.pt(),posWMuMu.eta(),posWMuMu.phi(),posWMuMu.mass());
                    genW2.SetPtEtaPhiM(negWMuMu.pt(),negWMuMu.eta(),negWMuMu.phi(),negWMuMu.mass());

                    gent1.SetPtEtaPhiM(topMuMu.pt(),topMuMu.eta(),topMuMu.phi(),topMuMu.mass());
                    gent2.SetPtEtaPhiM(antitopMuMu.pt(),antitopMuMu.eta(),antitopMuMu.phi(),antitopMuMu.mass());

                    genttbar = gent1 + gent2;
                    genPosLep.Boost(-genW1.BoostVector());
                    genW1.Boost(-gent1.BoostVector());


                    float gentheta1 = (genW1.Angle(genPosLep.Vect()));


                    genNegLep.Boost(-genW2.BoostVector());
                    genW2.Boost(-gent2.BoostVector());
                    float gentheta2 = (genW2.Angle(genNegLep.Vect()));
                    h_truthRecoMuMu->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoCos->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoMTT->Fill(ttbar.M(),genttbar.M());
                    h_truthRecoMuMu->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_truthRecoCos->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_truthMinusRecoMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }
        h_NBJetsMuMu->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);


    }


    //////DiElectron:
    if(isDiElectron && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negEl.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posEl.p4()) < 0.5) ){
        // cout << "IM HERE! ElEl" << endl;
        TLorentzVector W1;
        TLorentzVector W2;
        TLorentzVector t1;
        TLorentzVector t2;
        TLorentzVector ttbar;
        TLorentzVector lepPos;
        TLorentzVector lepNeg;
        TLorentzVector BJet;
        TLorentzVector BBJet;
        TLorentzVector dilepton;
        lepPos.Clear();
        lepNeg.Clear();

        lepPos.SetPtEtaPhiE(posEl.pt(),posEl.eta(),posEl.phi(),posEl.energy());
        lepNeg.SetPtEtaPhiE(negEl.pt(),negEl.eta(),negEl.phi(),negEl.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());

        //        amwtsolver->SetConstraints(met.px(),met.py());

        //        TtFullLepKinSolver::NeutrinoSolution nuSol=amwtsolver->getNuSolution(lepPos,lepNeg,BJet,BBJet);

        amwtSolver->SetConstraints(met.px(),met.py());
        TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);
        // cout << nuSol.neutrino.p4() << " neutrino "<< endl;
        // cout << nuSol.neutrinoBar.p4() << " neutrinoBar "<< endl;
        dilepton = lepPos + lepNeg;
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 && dilepton.M() > 76)){
            // cout << "what 0" << endl;
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.Clear();
            nuBar.Clear();
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;
            // cout << "what 1" << endl;
            h_yWElEl->Fill(W1.Rapidity(),theWeight);
            h_yWElEl->Fill(W2.Rapidity(),theWeight);
            h_yWDiLep->Fill(W1.Rapidity(),theWeight);
            h_yWDiLep->Fill(W2.Rapidity(),theWeight);
            h_ptWElEl->Fill(W1.Pt(),theWeight);
            h_ptWElEl->Fill(W2.Pt(),theWeight);
            h_ptWDiLep->Fill(W1.Pt(),theWeight);
            h_ptWDiLep->Fill(W2.Pt(),theWeight);

            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            // cout << "what 2" << endl;
            h_yTElEl->Fill(t1.Rapidity(),theWeight);
            h_yTElEl->Fill(t2.Rapidity(),theWeight);
            h_yTDiLep->Fill(t1.Rapidity(),theWeight);
            h_yTDiLep->Fill(t2.Rapidity(),theWeight);
            h_ptTElEl->Fill(t1.Pt(),theWeight);
            h_ptTElEl->Fill(t2.Pt(),theWeight);
            h_ptTDiLep->Fill(t1.Pt(),theWeight);
            h_ptTDiLep->Fill(t2.Pt(),theWeight);
            h_mTTbarElEl->Fill(ttbar.M(),theWeight);
            h_TTbarM->Fill(ttbar.M(),theWeight);
            ++n_afterTop;
            RecoMTT.push_back(ttbar.M());





            // cout << "what 3" << endl;
            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            RecoCos.push_back(TMath::Cos(theta1));
            //            RecoCos.push_back(theta1);
            // cout << "what 3.1" << endl;
            h_cosElEl->Fill(TMath::Cos(theta1),theWeight);
            // cout << "what 3.2" << endl;
            h_cosDiLep->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            RecoCos.push_back(TMath::Cos(theta2));
            //            RecoCos.push_back(theta2);
            h_cosElEl->Fill(TMath::Cos(theta2),theWeight);
            h_cosDiLep->Fill(TMath::Cos(theta2),theWeight);
            // cout << "what 4" << endl;
            if(!isData){
                TLorentzVector genW1;
                TLorentzVector genW2;
                TLorentzVector gent1;
                TLorentzVector gent2;
                TLorentzVector genttbar;
                TLorentzVector genPosLep;
                TLorentzVector genNegLep;
                const pat::PackedGenParticle matchedPosEl = mse::getMatchedGenParticle(posEl, genColl,11);

                const reco::GenParticle posWElEl = mse::getMotherPacked(matchedPosEl);
                const reco::GenParticle topElEl = mse::getMother(posWElEl);
                const pat::PackedGenParticle matchedNegEl = mse::getMatchedGenParticle(negEl, genColl,11);
                const reco::GenParticle negWElEl = mse::getMotherPacked(matchedNegEl);
                const reco::GenParticle antitopElEl = mse::getMother(negWElEl);
                if(matchedPosEl.pt()> 0 && matchedNegEl.pt() > 0 && posWElEl.pdgId() == 24 && topElEl.pdgId() == 6 && negWElEl.pdgId() == -24 && antitopElEl.pdgId() == -6  )
                {
                    genPosLep.Clear();
                    genNegLep.Clear();
                    genW1.Clear();
                    genW2.Clear();
                    gent1.Clear();
                    gent2.Clear();
                    genttbar.Clear();
                    genPosLep.SetPtEtaPhiM(matchedPosEl.pt(),matchedPosEl.eta(),matchedPosEl.phi(),matchedPosEl.mass());
                    genNegLep.SetPtEtaPhiM(matchedNegEl.pt(),matchedNegEl.eta(),matchedNegEl.phi(),matchedNegEl.mass());


                    genW1.SetPtEtaPhiM(posWElEl.pt(),posWElEl.eta(),posWElEl.phi(),posWElEl.mass());
                    genW2.SetPtEtaPhiM(negWElEl.pt(),negWElEl.eta(),negWElEl.phi(),negWElEl.mass());

                    gent1.SetPtEtaPhiM(topElEl.pt(),topElEl.eta(),topElEl.phi(),topElEl.mass());
                    gent2.SetPtEtaPhiM(antitopElEl.pt(),antitopElEl.eta(),antitopElEl.phi(),antitopElEl.mass());

                    genttbar = gent1 + gent2;
                    genPosLep.Boost(-genW1.BoostVector());
                    genW1.Boost(-gent1.BoostVector());


                    float gentheta1 = (genW1.Angle(genPosLep.Vect()));

                    genNegLep.Boost(-genW2.BoostVector());
                    genW2.Boost(-gent2.BoostVector());
                    float gentheta2 = (genW2.Angle(genNegLep.Vect()));

                    h_truthRecoCos->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoMTT->Fill(ttbar.M(),genttbar.M());
                    h_truthRecoCos->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_truthMinusRecoMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }
        // cout << "what 5" << endl;
        h_NBJetsElEl->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);
        // cout << "what 6" << endl;
    }
    //////ElectronMuon:
    if(isElMu && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negMu.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posEl.p4()) < 0.5) ){
        // cout << "IM HERE!ElMu" << endl;
        TLorentzVector W1;
        TLorentzVector W2;
        TLorentzVector t1;
        TLorentzVector t2;
        TLorentzVector ttbar;
        TLorentzVector lepPos;
        TLorentzVector lepNeg;
        TLorentzVector BJet;
        TLorentzVector BBJet;
        TLorentzVector dilepton;
        lepPos.SetPtEtaPhiE(posEl.pt(),posEl.eta(),posEl.phi(),posEl.energy());
        lepNeg.SetPtEtaPhiE(negMu.pt(),negMu.eta(),negMu.phi(),negMu.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());


        //        amwtsolver->SetConstraints(met.px(),met.py());

        //        TtFullLepKinSolver::NeutrinoSolution nuSol=amwtsolver->getNuSolution(lepPos,lepNeg,BJet,BBJet);
        amwtSolver->SetConstraints(met.px(),met.py());
        TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);
        // cout << nuSol.neutrino.p4() << " neutrino "<< endl;
        dilepton = lepPos + lepNeg;
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 && dilepton.M() > 76)){
            // cout << "what 0" << endl;
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;
            // cout << "what 1" << endl;
            h_yWElMu->Fill(W1.Rapidity(),theWeight);
            h_yWElMu->Fill(W2.Rapidity(),theWeight);
            h_yWDiLep->Fill(W1.Rapidity(),theWeight);
            h_yWDiLep->Fill(W2.Rapidity(),theWeight);
            h_ptWElMu->Fill(W1.Pt(),theWeight);
            h_ptWElMu->Fill(W2.Pt(),theWeight);
            h_ptWDiLep->Fill(W1.Pt(),theWeight);
            h_ptWDiLep->Fill(W2.Pt(),theWeight);

            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            // cout << "what 2" << endl;
            h_yTElMu->Fill(t1.Rapidity(),theWeight);
            h_yTElMu->Fill(t2.Rapidity(),theWeight);
            h_yTDiLep->Fill(t1.Rapidity(),theWeight);
            h_yTDiLep->Fill(t2.Rapidity(),theWeight);
            h_ptTElMu->Fill(t1.Pt(),theWeight);
            h_ptTElMu->Fill(t2.Pt(),theWeight);
            h_ptTDiLep->Fill(t1.Pt(),theWeight);
            h_ptTDiLep->Fill(t2.Pt(),theWeight);
            h_mTTbarElMu->Fill(ttbar.M(),theWeight);
            h_TTbarM->Fill(ttbar.M(),theWeight);

            ++n_afterTop;
            RecoMTT.push_back(ttbar.M());



            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            RecoCos.push_back(TMath::Cos(theta1));
            //            RecoCos.push_back(theta1);
            h_cosElMu->Fill(TMath::Cos(theta1),theWeight);
            h_cosDiLep->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            RecoCos.push_back(TMath::Cos(theta2));
            //            RecoCos.push_back(theta2);
            h_cosElMu->Fill(TMath::Cos(theta2),theWeight);
            h_cosDiLep->Fill(TMath::Cos(theta2),theWeight);

            if(!isData){
                TLorentzVector genW1;
                TLorentzVector genW2;
                TLorentzVector gent1;
                TLorentzVector gent2;
                TLorentzVector genttbar;
                TLorentzVector genPosLep;
                TLorentzVector genNegLep;
                const pat::PackedGenParticle matchedPosEMEl = mse::getMatchedGenParticle(posEl, genColl,11);

                const reco::GenParticle posWElMu = mse::getMotherPacked(matchedPosEMEl);
                const reco::GenParticle topElMu = mse::getMother(posWElMu);

                const pat::PackedGenParticle matchedNegEMMu = mse::getMatchedGenParticle(negMu, genColl,13);
                const reco::GenParticle negWElMu = mse::getMotherPacked(matchedNegEMMu);
                const reco::GenParticle antitopElMu = mse::getMother(negWElMu);
                if(matchedPosEMEl.pt() > 0 && matchedNegEMMu.pt() > 0 && posWElMu.pdgId() == 24 && topElMu.pdgId() == 6 && negWElMu.pdgId() == -24 && antitopElMu.pdgId() == -6  )
                {
                    genPosLep.Clear();
                    genNegLep.Clear();
                    genW1.Clear();
                    genW2.Clear();
                    gent1.Clear();
                    gent2.Clear();
                    genttbar.Clear();
                    genPosLep.SetPtEtaPhiM(matchedPosEMEl.pt(),matchedPosEMEl.eta(),matchedPosEMEl.phi(),matchedPosEMEl.mass());
                    genNegLep.SetPtEtaPhiM(matchedNegEMMu.pt(),matchedNegEMMu.eta(),matchedNegEMMu.phi(),matchedNegEMMu.mass());


                    genW1.SetPtEtaPhiM(posWElMu.pt(),posWElMu.eta(),posWElMu.phi(),posWElMu.mass());
                    genW2.SetPtEtaPhiM(negWElMu.pt(),negWElMu.eta(),negWElMu.phi(),negWElMu.mass());

                    gent1.SetPtEtaPhiM(topElMu.pt(),topElMu.eta(),topElMu.phi(),topElMu.mass());
                    gent2.SetPtEtaPhiM(antitopElMu.pt(),antitopElMu.eta(),antitopElMu.phi(),antitopElMu.mass());


                    genttbar = gent1 + gent2;
                    genPosLep.Boost(-genW1.BoostVector());
                    genW1.Boost(-gent1.BoostVector());


                    float gentheta1 = (genW1.Angle(genPosLep.Vect()));

                    genNegLep.Boost(-genW2.BoostVector());
                    genW2.Boost(-gent2.BoostVector());
                    float gentheta2 = (genW2.Angle(genNegLep.Vect()));
                    h_truthRecoCos->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoMTT->Fill(ttbar.M(),genttbar.M());
                    h_truthRecoCos->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_truthMinusRecoMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }
        h_NBJetsElMu->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);
    }

    if(isMuEl && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negEl.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posMu.p4()) < 0.5) ){
        // cout << "IM HERE! MuEl" << endl;
        TLorentzVector W1;
        TLorentzVector W2;
        TLorentzVector t1;
        TLorentzVector t2;
        TLorentzVector ttbar;
        TLorentzVector lepPos;
        TLorentzVector lepNeg;
        TLorentzVector BJet;
        TLorentzVector BBJet;
        TLorentzVector dilepton;
        lepPos.SetPtEtaPhiE(posMu.pt(),posMu.eta(),posMu.phi(),posMu.energy());
        lepNeg.SetPtEtaPhiE(negEl.pt(),negEl.eta(),negEl.phi(),negEl.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());


        //        amwtsolver->SetConstraints(met.px(),met.py());

        //        TtFullLepKinSolver::NeutrinoSolution nuSol=amwtsolver->getNuSolution(lepPos,lepNeg,BJet,BBJet);
        amwtSolver->SetConstraints(met.px(),met.py());
        TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);
        // cout << nuSol.neutrino.p4() << " neutrino "<< endl;
        dilepton = lepPos + lepNeg;
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 && dilepton.M() > 76)){
            // cout << "what 0" << endl;
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;
            // cout << "what 1" << endl;
            h_yWElMu->Fill(W1.Rapidity(),theWeight);
            h_yWElMu->Fill(W2.Rapidity(),theWeight);
            h_yWDiLep->Fill(W1.Rapidity(),theWeight);
            h_yWDiLep->Fill(W2.Rapidity(),theWeight);
            h_ptWElMu->Fill(W1.Pt(),theWeight);
            h_ptWElMu->Fill(W2.Pt(),theWeight);
            h_ptWDiLep->Fill(W1.Pt(),theWeight);
            h_ptWDiLep->Fill(W2.Pt(),theWeight);

            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            // cout << "what 2" << endl;
            h_yTElMu->Fill(t1.Rapidity(),theWeight);
            h_yTElMu->Fill(t2.Rapidity(),theWeight);
            h_yTDiLep->Fill(t1.Rapidity(),theWeight);
            h_yTDiLep->Fill(t2.Rapidity(),theWeight);
            h_ptTElMu->Fill(t1.Pt(),theWeight);
            h_ptTElMu->Fill(t2.Pt(),theWeight);
            h_ptTDiLep->Fill(t1.Pt(),theWeight);
            h_ptTDiLep->Fill(t2.Pt(),theWeight);
            h_mTTbarElMu->Fill(ttbar.M(),theWeight);
            h_TTbarM->Fill(ttbar.M(),theWeight);

            ++n_afterTop;
            RecoMTT.push_back(ttbar.M());



            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            RecoCos.push_back(TMath::Cos(theta1));
            //            RecoCos.push_back(theta1);
            h_cosElMu->Fill(TMath::Cos(theta1),theWeight);
            h_cosDiLep->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            RecoCos.push_back(TMath::Cos(theta2));
            //            RecoCos.push_back(theta1);
            h_cosElMu->Fill(TMath::Cos(theta2),theWeight);
            h_cosDiLep->Fill(TMath::Cos(theta2),theWeight);

            if(!isData){
                TLorentzVector genW1;
                TLorentzVector genW2;
                TLorentzVector gent1;
                TLorentzVector gent2;
                TLorentzVector genttbar;
                TLorentzVector genPosLep;
                TLorentzVector genNegLep;
                const pat::PackedGenParticle matchedPosEMMu = mse::getMatchedGenParticle(posMu, genColl,13);

                const reco::GenParticle posWElMu2 = mse::getMotherPacked(matchedPosEMMu);
                const reco::GenParticle topElMu2 = mse::getMother(posWElMu2);
                const pat::PackedGenParticle matchedNegEMEl = mse::getMatchedGenParticle(negEl, genColl,11);
                const reco::GenParticle negWElMu2 = mse::getMotherPacked(matchedNegEMEl);
                const reco::GenParticle antitopElMu2 = mse::getMother(negWElMu2);
                if(matchedPosEMMu.pt()> 0 && matchedNegEMEl.pt() > 0 && posWElMu2.pdgId() == 24 && topElMu2.pdgId() == 6 && negWElMu2.pdgId() == -24 && antitopElMu2.pdgId() == -6  )
                {
                    genPosLep.Clear();
                    genNegLep.Clear();
                    genW1.Clear();
                    genW2.Clear();
                    gent1.Clear();
                    gent2.Clear();
                    genttbar.Clear();
                    genPosLep.SetPtEtaPhiM(matchedPosEMMu.pt(),matchedPosEMMu.eta(),matchedPosEMMu.phi(),matchedPosEMMu.mass());
                    genNegLep.SetPtEtaPhiM(matchedNegEMEl.pt(),matchedNegEMEl.eta(),matchedNegEMEl.phi(),matchedNegEMEl.mass());


                    genW1.SetPtEtaPhiM(posWElMu2.pt(),posWElMu2.eta(),posWElMu2.phi(),posWElMu2.mass());
                    genW2.SetPtEtaPhiM(negWElMu2.pt(),negWElMu2.eta(),negWElMu2.phi(),negWElMu2.mass());

                    gent1.SetPtEtaPhiM(topElMu2.pt(),topElMu2.eta(),topElMu2.phi(),topElMu2.mass());
                    gent2.SetPtEtaPhiM(antitopElMu2.pt(),antitopElMu2.eta(),antitopElMu2.phi(),antitopElMu2.mass());


                    genttbar = gent1 + gent2;
                    genPosLep.Boost(-genW1.BoostVector());
                    genW1.Boost(-gent1.BoostVector());


                    float gentheta1 = (genW1.Angle(genPosLep.Vect()));

                    genNegLep.Boost(-genW2.BoostVector());
                    genW2.Boost(-gent2.BoostVector());
                    float gentheta2 = (genW2.Angle(genNegLep.Vect()));
                    h_truthRecoCos->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoMTT->Fill(ttbar.M(),genttbar.M());
                    h_truthRecoCos->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_truthMinusRecoCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_truthMinusRecoMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }
        h_NBJetsElMu->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);
    }







    t_outTree->Fill();


}

void MiniAnalyzer::beginJob() {
    TH1::SetDefaultSumw2();

    h_cosMuMu = fs->make<TH1F>("h_cosMuMu",";cos(#theta);",100,-1,1);
    h_cosElEl = fs->make<TH1F>("h_cosElEl",";cos(#theta);",100,-1,1);
    h_cosElMu = fs->make<TH1F>("h_cosElMu",";cos(#theta);",100,-1,1);
    h_cosDiLep = fs->make<TH1F>("h_cosDiLep",";cos(#theta);",100,-1,1);
    h_cosGenElMu = fs->make<TH1F>("h_cosGenElMu",";cos(#theta);",100,-1,1);
    h_cosGenMuMu = fs->make<TH1F>("h_cosGenMuMu",";cos(#theta);",100,-1,1);
    h_cosGenElEl = fs->make<TH1F>("h_cosGenElEl",";cos(#theta);",100,-1,1);
    h_cosGen = fs->make<TH1F>("h_cosGen",";cos(#theta);",100,-1,1);
    h_GenTTbarM = fs->make<TH1F>("h_GenTTbarM",";M_{TTbar};",100,90,1300);
    h_etaMu = fs->make<TH1F>("h_EtaMu",";#eta_{l};",100,-3,3);
    h_TTbarM = fs->make<TH1F>("h_TTbarM",";M_{TTbar};",100,0,1300);
    h_PtMu = fs->make<TH1F>("h_PtMu",";Mu_{Pt};",100,0.,200.);
    h_NPV = fs->make<TH1F>("h_NPV",";NPV;",100,0,30);
    h_NBJetsMuMu = fs->make<TH1F>("h_NBJetsMuMu",";N_{BJets};",30,0,8);
    h_NBJetsElEl = fs->make<TH1F>("h_NBJetsElEl",";N_{BJets};",30,0,8);
    h_NBJetsElMu = fs->make<TH1F>("h_NBJetsElMu",";N_{BJets};",30,0,8);
    h_NBJetsDiLep = fs->make<TH1F>("h_NBJetsDiLep",";N_{BJets};",30,0,8);
    h_yTMuMu = fs->make<TH1F>("h_yTMuMu",";y^{Top};",100,-3,3);
    h_yTElEl = fs->make<TH1F>("h_yTElEl",";y^{Top};",100,-3,3);
    h_yTElMu = fs->make<TH1F>("h_yTElMu",";y^{Top};",100,-3,3);
    h_yTDiLep = fs->make<TH1F>("h_yTDiLep",";y^{Top};",100,-3,3);
    h_yWMuMu = fs->make<TH1F>("h_yWMuMu",";y^{W};",100,-3,3);
    h_yWElEl = fs->make<TH1F>("h_yWElEl",";y^{W};",100,-3,3);
    h_yWElMu = fs->make<TH1F>("h_yWElMu",";y^{W};",100,-3,3);
    h_yWDiLep = fs->make<TH1F>("h_yWDiLep",";y^{W};",100,-3,3);
    h_ptTMuMu= fs->make<TH1F>("h_ptTMuMu",";Pt^{Top};",100,0.,300.);
    h_ptTElEl= fs->make<TH1F>("h_ptTElEl",";Pt^{Top};",100,0.,300.);
    h_ptTElMu= fs->make<TH1F>("h_ptTElMu",";Pt^{Top};",100,0.,300.);
    h_ptTDiLep= fs->make<TH1F>("h_ptTDiLep",";Pt^{Top};",100,0.,300.);
    h_ptWMuMu= fs->make<TH1F>("h_ptWMuMu",";Pt^{W};",100,0.,1000.);
    h_ptWElEl= fs->make<TH1F>("h_ptWElEl",";Pt^{W};",100,0.,1000.);
    h_ptWElMu= fs->make<TH1F>("h_ptWElMu",";Pt^{W};",100,0.,1000.);
    h_ptWDiLep= fs->make<TH1F>("h_ptWDiLep",";Pt^{W};",100,0.,1000.);
    h_mTTbarMuMu =fs->make<TH1F>("h_mTTbarMuMu",";M^{tt};",100,0,1300);
    h_mTTbarElEl =fs->make<TH1F>("h_mTTbarElEl",";M^{tt};",100,0,1300);
    h_mTTbarElMu =fs->make<TH1F>("h_mTTbarElMu",";M^{tt};",100,0,1300);



    h_ALS_etaLMuMu = fs->make<TH1F>("h_ALS_EtaLepMuMu",";#eta_{l};",100,-3,3);
    h_ALS_etaLElEl = fs->make<TH1F>("h_ALS_EtaLepElEl",";#eta_{l};",100,-3,3);
    h_ALS_etaLElMu = fs->make<TH1F>("h_ALS_EtaLepElMu",";#eta_{l};",100,-3,3);
    h_ALS_etaLDiLep = fs->make<TH1F>("h_ALS_EtaLepDiLep",";#eta_{l};",100,-3,3);
    h_ALS_ptLMuMu = fs->make<TH1F>("h_ALS_PtLepMuMu",";Pt_{l};",100,0.,300.);
    h_ALS_ptLElMu = fs->make<TH1F>("h_ALS_PtLepElMu",";Pt_{l};",100,0.,300.);
    h_ALS_ptLElEl = fs->make<TH1F>("h_ALS_PtLepElEl",";Pt_{l};",100,0.,300.);
    h_ALS_ptLDiLep = fs->make<TH1F>("h_ALS_PtLepDiLep",";Pt_{l};",100,0.,300.);

    h_AJS_etaLepMuMu = fs->make<TH1F>("h_AJS_EtaLepMuMu",";#eta_{l};",100,-3,3);
    h_AJS_etaLepElEl = fs->make<TH1F>("h_AJS_EtaLepElEl",";#eta_{l};",100,-3,3);
    h_AJS_etaLepElMu = fs->make<TH1F>("h_AJS_EtaLepElMu",";#eta_{l};",100,-3,3);
    h_AJS_etaLepDiLep = fs->make<TH1F>("h_AJS_EtaLepDiLep",";#eta_{l};",100,-3,3);
    h_AJS_ptLepMuMu = fs->make<TH1F>("h_AJS_PtLepMuMu",";Pt_{l};",100,0.,300.);
    h_AJS_ptLepElMu = fs->make<TH1F>("h_AJS_PtLepElMu",";Pt_{l};",100,0.,300.);
    h_AJS_ptLepElEl = fs->make<TH1F>("h_AJS_PtLepElEl",";Pt_{l};",100,0.,300.);
    h_AJS_ptLepDiLep = fs->make<TH1F>("h_AJS_PtLepDiLep",";Pt_{l};",100,0.,300.);

    h_ABS_etaLepMuMu = fs->make<TH1F>("h_ABS_EtaLepMuMu",";#eta_{l};",100,-3,3);
    h_ABS_etaLepElEl = fs->make<TH1F>("h_ABS_EtaLepElEl",";#eta_{l};",100,-3,3);
    h_ABS_etaLepElMu = fs->make<TH1F>("h_ABS_EtaLepElMu",";#eta_{l};",100,-3,3);
    h_ABS_etaLepDiLep = fs->make<TH1F>("h_ABS_EtaLepDiLep",";#eta_{l};",100,-3,3);
    h_ABS_ptLepMuMu = fs->make<TH1F>("h_ABS_PtLepMuMu",";Pt_{l};",100,0.,300.);
    h_ABS_ptLepElMu = fs->make<TH1F>("h_ABS_PtLepElMu",";Pt_{l};",100,0.,300.);
    h_ABS_ptLepElEl = fs->make<TH1F>("h_ABS_PtLepElEl",";Pt_{l};",100,0.,300.);
    h_ABS_ptLepDiLep = fs->make<TH1F>("h_ABS_PtLepDiLep",";Pt_{l};",100,0.,300.);
    h_ABS_mLepMuMu = fs->make<TH1F>("h_ABS_mLepMuMu",";M_{ll};",100,0.,300.);
    h_ABS_mLepElEl = fs->make<TH1F>("h_ABS_mLepElEl",";M_{ll};",100,0.,300.);
    h_ABS_mLepElMu = fs->make<TH1F>("h_ABS_mLepElMu",";M_{ll};",100,0.,300.);
    h_ABS_mLepDiLep = fs->make<TH1F>("h_ABS_mLepDiLep",";M_{ll};",100,0.,300.);
    h_ABS_METMuMu = fs->make<TH1F>("h_ABS_METMuMu",";MET;",100,0.,300.);
    h_ABS_METElEl = fs->make<TH1F>("h_ABS_METElEl",";MET;",100,0.,300.);
    h_ABS_METElMu = fs->make<TH1F>("h_ABS_METElMu",";MET;",100,0.,300.);
    h_ABS_METDiLep = fs->make<TH1F>("h_ABS_METDiLep",";MET;",100,0.,300.);
    h_ABS_NBJetsMuMu = fs->make<TH1F>("h_ABS_NBJetsMuMu",";N_{BJets};",30,0,8);
    h_ABS_NBJetsElEl = fs->make<TH1F>("h_ABS_NBJetsElEl",";N_{BJets};",30,0,8);
    h_ABS_NBJetsElMu = fs->make<TH1F>("h_ABS_NBJetsElMu",";N_{BJets};",30,0,8);
    h_ABS_NBJetsDiLep = fs->make<TH1F>("h_ABS_NBJetsDiLep",";N_{BJets};",30,0,8);
    h_ABS_NJetsMuMu = fs->make<TH1F>("h_ABS_NJetsMuMu",";N_{BJets};",30,0,8);
    h_ABS_NJetsElEl = fs->make<TH1F>("h_ABS_NJetsElEl",";N_{BJets};",30,0,8);
    h_ABS_NJetsElMu = fs->make<TH1F>("h_ABS_NJetsElMu",";N_{BJets};",30,0,8);
    h_ABS_NJetsDiLep = fs->make<TH1F>("h_ABS_NJetsDiLep",";N_{BJets};",30,0,8);
    h_ABS_etaLeadingJetMuMu = fs->make<TH1F>("h_ABS_etaLeadingJetMuMu",";#eta_{l};",100,-3,3);
    h_ABS_etaLeadingJetElEl = fs->make<TH1F>("h_ABS_EtaLeadingJetElEl",";#eta_{l};",100,-3,3);
    h_ABS_etaLeadingJetElMu = fs->make<TH1F>("h_ABS_EtaLeadingJetElMu",";#eta_{l};",100,-3,3);
    h_ABS_etaLeadingJetDiLep = fs->make<TH1F>("h_ABS_EtaLeadingJetDiLep",";#eta_{l};",100,-3,3);
    h_ABS_ptLeadingJetMuMu = fs->make<TH1F>("h_ABS_PtLeadingJetMuMu",";Pt_{l};",100,0.,300.);
    h_ABS_ptLeadingJetElMu = fs->make<TH1F>("h_ABS_PtLeadingJetElMu",";Pt_{l};",100,0.,300.);
    h_ABS_ptLeadingJetElEl = fs->make<TH1F>("h_ABS_PtLeadingJetElEl",";Pt_{l};",100,0.,300.);
    h_ABS_ptLeadingJetDiLep = fs->make<TH1F>("h_ABS_ptLeadingJetDiLep",";Pt_{l};",100,0.,300.);

    h_AMS_ptLepMuMu = fs->make<TH1F>("h_AMS_PtLepMuMu",";Pt_{l};",100,0.,300.);
    h_AMS_ptLepElMu = fs->make<TH1F>("h_AMS_PtLepElMu",";Pt_{l};",100,0.,300.);
    h_AMS_ptLepElEl = fs->make<TH1F>("h_AMS_PtLepElEl",";Pt_{l};",100,0.,300.);
    h_AMS_ptLepDiLep = fs->make<TH1F>("h_AMS_PtLepDiLep",";Pt_{l};",100,0.,300.);
    h_AMS_mLepMuMu = fs->make<TH1F>("h_AMS_mLepMuMu",";M_{ll};",100,0.,300.);
    h_AMS_mLepElEl = fs->make<TH1F>("h_AMS_mLepElEl",";M_{ll};",100,0.,300.);
    h_AMS_mLepElMu = fs->make<TH1F>("h_AMS_mLepElMu",";M_{ll};",100,0.,300.);
    h_AMS_mLepDiLep = fs->make<TH1F>("h_AMS_mLepDiLep",";M_{ll};",100,0.,300.);
    h_AMS_METMuMu = fs->make<TH1F>("h_AMS_METMuMu",";MET;",100,0.,300.);
    h_AMS_METElEl = fs->make<TH1F>("h_AMS_METElEl",";MET;",100,0.,300.);
    h_AMS_METElMu = fs->make<TH1F>("h_AMS_METElMu",";MET;",100,0.,300.);
    h_AMS_METDiLep = fs->make<TH1F>("h_AMS_METDiLep",";MET;",100,0.,300.);

    h_truthRecoMuMu = fs->make<TH2F>("h_truthRecoMuMu","",400,-1.,1.,400,-1.,1.);
    h_truthRecoCos = fs->make<TH2F>("h_truthRecoCos","",400,-1.,1.,400,-1.,1.);
    h_truthRecoMTT = fs->make<TH2F>("h_truthRecoMTT","",400,0.,1000.,400,0.,1000.);
    h_truthMinusRecoCos = fs->make<TH1F>("h_truthMinusRecoCos","",100,-2,2);
    h_truthMinusRecoMTT = fs->make<TH1F>("h_truthMinusRecoMTT","",100,-400,400);



    //initialize the tree
    f_outFile->cd();
    t_outTree =  new TTree("tree","tr");
    t_outTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",&b_Mu1);
    t_outTree->Branch("SlimmedMuon_Pt",&SlimmedMuon_Pt);
//    t_outTree->Branch("TotalNumberOfEvents",&NEvent,"TotalNumberOfEvents/I");
//    t_outTree->Branch("NGoodvtx",&nGoodVtxs,"NGoodvtx/I");
//    t_outTree->Branch("RecoCos",&RecoCos);
//    if(!isData)
//    {
//        t_outTree->Branch("TruthCos",&TruthCos);
//        t_outTree->Branch("TruthMTT",&TruthMTT);
//        t_outTree->Branch("TruthCosMuMu",&TruthCosMuMu);
//        t_outTree->Branch("TruthMTTMuMu",&TruthMTTMuMu);
//    }
//    t_outTree->Branch("RecoMTT",&RecoMTT);
//    t_outTree->Branch("RecoCosMuMu",&RecoCosMuMu);

//    t_outTree->Branch("RecoMTTMuMu",&RecoMTTMuMu);
//    t_outTree->Branch("EvantsAfterVert",&n_afterVertex,"EvantsAfterVert/I");
//    t_outTree->Branch("EvantsAfterHLT",&n_afterHLT,"EvantsAfterHLT/I");
//    t_outTree->Branch("EvantsAfterDiLep",&n_afterDiLepton,"EvantsAfterDiLep/I");
//    t_outTree->Branch("EvantsAfterDiMu",&n_afterDiMu,"EvantsAfterDiMu/I");
//    t_outTree->Branch("EvantsAfterDiEl",&n_afterDiEl,"EvantsAfterDiEl/I");
//    t_outTree->Branch("EvantsAfterElMu",&n_afterElMu,"EvantsAfterElMu/I");
//    t_outTree->Branch("EvantsAfter2Jets",&n_after2Jets,"EvantsAfter2Jets/I");
//    t_outTree->Branch("EvantsAfter2BJets",&n_after2BJets,"EvantsAfter2BJets/I");
//    t_outTree->Branch("EvantsAfterMet",&n_afterMet,"EvantsAfterMet/I");
//    t_outTree->Branch("EvantsAfterTop",&n_afterTop,"EvantsAfterTop/I");



}
void MiniAnalyzer::endJob() {
    f_outFile->cd();

    t_outTree->Write();
    f_outFile->Close();
    amwtSolver->writeOut();
}


//____________________________________________________________________________
bool MiniAnalyzer::IsSoftMuon(const pat::Muon& mu , const reco::Vertex& vertex) {
    if (!(muon::isGoodMuon(mu, muon::TMOneStationTight))) return false;
    if (!(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)) return false;
    if (!(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0)) return false;
    if (!(mu.innerTrack()->quality(reco::TrackBase::highPurity))) return false;
    if (!((fabs(mu.innerTrack()->dxy(vertex.position())) < 0.3) && (fabs(mu.innerTrack()->dz(vertex.position())) < 20.))) return false;
    return true;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//____________________________________________________________________________
void MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
void MiniAnalyzer::beginRun(const Run & iRun, const EventSetup &){


}
void MiniAnalyzer::endRun(const Run & iRun, const EventSetup &){
    if(!isData){
    edm::Handle<LHERunInfoProduct> run;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    iRun.getByToken(lheInfo_,run);
//    iRun.getByLabel( "externalLHEProducer", run );
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());

    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      std::cout << iter->tag() << std::endl;
      std::vector<std::string> lines = iter->lines();
      for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
       std::cout << lines.at(iLine);
      }
    }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
