#include "MiniAnalyzer.h"



MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(mayConsume<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(mayConsume<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(mayConsume<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
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
    ptResData( (iConfig.getParameter<string>("ptResData"))),
    phiResData( (iConfig.getParameter<string>("phiResData")) ),
    sfResData( (iConfig.getParameter<string>("sfResData")) ),
    outfileName((iConfig.getParameter<string>("outFileName"))),
    triggerResluts_(mayConsume<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
    lheInfo_(mayConsume<LHERunInfoProduct,InRun>(iConfig.getParameter<edm::InputTag>("externalLHEProducer"))),
    isPythia((iConfig.getParameter<bool>("isPythia"))),
    triggerFilters_(mayConsume<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerFilters"))),
    DiMu((iConfig.getParameter<bool>("DiMu"))),
    DiEl((iConfig.getParameter<bool>("DiEl"))),
    ElMu((iConfig.getParameter<bool>("ElMu"))),
    egammaSF( (iConfig.getParameter<string>("egammaSF"))),
    isRunGH((iConfig.getParameter<bool>("isRunGH"))),
    muonISOSF( (iConfig.getParameter<string>("muonISOSF"))),
    muonIDSF( (iConfig.getParameter<string>("muonIDSF"))),
    egammaTkSF( (iConfig.getParameter<string>("egammaTkSF"))),
    btagSf( (iConfig.getParameter<string>("btagSf"))),
    muonTkSF( (iConfig.getParameter<string>("muonTkSF"))),
    isSingleMuon((iConfig.getParameter<bool>("isSingleMuon"))),
    isSingleElectron((iConfig.getParameter<bool>("isSingleElectron"))),
    isRunH((iConfig.getParameter<bool>("isRunH"))),
    eleIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMap"))),
    ee_sf( (iConfig.getParameter<string>("ee_sf"))),
    em_sf( (iConfig.getParameter<string>("em_sf"))),
    mm_sf( (iConfig.getParameter<string>("mm_sf"))),
    PileupSrc_(mayConsume<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PileupSrc"))),
    genJetToken_(mayConsume<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("SlimmedGenJets"))),
    s_pileup_data( (iConfig.getParameter<string>("pileup_data"))),
    s_pileup_mc( (iConfig.getParameter<string>("pileup_mc"))),
    eleIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >
                           (iConfig.getParameter<edm::InputTag>("eleIdFullInfoMap"))),
    verboseIdFlag_(iConfig.getParameter<bool>("eleIdVerbose"))





{
    beamSpotToken_    = consumes<reco::BeamSpot>
            (iConfig.getParameter <edm::InputTag>
             ("beamSpot"));
    conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >
            (iConfig.getParameter<edm::InputTag>
             ("conversionsMiniAOD"));

    // initializing the solver
    if(!isData)
    {
        amwtSolver = new TtAMWTSolver(isData,172.44,172.50,0.1,80.4,4.8,ptRes,phiRes,sfRes);
        ptResol  = new JME::JetResolution(ptRes);
        phiResol = new JME::JetResolution(phiRes);
        jerSF = new JME::JetResolutionScaleFactor(sfRes);
    }
    else
    {
        amwtSolver = new TtAMWTSolver(isData,172.44,172.50,0.1,80.4,4.8,ptResData,phiResData,sfResData);
        ptResol  = new JME::JetResolution(ptResData);
        phiResol = new JME::JetResolution(phiResData);
        jerSF = new JME::JetResolutionScaleFactor(sfResData);
    }
    //    std::vector<double> nupars= {30.7137,56.2880,23.0744,59.1015,24.9145};
    //    amwtsolver = new TtFullLepKinSolver(172.5,172.5,10,nupars,80.4,4.8);

    // Initializing  output root file;
    f_outFile = TFile::Open(outfileName.c_str(),"RECREATE");
    // cout << outfileName.c_str() << endl;
    if(!isData)
    {
        f_egammaSF = new TFile(egammaSF.c_str());
        f_egammaTkSF = new TFile(egammaTkSF.c_str());

        f_muonIDSF = new TFile(muonIDSF.c_str());
        f_muonISOSF = new TFile(muonISOSF.c_str());

        f_muonTkSF = new TFile(muonTkSF.c_str());
        f_ee_sf = new TFile(ee_sf.c_str());
        f_em_sf = new TFile(em_sf.c_str());
        f_mm_sf = new TFile(mm_sf.c_str());

        // setup calibration + reader

        BTagCalibration calib("csvv2", btagSf);
        reader = BTagCalibrationReader(BTagEntry::OP_LOOSE,  // operating point
                                       "central",             // central sys type
        {"up", "down"}
                                       );      // other sys types

        reader.load(calib,                // calibration instance
                    BTagEntry::FLAV_B,    // btag flavour
                    "comb");              // measurement type

    }
    if(!isData)
    {
        LumiWeights_ = edm::LumiReWeighting(s_pileup_mc,s_pileup_data,"pileup","pileup");

    }


}

MiniAnalyzer::~MiniAnalyzer()
{
    delete f_outFile;
}
double MiniAnalyzer::SF(double x)
{
    double SF = exp(0.0615 + (-0.0005)*x);
    return SF;
}

void MiniAnalyzer::printCutFlowResult(vid::CutFlowResult &cutflow)
{
    printf("    CutFlow name= %s    decision is %d\n",
           cutflow.cutFlowName().c_str(),
           (int) cutflow.cutFlowPassed());
    int ncuts = cutflow.cutFlowSize();
    printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
    for(int icut = 0; icut<ncuts; icut++){
        printf("  %2d      %50s    %d        %f          %d\n", icut,
               cutflow.getNameAtIndex(icut).c_str(),
               (int)cutflow.isCutMasked(icut),
               cutflow.getValueCutUpon(icut),
               (int)cutflow.getCutResultByIndex(icut));
    }
}
void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


    using namespace edm;
    h_Nevents->Fill(1);

    t_Run   = iEvent.id().run();
    t_Event = iEvent.id().event();
    t_lumi  = iEvent.luminosityBlock();
    t_bunch = iEvent.bunchCrossing();

    /////////////FILLINGINPUTS/////////////////////
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<edm::View<reco::GsfElectron>> electrons;
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
    // Get the electron ID data from the event stream.
    // Note: this implies that the VID ID modules have been run upstream.
    // If you need more info, check with the EGM group.
    edm::Handle<edm::ValueMap<bool> > ele_id_decisions;
    iEvent.getByToken(eleIdMapToken_,ele_id_decisions);
    // Get the beam spot
    edm::Handle<reco::BeamSpot> theBeamSpot;
    iEvent.getByToken(beamSpotToken_,theBeamSpot);
    edm::Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(conversionsMiniAODToken_, conversions);
    edm::Handle<edm::ValueMap<vid::CutFlowResult> > ele_id_cutflow_data;
    iEvent.getByToken(eleIdFullInfoMapToken_,ele_id_cutflow_data);
    t_num_PV = -999;
    t_PU_BunchCrossing =-999;
    t_num_PU_gen_vertices = -999;
    t_num_PU_vertices = -999;
    t_PU_weight = 1;
    t_PU_weight_secondWay = 1;
    edm::Handle<reco::GenJetCollection> genJets;

    if(!isData)
    {
        edm::Handle<vector<PileupSummaryInfo> > PupInfo;
        iEvent.getByToken(PileupSrc_, PupInfo);
        if(!PupInfo.isValid()) return;
        float Tnpv = -1;
        for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin(); cand != PupInfo->end(); ++ cand )
        {
            t_num_PU_vertices = cand->getTrueNumInteractions();

            t_PU_BunchCrossing = cand->getBunchCrossing();
            t_num_PU_gen_vertices = cand->getPU_NumInteractions();




            int BX = cand->getBunchCrossing();

            if(BX == 0) {
                Tnpv = cand->getTrueNumInteractions();
                continue;
            }

            t_PU_weight = LumiWeights_.weight( Tnpv );


        }
        //        const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
        //        t_PU_weight_secondWay = LumiWeights_.weight( (*iEventB) );
        iEvent.getByToken(prunedGenToken_,pruned);
        iEvent.getByToken(packedGenToken_,packed);
        iEvent.getByToken( genEvtInfo_, genEvtInfo );
        iEvent.getByToken( lheEvtInfo_, lhEvtInfo);
        iEvent.getByToken(ttgenEvt_, genEvent);
        genColl = *packed;

        iEvent.getByToken(genJetToken_, genJets);

    }

    edm::Handle<double> rhoHandle_;
    iEvent.getByToken( rho_,rhoHandle_);
    float rho = *rhoHandle_;
    t_rho = rho;
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
    SF_posMu = 1;
    SF_negMu = 1;
    SF_posEl = 1;
    SF_negEl = 1;
    w_bjets.clear();
    w_bjet = 1;
    w_posEl = 1;
    w_negEl = 1;
    w_posMu = 1;
    w_negMu = 1;
    w_ee = 1;
    w_em = 1;
    w_mm = 1;
    w_top = 1;
    w_mc = 1;



    ////Tree variables
    t_antinu.clear();
    t_nu.clear();
    t_antitop.clear();
    t_top.clear();
    t_antiW.clear();
    t_W.clear();
    t_bJets.clear();
    t_cos.clear();
    t_Leptons.clear();
    /////////////////////0b HLT trigger/////////////////////
    edm::Handle<edm::TriggerResults> trigResults; //our trigger result object

    iEvent.getByToken(triggerResluts_,trigResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
    std::string Mu1,Mu2,Mu3,Mu4,MuMu1,MuMu2,MuMu3,MuMu4,El1,El2,El3,ElEl1,ElEl2,ElMu1,ElMu2,ElMu3,ElMu4,ElMu5,ElMu6;

    if(isData)
    {
        // Double Muon  36.811 fb^-1 Mu1 or Mu2
        Mu1="HLT_IsoMu24_v";
        Mu2="HLT_IsoTkMu24_v";

        Mu3 = "HLT_IsoMu22_eta2p1_v";
        Mu4 = "HLT_IsoTkMu22_eta2p1_v";
        MuMu1 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v";
        MuMu2 = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v";
        MuMu3 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
        MuMu4 = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
        El1="HLT_Ele32_eta2p1_WPTight_Gsf_v";
        El2="HLT_Ele27_WPTight_Gsf_v";
        El3="HLT_Ele25_eta2p1_WPTight_Gsf_v";

        // Double Electron 36.615 fb^-1
        ElEl1="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

        // Double Electron 36.811 fb^-1
        ElEl2="HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v";



        // El Mu
        ElMu1="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
        ElMu2="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
        ElMu3="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
        ElMu4="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v";
        ElMu5="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
        ElMu6="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v";
        //cout << "safe here" << endl;
    }
    else
    {
        // Double Muon  36.811 fb^-1 Mu1 or Mu2
        Mu1="HLT_IsoMu24_v";
        Mu2="HLT_IsoTkMu24_v";

        Mu3 = "HLT_IsoMu22_eta2p1_v";
        Mu4 = "HLT_IsoTkMu22_eta2p1_v";
        MuMu1 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v";
        MuMu2 = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v";
        MuMu3 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
        MuMu4 = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
        El1="HLT_Ele32_eta2p1_WPTight_Gsf_v";
        El2="HLT_Ele27_WPTight_Gsf_v";
        El3="HLT_Ele25_eta2p1_WPTight_Gsf_v";

        // Double Electron 36.615 fb^-1
        ElEl1="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

        // Double Electron 36.811 fb^-1
        ElEl2="HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v";



        // El Mu
        ElMu1="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
        ElMu2="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
        ElMu3="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
        ElMu4="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v";
        ElMu5="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
        ElMu6="HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v";
        //cout << "safe here" << endl;


    }

    for (unsigned int i = 0, n = trigResults->size(); i < n; ++i)
    {
        std::string nameHLT,st;
        nameHLT = trigNames.triggerName(i);
        st = nameHLT.substr(0, nameHLT.size()-1);

        if(st.compare(Mu1) == 0)
        {
            b_Mu1 = trigResults->accept(i);
            // if(b_Mu1) cout << st << endl;
        }
        if(st.compare(Mu2) == 0)
        {
            b_Mu2 = trigResults->accept(i);
            //if(b_Mu2) cout << st << endl;
        }
        if(st.compare(Mu3) == 0)
        {
            b_Mu3 = trigResults->accept(i);
            //if(b_Mu3) cout << st << endl;
        }
        if(st.compare(Mu4) == 0)
        {
            b_Mu4 = trigResults->accept(i);
            //if(b_Mu4) cout << st << endl;
        }
        if(st.compare(MuMu1) == 0)
        {
            b_MuMu1 = trigResults->accept(i);
            //            cout << st << endl;
        }
        if(st.compare(MuMu2) == 0)
        {
            b_MuMu2 = trigResults->accept(i);
            //            cout << st << endl;
        }
        if(st.compare(MuMu3) == 0)
        {
            b_MuMu3 = trigResults->accept(i);
            //            cout << st << endl;
        }
        if(st.compare(MuMu4) == 0)
        {
            b_MuMu4 = trigResults->accept(i);
            //            cout << st << endl;
        }
        if(st.compare(ElEl1) == 0)
        {
            b_ElEl1 = trigResults->accept(i);
            //            cout << st << endl;
        }
        if(st.compare(ElEl2) == 0)
        {
            b_ElEl2 = trigResults->accept(i);
            //            cout << st << endl;
        }
        if(st.compare(El1) == 0)
        {
            b_El1 = trigResults->accept(i);
            //if(b_El1) cout << st << endl;
        }
        if(st.compare(El2) == 0)
        {
            b_El2 = trigResults->accept(i);
            //if(b_El2) cout << st << endl;
        }
        if(st.compare(El3) == 0)
        {
            b_El3 = trigResults->accept(i);
            //if(b_El3) cout << st << endl;
        }
        if(st.compare(ElMu1) == 0)
        {
            b_ElMu1 = trigResults->accept(i);
            //if(b_ElMu1) cout << st << endl;
        }
        if(st.compare(ElMu2) == 0)
        {
            b_ElMu2 = trigResults->accept(i);
            //if(b_ElMu2) cout << st << endl;
        }
        if(st.compare(ElMu3) == 0)
        {
            b_ElMu3 = trigResults->accept(i);
            //if(b_ElMu3) cout << st << endl;
        }
        if(st.compare(ElMu4) == 0)
        {
            b_ElMu4 = trigResults->accept(i);
            //if(b_ElMu4) cout << st << endl;
        }
        if(st.compare(ElMu5) == 0)
        {
            b_ElMu5 = trigResults->accept(i);
            //if(b_ElMu5) cout << st << endl;
        }
        if(st.compare(ElMu6) == 0)
        {
            b_ElMu6 = trigResults->accept(i);
            //if(b_ElMu6) cout << st << endl;
        }
        //////////////PRINT ALL TRIGGER PASSES /////////////////////
        //        std::cout << "Trigger " << trigNames.triggerName(i) << "i: "<<i

        //                  <<": " << (trigResults->accept(i) ? "PASS" : "fail (or not run)")
        //                 << std::endl;
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

        if(!isPythia)
        {

            float EventWeight = 1.0;
            if(!genEvtInfo.isValid()) return;
            EventWeight = genEvtInfo->weight();
            //std::cout<<"mc_weight = "<< genEvtInfo->weight() <<std::endl;
            float mc_weight = ( EventWeight > 0 ) ? 1 : -1;
            //std::cout<<"mc_weight = "<< mc_weight <<std::endl;
            theWeight = mc_weight;
        }
        //int whichWeight = 1001;
        //if(!isPythia) theWeight *= lhEvtInfo->weights()[whichWeight].wgt/lhEvtInfo->originalXWGTUP();

    }
    w_mc = theWeight;
    ////////////0a total number of events /////////////////
    h_Nevents_weighted->Fill(1,theWeight);

    if(!isData)
    {
        if(ElMu && !(b_Mu1 || b_Mu2 || b_El2 || b_ElMu1 || b_ElMu2 || b_ElMu3 || b_ElMu4 )) return;
        if(DiEl && !(b_El2 || b_ElEl1)) return;
        if(DiMu && !(b_Mu1 || b_Mu2  || b_MuMu1 || b_MuMu3 || b_MuMu4 ) ) return;

    }
    if(DiEl && isData)
    {
        if(!isSingleElectron && !isSingleMuon)
        {

            if(!b_ElEl1) return;

        }
        if(isSingleElectron)
        {


            if(!b_El2 || b_ElEl1 )return;

        }

    }

    if(DiMu && isData)
    {
        if(!isSingleElectron && !isSingleMuon)
        {
            if(isRunH)
            {
                if(!(b_MuMu3 || b_MuMu4 ) ) return;
            }
            else
            {
                if(!(b_MuMu1 ) ) return;
            }
        }

        if(isSingleMuon)
        {

            if(isRunH)
            {
                if(!(b_Mu1|| b_Mu2) || (  b_MuMu3 || b_MuMu4 ) ) return;
            }
            else
            {
                if(!(b_Mu1|| b_Mu2) ||  b_MuMu1    ) return;

            }
        }
    }

    if(ElMu && isData)
    {
        if(!isSingleElectron && !isSingleMuon)
        {
            if(isRunH)
            {
                if(!( b_ElMu2  || b_ElMu4 )) return;
            }
            else
            {
                if(!( b_ElMu1  || b_ElMu3  )) return;
            }
        }
        if(isSingleElectron)
        {

            if(isRunH)
            {
                if(!b_El2 || (  b_ElMu2  || b_ElMu4) )return;
            }
            else
            {
                if(!b_El2 || (b_ElMu1 || b_ElMu3 ) )return;

            }
        }
        if(isSingleMuon)
        {

            if(isRunH)
            {
                if(!(b_Mu1 || b_Mu2)   || (b_El2  || b_ElMu2 || b_ElMu4) )return;
            }
            else
            {
                if(!(b_Mu1 || b_Mu2)   || (b_El2 || b_ElMu1 || b_ElMu3 ) )return;

            }
        }
    }
    ///Number of events after trigger
    h_Nevents_AT->Fill(1,theWeight);;

    /////////////////////0c Event Filters //////////////////


    ///Vertex Filter
    if (vertices->empty()) return;
    VertexCollection::const_iterator PV = vertices->begin();
    bool isFake = (PV->chi2()==0 && PV->ndof()==0);
    if(isFake) return;
    if(!(PV->ndof() >= 4.)) return;
    if(!(fabs(PV->position().Rho()) < 2.0)) return;
    if(!(fabs(PV->position().Z()) < 24.0)) return;
    if (vertices->empty()) return;

    //        VertexCollection::const_iterator PV = vertices->end();
    //        for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV)
    //        {
    //            bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    //            if ( !(isFake)
    //                    && vtx->ndof() >= 4. && fabs(vtx->position().Rho()) < 2.0
    //                    && fabs(vtx->position().Z()) < 24.0)
    //           {
    //               PV = vtx;
    //               break;
    //            }
    //       }
    //       if (PV==vertices->end()) return;





    ///Noise filters
    edm::Handle<edm::TriggerResults> triggerFilters; // Flag and filters
    iEvent.getByToken(triggerFilters_,triggerFilters);
    const edm::TriggerNames& filterNames = iEvent.triggerNames(*triggerFilters);
    for (unsigned int i = 0, n = triggerFilters->size(); i < n; ++i)
    {
        std::string nameFilter;
        nameFilter = filterNames.triggerName(i);
        if(nameFilter.compare("Flag_globalTightHalo2016Filter") == 0)
        {
            Flag_globalTightHalo2016Filter = triggerFilters->accept(i);
            //            cout << nameFilter <<endl;
        }
        if(nameFilter.compare("Flag_HBHENoiseFilter") == 0)
        {
            Flag_HBHENoiseFilter = triggerFilters->accept(i);
            //            cout << nameFilter <<endl;
        }
        if(nameFilter.compare("Flag_HBHENoiseIsoFilter") == 0)
        {
            Flag_HBHENoiseIsoFilter = triggerFilters->accept(i);
            //            cout << nameFilter <<endl;
        }
        if(nameFilter.compare("Flag_EcalDeadCellTriggerPrimitiveFilter") == 0)
        {
            Flag_EcalDeadCellTriggerPrimitiveFilter = triggerFilters->accept(i);
            //            cout << nameFilter <<endl;
        }
        if(nameFilter.compare("Flag_BadPFMuonFilter") == 0)
        {
            Flag_BadPFMuonFilter = triggerFilters->accept(i);
            //            cout << nameFilter <<endl;
        }
        if(nameFilter.compare("Flag_BadChargedCandidateFilter") == 0)
        {
            Flag_BadChargedCandidateFilter = triggerFilters->accept(i);
            //            cout << nameFilter <<endl;
        }
        if(nameFilter.compare("Flag_eeBadScFilter") == 0)
        {
            Flag_eeBadScFilter = triggerFilters->accept(i);
            //            cout << nameFilter <<endl;
        }

        //        if(nameFilter.find("Flag") != std::string::npos) cout << nameFilter <<endl;
    }

    if(!( Flag_globalTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter )) return;


    if(isData && !(Flag_eeBadScFilter) ) return; /*&& Flag_BadChargedCandidateFilter && Flag_BadPFMuonFilter*/





    ///Number of events after vertex filter
    h_Nevents_AVS->Fill(1,theWeight);;


    ///////////////1 Dilepton pair choice //////

    //////////////MUON////////////////
    /// finding pos and neg muon with highest pt
    /// tight muon

    pat::GenericParticle posMu;
    pat::GenericParticle negMu;
    pat::GenericParticle posEl;
    pat::GenericParticle negEl;
    vector<pat::GenericParticle> MyMuons;
    vector<pat::GenericParticle> MyElectrons;
    vector<MyLepton> MyLeptons;


    for (pat::MuonCollection::const_iterator mup = muons->begin(); mup != muons->end(); ++mup)
    {

        //        if( !(mup->charge() > 0)) continue;
        if( !(mup->pt() > 20.0 )) continue;
        if( !(fabs(mup->eta()) <= 2.4 )) continue;


        if(!(mup->isTightMuon(*PV))) continue;
        coriso = 9999;
        //        cout << "isolation before " <<coriso <<endl;
        if(!(mup->isIsolationValid())) continue;
        if( mup->isIsolationValid())
        {
            reco::MuonPFIsolation pfR04 = mup->pfIsolationR04();
            coriso = pfR04.sumChargedHadronPt + std::max(0., pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
        }
        //        cout << "isolation after " <<coriso <<endl;

        if (!(coriso/mup->pt() <  0.15)) continue;
        MyMuons.push_back(*mup);
        MyLepton lep;
        lep.pat_particle = *mup;
        lep.flavour = 13;
        lep.type = "muon";
        MyLeptons.push_back(lep);
        //        if(mup->pt() > posMu.pt() ) posMu = *mup;
        //        // cout << mup->pt() <<endl;
        //        // cout << posMu.pt() << endl;
        h_PtMu->Fill(mup->pt(),theWeight);

        h_etaMu->Fill(mup->eta(),theWeight);
        //        isPosMu = true;

    }



    //    // cout << posMu.charge() << " charges " << negMu.charge() << endl;
    ////////////////////ELECTRONS//////////////////////////
    /// electron identification
    /////tight electrons

    for (size_t i = 0; i < electrons->size(); ++i)
    {

        const auto elp = electrons->ptrAt(i);

        if( !( elp->pt() > 20 ) ) continue;
        if( !( fabs(elp->eta() ) < 2.4)) continue;
        // Conversion rejection
        bool passConvVeto = !ConversionTools::hasMatchedConversion(*elp,
                                                                   conversions,
                                                                   theBeamSpot->position());
        if(!passConvVeto) continue;
        bool isPassEleId = (*ele_id_decisions)[elp];
        if(!isPassEleId)continue;
        if(fabs(elp->superCluster()->eta()) <= 1.479)
        {
            if(fabs(elp->superCluster()->eta()) >= 1.4442) continue;
            //impact parameters
            if( !( fabs((-1)*elp->gsfTrack()->dxy(PV->position())) < 0.05 )) continue;
            if(!( fabs(elp->gsfTrack()->dz(PV->position())) < 0.10)) continue;
        }
        if(fabs(elp->superCluster()->eta()) > 1.479)
        {
            if(fabs(elp->superCluster()->eta()) <= 1.5660) continue;
            //impact parameters
            if( !( fabs((-1)*elp->gsfTrack()->dxy(PV->position()))  < 0.1 )) continue;
            if(!( fabs(elp->gsfTrack()->dz(PV->position())) < 0.20)) continue;
        }
        if(verboseIdFlag_){
            vid::CutFlowResult fullCutFlowData = (*ele_id_cutflow_data)[elp];
            //
            // Full printout
            //
            printf("\nDEBUG CutFlow, full info for cand with pt=%f:\n", elp->pt());
            printCutFlowResult(fullCutFlowData);
            //
            // Example of how to find the ID decision with one cut removed,
            // this could be needed for N-1 studies.
            //
            const int cutIndexToMask = 4;
            // Here we masked the cut by cut index, but you can also do it by cut name string.
            vid::CutFlowResult maskedCutFlowData
                    = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
            printf("DEBUG CutFlow, the result with cut %s masked out\n",
                   maskedCutFlowData.getNameAtIndex(cutIndexToMask).c_str());
            printCutFlowResult(maskedCutFlowData);
        }
        MyElectrons.push_back(*elp);
        MyLepton lep;
        lep.pat_particle = *elp;
        lep.flavour = 11;
        lep.type = "electron";
        MyLeptons.push_back(lep);
    }

    bool my_mumu = false;
    bool my_ee = false;
    bool my_emu = false;

    t_Leptons.clear();
    t_Leptons_e.clear();
    t_Leptons_eta.clear();
    t_Leptons_phi.clear();
    t_Leptons_pt.clear();
    t_Leptons_charge.clear();
    t_Leptons_etaSC.clear();
    t_Leptons_d0.clear();
    t_Leptons_dz.clear();
    t_Leptons_type.clear();

    t_gen_Leptons_e.clear();
    t_gen_Leptons_eta.clear();
    t_gen_Leptons_phi.clear();
    t_gen_Leptons_pt.clear();
    t_gen_Leptons_id.clear();
    t_gen_Leptons_charge.clear();

    t_gen_daughters_pt.clear();
    t_gen_daughters_phi.clear();
    t_gen_daughters_eta.clear();
    t_gen_daughters_e.clear();
    t_gen_daughters_id.clear();

    t_gen_mothers_e.clear();
    t_gen_mothers_eta.clear();
    t_gen_mothers_id.clear();
    t_gen_mothers_phi.clear();
    t_gen_mothers_pt.clear();


bool b_thisway = true;
if(b_thisway){
        MyLepton lep1;
        MyLepton lep2;
        std::sort(MyLeptons.begin(),MyLeptons.end(),by_pt());
        if(MyLeptons.size() > 0 ){
            lep1 = MyLeptons[0];
        }
        if(MyLeptons.size() > 1){
            lep2 = MyLeptons[1];
        }
        if(ElMu)
        {
            if(MyLeptons.size()>=1 && lep1.pat_particle.pt() > 25 && ((lep1.flavour == 11 && lep2.flavour == 13) || (lep2.flavour == 11 && lep1.flavour == 13))
                    && (lep1.pat_particle.p4() + lep2.pat_particle.p4()).M() > 20 && lep1.pat_particle.charge()*lep2.pat_particle.charge() == -1 ) my_emu = true ;
        }
        else if(DiEl)
        {
            if(MyLeptons.size()>=1 && lep1.pat_particle.pt() > 25 && ((lep1.flavour == 11 && lep2.flavour == 11) || (lep2.flavour == 11 && lep1.flavour == 11))
                    && (lep1.pat_particle.p4() + lep2.pat_particle.p4()).M() > 20 && lep1.pat_particle.charge()*lep2.pat_particle.charge() == -1 ) my_ee = true ;
        }
        else if(DiMu)
        {
            if(MyLeptons.size()>=1 && lep1.pat_particle.pt() > 25 && ((lep1.flavour == 13 && lep2.flavour == 13) || (lep2.flavour == 13 && lep1.flavour == 13))
                    && (lep1.pat_particle.p4() + lep2.pat_particle.p4()).M() > 20 && lep1.pat_particle.charge()*lep2.pat_particle.charge() == -1 ) my_mumu = true ;
        }
        else
        {
            return;
        }
        if(my_ee)
        {
            if(lep1.pat_particle.charge() > 0)
            {
                posEl = lep1.pat_particle;
                isPosEl = true;
                t_Leptons.push_back(posEl.p4());

                t_Leptons_pt.push_back(posEl.pt());
                t_Leptons_eta.push_back(posEl.eta());
                t_Leptons_phi.push_back(posEl.phi());
                t_Leptons_e.push_back(posEl.energy());
                t_Leptons_charge.push_back(posEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(posEl.superCluster()->eta());
                t_Leptons_d0.push_back(posEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(posEl.gsfTrack()->dz(PV->position()));


                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());
                }
                negEl = lep2.pat_particle;
                isNegEl = true;
                t_Leptons.push_back(negEl.p4());

                t_Leptons_pt.push_back(negEl.pt());
                t_Leptons_eta.push_back(negEl.eta());
                t_Leptons_phi.push_back(negEl.phi());
                t_Leptons_e.push_back(negEl.energy());
                t_Leptons_charge.push_back(negEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(negEl.superCluster()->eta());
                t_Leptons_d0.push_back(negEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(negEl.gsfTrack()->dz(PV->position()));

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }

            }
            else
            {

                posEl = lep2.pat_particle;
                isPosEl = true;
                t_Leptons.push_back(posEl.p4());

                t_Leptons_pt.push_back(posEl.pt());
                t_Leptons_eta.push_back(posEl.eta());
                t_Leptons_phi.push_back(posEl.phi());
                t_Leptons_e.push_back(posEl.energy());
                t_Leptons_charge.push_back(posEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(posEl.superCluster()->eta());
                t_Leptons_d0.push_back(posEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(posEl.gsfTrack()->dz(PV->position()));


                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());
                }
                negEl = lep1.pat_particle;
                isNegEl = true;
                t_Leptons.push_back(negEl.p4());

                t_Leptons_pt.push_back(negEl.pt());
                t_Leptons_eta.push_back(negEl.eta());
                t_Leptons_phi.push_back(negEl.phi());
                t_Leptons_e.push_back(negEl.energy());
                t_Leptons_charge.push_back(negEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(negEl.superCluster()->eta());
                t_Leptons_d0.push_back(negEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(negEl.gsfTrack()->dz(PV->position()));

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }


            }

        }
        else if(my_mumu)
        {
            if(lep1.pat_particle.charge() > 0)
            {
                posMu = lep1.pat_particle;
                isPosMu = true;
                t_Leptons.push_back(posMu.p4());

                t_Leptons_pt.push_back(posMu.pt());
                t_Leptons_eta.push_back(posMu.eta());
                t_Leptons_phi.push_back(posMu.phi());
                t_Leptons_e.push_back(posMu.energy());
                t_Leptons_charge.push_back(posMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());

                }
                negMu = lep2.pat_particle;
                isNegMu = true;
                t_Leptons.push_back(negMu.p4());

                t_Leptons_pt.push_back(negMu.pt());
                t_Leptons_eta.push_back(negMu.eta());
                t_Leptons_phi.push_back(negMu.phi());
                t_Leptons_e.push_back(negMu.energy());
                t_Leptons_charge.push_back(negMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }

            }
            else
            {
                posMu = lep2.pat_particle;
                isPosMu = true;
                t_Leptons.push_back(posMu.p4());

                t_Leptons_pt.push_back(posMu.pt());
                t_Leptons_eta.push_back(posMu.eta());
                t_Leptons_phi.push_back(posMu.phi());
                t_Leptons_e.push_back(posMu.energy());
                t_Leptons_charge.push_back(posMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());

                }
                negMu = lep1.pat_particle;
                isNegMu = true;
                t_Leptons.push_back(negMu.p4());

                t_Leptons_pt.push_back(negMu.pt());
                t_Leptons_eta.push_back(negMu.eta());
                t_Leptons_phi.push_back(negMu.phi());
                t_Leptons_e.push_back(negMu.energy());
                t_Leptons_charge.push_back(negMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }
            }
        }
        else if(my_emu)
        {
            if(lep1.pat_particle.charge() > 0 && lep1.flavour == 11)
            {
                posEl = lep1.pat_particle;
                isPosEl = true;
                t_Leptons.push_back(posEl.p4());

                t_Leptons_pt.push_back(posEl.pt());
                t_Leptons_eta.push_back(posEl.eta());
                t_Leptons_phi.push_back(posEl.phi());
                t_Leptons_e.push_back(posEl.energy());
                t_Leptons_charge.push_back(posEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(posEl.superCluster()->eta());
                t_Leptons_d0.push_back(posEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(posEl.gsfTrack()->dz(PV->position()));


                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());
                }
                negMu = lep2.pat_particle;
                isNegMu = true;
                t_Leptons.push_back(negMu.p4());

                t_Leptons_pt.push_back(negMu.pt());
                t_Leptons_eta.push_back(negMu.eta());
                t_Leptons_phi.push_back(negMu.phi());
                t_Leptons_e.push_back(negMu.energy());
                t_Leptons_charge.push_back(negMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }


            }
            else if(lep1.pat_particle.charge() > 0 && lep1.flavour == 13)
            {
                posMu = lep1.pat_particle;
                isPosMu = true;
                t_Leptons.push_back(posMu.p4());

                t_Leptons_pt.push_back(posMu.pt());
                t_Leptons_eta.push_back(posMu.eta());
                t_Leptons_phi.push_back(posMu.phi());
                t_Leptons_e.push_back(posMu.energy());
                t_Leptons_charge.push_back(posMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());

                }
                negEl = lep2.pat_particle;
                isNegEl = true;
                t_Leptons.push_back(negEl.p4());

                t_Leptons_pt.push_back(negEl.pt());
                t_Leptons_eta.push_back(negEl.eta());
                t_Leptons_phi.push_back(negEl.phi());
                t_Leptons_e.push_back(negEl.energy());
                t_Leptons_charge.push_back(negEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(negEl.superCluster()->eta());
                t_Leptons_d0.push_back(negEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(negEl.gsfTrack()->dz(PV->position()));

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }



            }
            else if(lep2.pat_particle.charge() > 0 && lep2.flavour == 13)
            {
                posMu = lep2.pat_particle;
                isPosMu = true;
                t_Leptons.push_back(posMu.p4());

                t_Leptons_pt.push_back(posMu.pt());
                t_Leptons_eta.push_back(posMu.eta());
                t_Leptons_phi.push_back(posMu.phi());
                t_Leptons_e.push_back(posMu.energy());
                t_Leptons_charge.push_back(posMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());

                }
                negEl = lep1.pat_particle;
                isNegEl = true;
                t_Leptons.push_back(negEl.p4());

                t_Leptons_pt.push_back(negEl.pt());
                t_Leptons_eta.push_back(negEl.eta());
                t_Leptons_phi.push_back(negEl.phi());
                t_Leptons_e.push_back(negEl.energy());
                t_Leptons_charge.push_back(negEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(negEl.superCluster()->eta());
                t_Leptons_d0.push_back(negEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(negEl.gsfTrack()->dz(PV->position()));

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }



            }
            if(lep2.pat_particle.charge() > 0 && lep2.flavour == 11)
            {
                posEl = lep2.pat_particle;
                isPosEl = true;
                t_Leptons.push_back(posEl.p4());

                t_Leptons_pt.push_back(posEl.pt());
                t_Leptons_eta.push_back(posEl.eta());
                t_Leptons_phi.push_back(posEl.phi());
                t_Leptons_e.push_back(posEl.energy());
                t_Leptons_charge.push_back(posEl.charge());
                t_Leptons_type.push_back("electron");
                t_Leptons_etaSC.push_back(posEl.superCluster()->eta());
                t_Leptons_d0.push_back(posEl.gsfTrack()->dxy(PV->position()));
                t_Leptons_dz.push_back(posEl.gsfTrack()->dz(PV->position()));


                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posEl, genColl,11);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());
                }
                negMu = lep1.pat_particle;
                isNegMu = true;
                t_Leptons.push_back(negMu.p4());

                t_Leptons_pt.push_back(negMu.pt());
                t_Leptons_eta.push_back(negMu.eta());
                t_Leptons_phi.push_back(negMu.phi());
                t_Leptons_e.push_back(negMu.energy());
                t_Leptons_charge.push_back(negMu.charge());
                t_Leptons_type.push_back("muon");

                if(!isData)
                {
                    const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negMu, genColl,13);

                    const reco::GenParticle daughter = mse::getMotherPacked(matched);
                    const reco::GenParticle mom = mse::getMother(daughter);

                    t_gen_Leptons_pt.push_back(matched.pt());
                    t_gen_Leptons_eta.push_back(matched.eta());
                    t_gen_Leptons_phi.push_back(matched.phi());
                    t_gen_Leptons_e.push_back(matched.energy());
                    t_gen_Leptons_id.push_back(matched.pdgId());
                    t_gen_Leptons_charge.push_back(matched.charge());

                    t_gen_daughters_pt.push_back(daughter.pt());
                    t_gen_daughters_eta.push_back(daughter.eta());
                    t_gen_daughters_phi.push_back(daughter.phi());
                    t_gen_daughters_e.push_back(daughter.energy());
                    t_gen_daughters_id.push_back(daughter.pdgId());

                    t_gen_mothers_pt.push_back(mom.pt());
                    t_gen_mothers_eta.push_back(mom.eta());
                    t_gen_mothers_phi.push_back(mom.phi());
                    t_gen_mothers_e.push_back(mom.energy());
                    t_gen_mothers_id.push_back(mom.pdgId());



                }


            }

        }
        else
        {
            return;
        }
}
else{
    //cout << MyElectrons.size() << " Electron size" <<endl;

    if((MyMuons.size() + MyElectrons.size()) >=2)
    {
        // select leptons
        if(MyElectrons.size() == 0) my_mumu = true;
        else if (MyMuons.size() == 0) my_ee = true;
        else if (MyElectrons.size() == 1 )
        {
            if(MyMuons.size() == 1) my_emu = true;
            else if (MyElectrons[0].pt() > MyMuons[1].pt()) my_emu = true;
            else my_mumu = true;
        }
        else if (MyElectrons.size() > 1)
        {
            if(MyElectrons[1].pt() > MyMuons[0].pt() ) my_ee = true;
            else if (MyMuons.size() == 1) my_emu = true;
            else if (MyElectrons[0].pt() > MyMuons[1].pt() ) my_emu = true;
            else my_mumu = true;
        }
        if(my_ee)
        {
            if(MyElectrons[0].charge() != MyElectrons[1].charge())
            {
                isDiElectron = true;
                if(MyElectrons[0].charge() > 0)
                {
                    posEl = MyElectrons[0];
                    isPosEl = true;
                    t_Leptons.push_back(posEl.p4());

                    t_Leptons_pt.push_back(posEl.pt());
                    t_Leptons_eta.push_back(posEl.eta());
                    t_Leptons_phi.push_back(posEl.phi());
                    t_Leptons_e.push_back(posEl.energy());
                    t_Leptons_charge.push_back(posEl.charge());
                    t_Leptons_type.push_back("electron");
                    t_Leptons_etaSC.push_back(posEl.superCluster()->eta());
                    t_Leptons_d0.push_back(posEl.gsfTrack()->dxy(PV->position()));
                    t_Leptons_dz.push_back(posEl.gsfTrack()->dz(PV->position()));


                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posEl, genColl,11);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());
                    }
                    negEl = MyElectrons[1];
                    isNegEl = true;
                    t_Leptons.push_back(negEl.p4());

                    t_Leptons_pt.push_back(negEl.pt());
                    t_Leptons_eta.push_back(negEl.eta());
                    t_Leptons_phi.push_back(negEl.phi());
                    t_Leptons_e.push_back(negEl.energy());
                    t_Leptons_charge.push_back(negEl.charge());
                    t_Leptons_type.push_back("electron");
                    t_Leptons_etaSC.push_back(negEl.superCluster()->eta());
                    t_Leptons_d0.push_back(negEl.gsfTrack()->dxy(PV->position()));
                    t_Leptons_dz.push_back(negEl.gsfTrack()->dz(PV->position()));

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negEl, genColl,11);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());



                    }

                }
                else
                {

                    posEl = MyElectrons[1];
                    isPosEl = true;
                    t_Leptons.push_back(posEl.p4());

                    t_Leptons_pt.push_back(posEl.pt());
                    t_Leptons_eta.push_back(posEl.eta());
                    t_Leptons_phi.push_back(posEl.phi());
                    t_Leptons_e.push_back(posEl.energy());
                    t_Leptons_charge.push_back(posEl.charge());
                    t_Leptons_type.push_back("electron");
                    t_Leptons_etaSC.push_back(posEl.superCluster()->eta());
                    t_Leptons_d0.push_back(posEl.gsfTrack()->dxy(PV->position()));
                    t_Leptons_dz.push_back(posEl.gsfTrack()->dz(PV->position()));


                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posEl, genColl,11);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());
                    }
                    negEl = MyElectrons[0];
                    isNegEl = true;
                    t_Leptons.push_back(negEl.p4());

                    t_Leptons_pt.push_back(negEl.pt());
                    t_Leptons_eta.push_back(negEl.eta());
                    t_Leptons_phi.push_back(negEl.phi());
                    t_Leptons_e.push_back(negEl.energy());
                    t_Leptons_charge.push_back(negEl.charge());
                    t_Leptons_type.push_back("electron");
                    t_Leptons_etaSC.push_back(negEl.superCluster()->eta());
                    t_Leptons_d0.push_back(negEl.gsfTrack()->dxy(PV->position()));
                    t_Leptons_dz.push_back(negEl.gsfTrack()->dz(PV->position()));

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negEl, genColl,11);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());



                    }


                }

            }
        }
        else if (my_emu)
        {
            if(MyElectrons[0].charge() != MyMuons[0].charge() )
            {
                isElMu = true ;
                isMuEl = true;
                if(MyElectrons[0].charge() > 0)
                {
                    posEl = MyElectrons[0];
                    isPosEl = true;
                    t_Leptons.push_back(posEl.p4());

                    t_Leptons_pt.push_back(posEl.pt());
                    t_Leptons_eta.push_back(posEl.eta());
                    t_Leptons_phi.push_back(posEl.phi());
                    t_Leptons_e.push_back(posEl.energy());
                    t_Leptons_charge.push_back(posEl.charge());
                    t_Leptons_type.push_back("electron");
                    t_Leptons_etaSC.push_back(posEl.superCluster()->eta());
                    t_Leptons_d0.push_back(posEl.gsfTrack()->d0());
                    t_Leptons_dz.push_back(posEl.gsfTrack()->dz());


                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posEl, genColl,11);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());
                    }
                    negMu = MyMuons[0];
                    isNegMu = true;
                    t_Leptons.push_back(negMu.p4());

                    t_Leptons_pt.push_back(negMu.pt());
                    t_Leptons_eta.push_back(negMu.eta());
                    t_Leptons_phi.push_back(negMu.phi());
                    t_Leptons_e.push_back(negMu.energy());
                    t_Leptons_charge.push_back(negMu.charge());
                    t_Leptons_type.push_back("muon");

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negMu, genColl,13);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());



                    }


                }
                else
                {
                    posMu = MyMuons[0];
                    isPosMu = true;
                    t_Leptons.push_back(posMu.p4());

                    t_Leptons_pt.push_back(posMu.pt());
                    t_Leptons_eta.push_back(posMu.eta());
                    t_Leptons_phi.push_back(posMu.phi());
                    t_Leptons_e.push_back(posMu.energy());
                    t_Leptons_charge.push_back(posMu.charge());
                    t_Leptons_type.push_back("muon");

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posMu, genColl,13);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());

                    }
                    negEl = MyElectrons[0];
                    isNegEl = true;
                    t_Leptons.push_back(negEl.p4());

                    t_Leptons_pt.push_back(negEl.pt());
                    t_Leptons_eta.push_back(negEl.eta());
                    t_Leptons_phi.push_back(negEl.phi());
                    t_Leptons_e.push_back(negEl.energy());
                    t_Leptons_charge.push_back(negEl.charge());
                    t_Leptons_type.push_back("electron");
                    t_Leptons_etaSC.push_back(negEl.superCluster()->eta());
                    t_Leptons_d0.push_back(negEl.gsfTrack()->d0());
                    t_Leptons_dz.push_back(negEl.gsfTrack()->dz());

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negEl, genColl,11);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());



                    }



                }
            }
        }
        else if (my_mumu)
        {
            if(MyMuons[0].charge() != MyMuons[1].charge() )
            {
                isDiMuon = true;
                if(MyMuons[0].charge() > 0)
                {
                    posMu = MyMuons[0];
                    isPosMu = true;
                    t_Leptons.push_back(posMu.p4());

                    t_Leptons_pt.push_back(posMu.pt());
                    t_Leptons_eta.push_back(posMu.eta());
                    t_Leptons_phi.push_back(posMu.phi());
                    t_Leptons_e.push_back(posMu.energy());
                    t_Leptons_charge.push_back(posMu.charge());
                    t_Leptons_type.push_back("muon");

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posMu, genColl,13);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());

                    }
                    negMu = MyMuons[1];
                    isNegMu = true;
                    t_Leptons.push_back(negMu.p4());

                    t_Leptons_pt.push_back(negMu.pt());
                    t_Leptons_eta.push_back(negMu.eta());
                    t_Leptons_phi.push_back(negMu.phi());
                    t_Leptons_e.push_back(negMu.energy());
                    t_Leptons_charge.push_back(negMu.charge());
                    t_Leptons_type.push_back("muon");

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negMu, genColl,13);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());



                    }

                }
                else
                {
                    posMu = MyMuons[1];
                    isPosMu = true;
                    t_Leptons.push_back(posMu.p4());

                    t_Leptons_pt.push_back(posMu.pt());
                    t_Leptons_eta.push_back(posMu.eta());
                    t_Leptons_phi.push_back(posMu.phi());
                    t_Leptons_e.push_back(posMu.energy());
                    t_Leptons_charge.push_back(posMu.charge());
                    t_Leptons_type.push_back("muon");

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(posMu, genColl,13);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());

                    }
                    negMu = MyMuons[0];
                    isNegMu = true;
                    t_Leptons.push_back(negMu.p4());

                    t_Leptons_pt.push_back(negMu.pt());
                    t_Leptons_eta.push_back(negMu.eta());
                    t_Leptons_phi.push_back(negMu.phi());
                    t_Leptons_e.push_back(negMu.energy());
                    t_Leptons_charge.push_back(negMu.charge());
                    t_Leptons_type.push_back("muon");

                    if(!isData)
                    {
                        const pat::PackedGenParticle matched = mse::getMatchedGenParticle(negMu, genColl,13);

                        const reco::GenParticle daughter = mse::getMotherPacked(matched);
                        const reco::GenParticle mom = mse::getMother(daughter);

                        t_gen_Leptons_pt.push_back(matched.pt());
                        t_gen_Leptons_eta.push_back(matched.eta());
                        t_gen_Leptons_phi.push_back(matched.phi());
                        t_gen_Leptons_e.push_back(matched.energy());
                        t_gen_Leptons_id.push_back(matched.pdgId());
                        t_gen_Leptons_charge.push_back(matched.charge());

                        t_gen_daughters_pt.push_back(daughter.pt());
                        t_gen_daughters_eta.push_back(daughter.eta());
                        t_gen_daughters_phi.push_back(daughter.phi());
                        t_gen_daughters_e.push_back(daughter.energy());
                        t_gen_daughters_id.push_back(daughter.pdgId());

                        t_gen_mothers_pt.push_back(mom.pt());
                        t_gen_mothers_eta.push_back(mom.eta());
                        t_gen_mothers_phi.push_back(mom.phi());
                        t_gen_mothers_e.push_back(mom.energy());
                        t_gen_mothers_id.push_back(mom.pdgId());



                    }
                }
            }
        }
    }
    }
    //cout << my_ee + my_emu + my_mumu << " what is the number ! " <<endl;
    t_Muons_size = 0;
    t_Electrons_size = 0;
    t_Leptons_size = 0;
    t_Leptons_size = MyElectrons.size() + MyMuons.size();
    t_Electrons_size = MyElectrons.size();
    t_Muons_size = MyMuons.size();
    //   if(DiMu && !isDiMuon) return;
    //   if(DiEl && !isDiElectron) return;
    //   if(ElMu && !(isElMu || isMuEl) ) return;


    //After Lepton Selection
    if(isPosMu && isNegMu) isDiMuon = true;
    if(isPosEl && isNegEl) isDiElectron = true;
    if( isPosEl && isNegMu)  isElMu = true;
    if( isNegEl && isPosMu ) isMuEl = true;
    if(isDiMuon || isDiElectron || isElMu ) isDiLeptonic = true;
    if(!isDiLeptonic)return;
    if(DiMu && !isDiMuon) return;
    if(DiEl && !isDiElectron) return;
    if(ElMu && !(isElMu || isMuEl) ) return;
    //    //++n_afterDiLepton;


    ///////////////////JETS//////////////////////////////
    /// \brief bjets
    ///tightLepVeto
    vector<pat::Jet> bjets;
    vector<pat::Jet> njets;
    t_Jets_e.clear();
    t_Jets_eta.clear();
    t_Jets_phi.clear();
    t_Jets_pt.clear();
    t_Jets_hadflav.clear();
    t_Jets_pt_resolution.clear();
    t_Jets_sf_up.clear();
    t_Jets_sf_down.clear();
    t_Jets_sf_nominal.clear();
    t_gen_Jets_e.clear();
    t_gen_Jets_eta.clear();
    t_gen_Jets_phi.clear();
    t_gen_Jets_pt.clear();
    t_gen_Jets_hadflav.clear();
    t_gen_Jets_prtnflav.clear();
    t_Jets_LBdiscriminator.clear();
    t_Jets_MBdiscriminator.clear();
    t_Jets_TBdiscriminator.clear();
    t_Jets_dis.clear();
    t_bJets_e.clear();
    t_bJets_eta.clear();
    t_bJets_phi.clear();
    t_bJets_pt.clear();
    //    cout << "here I jet 1!" <<endl;

    for(pat::JetCollection::const_iterator jet_it=jets->begin(); jet_it != jets->end(); ++jet_it)
    {
        if( !( jet_it->pt() > 30 )) continue;
        if( !( fabs(jet_it->eta()) < 2.4 )) continue;

        //        if( !(jet_it->neutralHadronEnergyFraction() < 0.90 && jet_it->neutralEmEnergyFraction() < 0.90 && (jet_it->chargedMultiplicity() + jet_it->neutralMultiplicity())> 1.
        //                && jet_it->muonEnergyFraction() < 0.8  && jet_it->chargedHadronEnergyFraction() > 0. && jet_it->chargedEmEnergyFraction() < 0.90 && jet_it->chargedMultiplicity() > 0.)) continue;
        auto NHF  = jet_it->neutralHadronEnergyFraction();
        auto NEMF = jet_it->neutralEmEnergyFraction();
        auto CHF  = jet_it->chargedHadronEnergyFraction();
        auto MUF  = jet_it->muonEnergyFraction();
        auto CEMF = jet_it->chargedEmEnergyFraction();
        auto NumConst = jet_it->chargedMultiplicity()+jet_it->neutralMultiplicity();
        auto NumNeutralParticles =jet_it->neutralMultiplicity();
        auto CHM      = jet_it->chargedMultiplicity();
        bool looseJetID = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && (abs(jet_it->eta())<=2.4 && CHF>=0 && CHM>0 && CEMF<=0.99)) ;
        if(!looseJetID)continue;
        if(isDiMuon)
        {
            LorentzVector jetp4(jet_it->p4());
            LorentzVector mupp4(posMu.p4());
            LorentzVector mump4(negMu.p4());
            double dr1 = ROOT::Math::VectorUtil::DeltaR(jetp4,mupp4);
            double dr2 = ROOT::Math::VectorUtil::DeltaR(jetp4,mump4);
            if(dr1 < 0.4 || dr2 < 0.4 ) continue;
        }

        if(isDiElectron)
        {
            LorentzVector jetp4(jet_it->p4());
            LorentzVector mupp4(posEl.p4());
            LorentzVector mump4(negEl.p4());
            double dr1 = ROOT::Math::VectorUtil::DeltaR(jetp4,mupp4);
            double dr2 = ROOT::Math::VectorUtil::DeltaR(jetp4,mump4);
            if(dr1 < 0.4 || dr2 < 0.4 ) continue;
        }
        if(isElMu)
        {
            LorentzVector jetp4(jet_it->p4());
            LorentzVector mupp4(posEl.p4());
            LorentzVector mump4(negMu.p4());
            double dr1 = ROOT::Math::VectorUtil::DeltaR(jetp4,mupp4);
            double dr2 = ROOT::Math::VectorUtil::DeltaR(jetp4,mump4);
            if(dr1 < 0.4 || dr2 < 0.4 ) continue;
        }
        if(isMuEl)
        {
            LorentzVector jetp4(jet_it->p4());
            LorentzVector mupp4(posMu.p4());
            LorentzVector mump4(negEl.p4());
            double dr1 = ROOT::Math::VectorUtil::DeltaR(jetp4,mupp4);
            double dr2 = ROOT::Math::VectorUtil::DeltaR(jetp4,mump4);
            if(dr1 < 0.4 || dr2 < 0.4 ) continue;
        }
        njets.push_back(*jet_it);
        t_Jets.push_back(jet_it->p4());
        t_Jets_pt.push_back(jet_it->pt());
        t_Jets_eta.push_back(jet_it->eta());
        t_Jets_phi.push_back(jet_it->phi());
        t_Jets_e.push_back(jet_it->energy());
        t_Jets_hadflav.push_back(jet_it->hadronFlavour());
        t_Jets_dis.push_back(jet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        t_Jets_LBdiscriminator.push_back(jet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.5426  );
        t_Jets_MBdiscriminator.push_back(jet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.8484  );
        t_Jets_TBdiscriminator.push_back(jet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.9535  );
        double jet_resolution = ptResol->getResolution({{JME::Binning::JetPt, jet_it->pt()}, {JME::Binning::JetEta, jet_it->eta()}, {JME::Binning::Rho, rho}});
        t_Jets_pt_resolution.push_back(jet_resolution);
        double jet_phi_resolution = phiResol->getResolution({{JME::Binning::JetPt, jet_it->pt()}, {JME::Binning::JetEta, jet_it->eta()}, {JME::Binning::Rho, rho}});
        t_Jets_phi_resolution.push_back(jet_phi_resolution);

        double jer_sf_nominal = jerSF->getScaleFactor({{JME::Binning::JetEta, jet_it->eta()}}, Variation::NOMINAL);
        double jer_sf_up = jerSF->getScaleFactor({{JME::Binning::JetEta, jet_it->eta()}}, Variation::UP);
        double jer_sf_down = jerSF->getScaleFactor({{JME::Binning::JetEta, jet_it->eta()}}, Variation::DOWN);
        t_Jets_sf_nominal.push_back(jer_sf_nominal);
        t_Jets_sf_up.push_back(jer_sf_up);
        t_Jets_sf_down.push_back(jer_sf_down);

        if(!isData)
        {
            LorentzVector jetp4(jet_it->p4());
            //                // Try to find a gen jet matching
            //                // dR < m_dR_max
            //                // dPt < m_dPt_max_factor * resolution

            double min_dR = std::numeric_limits<double>::infinity();

            reco::GenJet gen_jet;
            double m_dR_max = 0.2;
            double m_dPt_max_factor = 3;
            for (reco::GenJetCollection::const_iterator it_genJet = genJets->begin(); it_genJet != genJets->end();++it_genJet) {
                LorentzVector genp4(it_genJet->p4());
                double dR = ROOT::Math::VectorUtil::DeltaR(genp4, jetp4);

                if (dR > min_dR)
                    continue;

                if (dR < m_dR_max) {
                    double dPt = std::abs(it_genJet->pt() - jet_it->pt());
                    if (dPt > m_dPt_max_factor * jet_resolution)
                        continue;

                    min_dR = dR;
                    gen_jet = *it_genJet;
                }
            }

            //            const reco::GenJet gen_jet = mse::match(*jet_it,jet_resolution,*genJets);
            t_gen_Jets_pt.push_back(gen_jet.pt());
            t_gen_Jets_eta.push_back(gen_jet.eta());
            t_gen_Jets_phi.push_back(gen_jet.phi());
            t_gen_Jets_e.push_back(gen_jet.energy());

        }
        if( !(jet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.5426  )) continue;          //0.5426 is loose working point for btaging new value is 0.605
        bjets.push_back(*jet_it);
        t_bJets.push_back(jet_it->p4());
        t_bJets_pt.push_back(jet_it->pt());
        t_bJets_eta.push_back(jet_it->eta());
        t_bJets_phi.push_back(jet_it->phi());
        t_bJets_e.push_back(jet_it->energy());

    }
    t_Met_pt.clear();
    t_Met_px.clear();
    t_Met_py.clear();
    t_Met_sum.clear();
    const pat::MET &met = mets->front();
    t_Met_pt.push_back(met.pt());
    t_Met_px.push_back(met.px());
    t_Met_py.push_back(met.py());
    t_Met_sum.push_back(met.sumEt());
    //    cout << "here I met 1!" <<endl;


    //////////////////////Lepton SFs//////////////////////////////////////

    if(!isData)
    {
        int count = 0;
        for(vector<pat::Jet>::const_iterator bjet_it = njets.begin(); bjet_it != njets.end(); ++bjet_it)
        {

            double jet_scalefactor    = reader.eval_auto_bounds(
                        "central",
                        BTagEntry::FLAV_B,
                        bjet_it->eta(),
                        bjet_it->pt()
                        );
            //cout <<"btag: "<< jet_scalefactor << endl;
            w_bjets.push_back(jet_scalefactor);
            if(count <= 2)
            {
                w_bjet = w_bjet*jet_scalefactor;
                theWeight = theWeight*jet_scalefactor;
            }
            ++count;


        }

        TH2F* h2D_egammaSF = (TH2F*) f_egammaSF->Get("EGamma_SF2D");
        TH2F* h2D_egammaTkSF = (TH2F*) f_egammaTkSF->Get("EGamma_SF2D");

        if(isPosEl)
        {
            Int_t binx = h2D_egammaSF->GetXaxis()->FindBin(posEl.superCluster()->eta());
            Int_t biny = h2D_egammaSF->GetYaxis()->FindBin(posEl.pt());
            Int_t binx2 = h2D_egammaTkSF->GetXaxis()->FindBin(posEl.superCluster()->eta());
            Int_t biny2 = 1;
            //            (posEl.pt()<26) ? ( biny2 = h2D_egammaTkSF->GetYaxis()->FindBin(50)) : ( biny2 = h2D_egammaTkSF->GetYaxis()->FindBin(posEl.pt()) );

            //cout <<"pt: "<< posEl.pt() << " eta: " << posEl.superCluster()->eta() <<" SF: "<< h2D_egammaSF->GetBinContent(binx,biny)<<" & "<< h2D_egammaTkSF->GetBinContent(binx2,biny2) <<endl;
            SF_posEl = h2D_egammaSF->GetBinContent(binx,biny)*h2D_egammaTkSF->GetBinContent(binx2,biny2);
            w_posEl = w_posEl*SF_posEl;
            theWeight = theWeight*SF_posEl;

        }
        if(isNegEl)
        {
            Int_t binx = h2D_egammaSF->GetXaxis()->FindBin(negEl.superCluster()->eta());
            Int_t biny = h2D_egammaSF->GetYaxis()->FindBin(negEl.pt());
            Int_t binx2 = h2D_egammaTkSF->GetXaxis()->FindBin(negEl.superCluster()->eta());
            Int_t biny2 = 1;
            //            (negEl.pt() < 26) ? ( biny2 = h2D_egammaTkSF->GetYaxis()->FindBin(50)) : ( biny2 = h2D_egammaTkSF->GetYaxis()->FindBin(negEl.pt()) );
            SF_negEl = h2D_egammaSF->GetBinContent(binx,biny)*h2D_egammaTkSF->GetBinContent(binx2,biny2);
            w_negEl = w_negEl*SF_negEl;
            theWeight = theWeight*SF_negEl;

        }

        TH2F* h2D_muonIDSF = (TH2F*) f_muonIDSF->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
        TH2F* h2D_muonISOSF = (TH2F*) f_muonISOSF->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
        TGraph* g_muonTkSF = (TGraph*) f_muonTkSF->Get("ratio_eff_aeta_dr030e030_corr");

        if(isPosMu)
        {
            Int_t binx = h2D_muonIDSF->GetXaxis()->FindBin(fabs(posMu.eta()));
            Int_t biny = h2D_muonIDSF->GetYaxis()->FindBin(posMu.pt());
            Int_t binx2 = h2D_muonISOSF->GetXaxis()->FindBin(fabs(posMu.eta()));
            Int_t biny2 = h2D_muonISOSF->GetYaxis()->FindBin(fabs(posMu.pt()));
            (posMu.pt() > 119) ? biny = h2D_muonIDSF->GetYaxis()->FindBin(118):biny = h2D_muonIDSF->GetYaxis()->FindBin(posMu.pt());
            (posMu.pt() > 119) ? biny2 = h2D_muonISOSF->GetYaxis()->FindBin(118):biny2 = h2D_muonISOSF->GetYaxis()->FindBin(posMu.pt());
            Double_t TkSF = g_muonTkSF->Eval(fabs(posMu.eta()));
            //cout <<"Mupt: "<< posMu.pt() << " eta: " << fabs(posMu.eta()) <<" SF: "<< h2D_muonIDSF->GetBinContent(binx,biny)<<"& "<< h2D_muonISOSF->GetBinContent(binx2,biny2) <<" & "<< TkSF <<endl;
            SF_posMu = h2D_muonIDSF->GetBinContent(binx,biny)*h2D_muonISOSF->GetBinContent(binx2,biny2)*TkSF;
            w_posMu = w_posMu*SF_posMu;
            theWeight = theWeight*SF_posMu;

        }
        if(isNegMu)
        {
            Int_t binx = h2D_muonIDSF->GetXaxis()->FindBin(fabs(negMu.eta()));
            Int_t biny = h2D_muonIDSF->GetYaxis()->FindBin(negMu.pt());
            Int_t binx2 = h2D_muonISOSF->GetXaxis()->FindBin(fabs(negMu.eta()));
            Int_t biny2 = h2D_muonISOSF->GetYaxis()->FindBin(fabs(negMu.pt()));
            (negMu.pt() > 119) ? biny = h2D_muonIDSF->GetYaxis()->FindBin(118): biny = h2D_muonIDSF->GetYaxis()->FindBin(negMu.pt());
            (negMu.pt() > 119) ? biny2 = h2D_muonISOSF->GetYaxis()->FindBin(118): biny2 = h2D_muonISOSF->GetYaxis()->FindBin(negMu.pt());
            Double_t TkSF = g_muonTkSF->Eval(fabs(negMu.eta()));
            //cout <<"Mupt: "<< negMu.pt() << " eta: " << fabs(negMu.eta()) <<" SF: "<< h2D_muonIDSF->GetBinContent(binx,biny)<<"& "<< h2D_muonISOSF->GetBinContent(binx2,biny2)<< " & "<< TkSF <<endl;
            SF_negMu = h2D_muonIDSF->GetBinContent(binx,biny)*h2D_muonISOSF->GetBinContent(binx2,biny2)*TkSF;
            w_negMu = w_negMu*SF_negMu;
            theWeight = theWeight*SF_negMu;


        }

        if(DiEl)
        {
            TH2F* h2D_ee_sf = (TH2F*) f_ee_sf->Get("scalefactor_eta2d_with_syst");

            Int_t binx = h2D_ee_sf->GetXaxis()->FindBin(fabs(posEl.eta()));
            Int_t biny = h2D_ee_sf->GetYaxis()->FindBin(fabs(negEl.eta()));

            double d_ee_sf = h2D_ee_sf->GetBinContent(binx,biny);
            //            cout << posEl.eta() << " ee " << negEl.eta() << " sf " << d_ee_sf << endl;
            w_ee = w_ee*d_ee_sf;
            theWeight = theWeight*d_ee_sf;
        }
        if(ElMu)
        {
            TH2F* h2D_em_sf = (TH2F*) f_em_sf->Get("scalefactor_eta2d_with_syst");
            Int_t binx = 0;
            Int_t biny = 0;
            if(isElMu)
            {
                binx = h2D_em_sf->GetXaxis()->FindBin(fabs(posEl.eta()));
                biny = h2D_em_sf->GetYaxis()->FindBin(fabs(negMu.eta()));
            }
            if(isMuEl)
            {
                binx = h2D_em_sf->GetXaxis()->FindBin(fabs(posMu.eta()));
                biny = h2D_em_sf->GetYaxis()->FindBin(fabs(negEl.eta()));
            }
            double d_em_sf = h2D_em_sf->GetBinContent(binx,biny);
            //            cout << d_em_sf << " em " << endl;
            w_em = w_em*d_em_sf;
            theWeight = theWeight*d_em_sf;
        }
        if(DiMu)
        {
            TH2F* h2D_mm_sf = (TH2F*) f_mm_sf->Get("scalefactor_eta2d_with_syst");

            Int_t binx = h2D_mm_sf->GetXaxis()->FindBin(fabs(posMu.eta()));
            Int_t biny = h2D_mm_sf->GetYaxis()->FindBin(fabs(negMu.eta()));

            double d_mm_sf = h2D_mm_sf->GetBinContent(binx,biny);
            // cout << posMu.eta() << " ee " << negMu.eta() << " sf " << d_mm_sf << endl;
            w_mm = w_mm * d_mm_sf;
            theWeight = theWeight*d_mm_sf;
        }
        //            if(genEvent->isFullLeptonic(true))
        //            {
        //
        //            cout << "why not here!" << endl;
        //            if(genEvent->leptonicDecayTop(true) != NULL && genEvent->hadronicDecayTop(true) != NULL)
        //            {
        //
        //            cout << genEvent->leptonicDecayTop(true)->pt()<< endl;
        //            cout <<  genEvent->hadronicDecayTop(true)->pt() << endl;
        //            double topPtLepTrue = genEvent->leptonicDecayTop(true)->pt();
        //            double topPtHadTrue = genEvent->hadronicDecayTop(true)->pt();
        //            cout << "So far we are here" << endl;
        //            double w = sqrt(SF(topPtLepTrue)*SF(topPtHadTrue));
        //            w_top = w_top*w;
        //            theWeight=theWeight*w;
        //            }
        //            }

    }

    nGoodVtxs = 0;
    for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx)
    {
        bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
        if ( !(isFake) && vtx->ndof() >= 4. && fabs(vtx->position().Rho()) <= 2.0 && fabs(vtx->position().Z()) <= 24.) nGoodVtxs++;
    }
    t_num_PV = nGoodVtxs;
    h_NPV->Fill(nGoodVtxs,1);
   
    //cout << "IfilltheTree"<<endl;
    
    
    
    if(bjets.size() >= 2)
    {
           //////////////////////TOP RECO/////////////////
    //    TtDilepEvtSolution asol;

    //    asol.setGenEvt(genEvent);
    t_top_pt.clear();
    t_top_eta.clear();
    t_top_phi.clear();
    t_top_e.clear();
    
    t_antitop_pt.clear();
    t_antitop_eta.clear();
    t_antitop_phi.clear();
    t_antitop_e.clear();
    
    t_W_pt.clear();
    t_W_eta.clear();
    t_W_phi.clear();
    t_W_e.clear();
    
    t_antiW_pt.clear();
    t_antiW_eta.clear();
    t_antiW_phi.clear();
    t_antiW_e.clear();
    
    t_nu_pt.clear();
    t_nu_eta.clear();
    t_nu_phi.clear();
    t_nu_e.clear();
    
    t_antinu_pt.clear();
    t_antinu_eta.clear();
    t_antinu_phi.clear();
    t_antinu_e.clear();
    t_cos.clear();
    
    
    ////////DiMuon:
    if(isDiMuon)
    {
        //cout << "IM HERE!MuMu" << endl;
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

        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 )
        {
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

            t_antinu.push_back(LorentzVector(nuBar.Pt(),nuBar.Eta(),nuBar.Phi(),nuBar.M()));
            t_antinu_pt.push_back(nuBar.Pt());
            t_antinu_eta.push_back(nuBar.Eta());
            t_antinu_phi.push_back(nuBar.Phi());
            t_antinu_e.push_back(nuBar.E());
            t_nu.push_back(LorentzVector(nu.Pt(),nu.Eta(),nu.Phi(),nu.M()));
            t_nu_pt.push_back(nu.Pt());
            t_nu_eta.push_back(nu.Eta());
            t_nu_phi.push_back(nu.Phi());
            t_nu_e.push_back(nu.E());
            t_antiW.push_back(LorentzVector(W2.Pt(),W2.Eta(),W2.Phi(),W2.M()));
            t_antiW_pt.push_back(W2.Pt());
            t_antiW_eta.push_back(W2.Eta());
            t_antiW_phi.push_back(W2.Phi());
            t_antiW_e.push_back(W2.E());
            t_W.push_back(LorentzVector(W1.Pt(),W1.Eta(),W1.Phi(),W1.M()));
            t_W_pt.push_back(W1.Pt());
            t_W_eta.push_back(W1.Eta());
            t_W_phi.push_back(W1.Phi());
            t_W_e.push_back(W1.E());
            t_antitop.push_back(LorentzVector(t2.Pt(),t2.Eta(),t2.Phi(),t2.M()));
            t_antitop_pt.push_back(t2.Pt());
            t_antitop_eta.push_back(t2.Eta());
            t_antitop_phi.push_back(t2.Phi());
            t_antitop_e.push_back(t2.E());
            t_top.push_back(LorentzVector(t1.Pt(),t1.Eta(),t1.Phi(),t1.M()));
            t_top_pt.push_back(t1.Pt());
            t_top_eta.push_back(t1.Eta());
            t_top_phi.push_back(t1.Phi());
            t_top_e.push_back(t1.E());

            //++n_afterTop;
            h_Nevents_top->Fill(1,theWeight);
            h_NeventsMuMu_top->Fill(1,theWeight);
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


            t_cos.push_back(TMath::Cos(theta1));
            t_cos.push_back(TMath::Cos(theta2));


            if(!isData)
            {
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
                //                if(matchedPosMu.pt() != 0 )  cout << matchedPosMu.pt() << "muon match mother"<< posWMuMu.pdgId()<<"topMuMu Mother pdgId"<< topMuMu.pdgId() << endl;
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
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_RecoMinusTruthMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }

        h_Nevents_ATS->Fill(1,theWeight);
        h_NeventsMuMu_ATS->Fill(1,theWeight);
        h_ATS_ptLepMuMu->Fill(PosLep.Pt(),theWeight);
        h_ATS_ptLepMuMu->Fill(NegLep.Pt(),theWeight);
        h_ATS_etaLepMuMu->Fill(PosLep.Eta(),theWeight);
        h_ATS_etaLepMuMu->Fill(NegLep.Eta(),theWeight);
        h_ATS_mLepMuMu->Fill(dilepton.M(),theWeight);
        h_ATS_METMuMu->Fill(met.pt(),theWeight);
        h_ATS_NBJetsMuMu->Fill(bjets.size(),theWeight);
        h_ATS_NJetsMuMu->Fill(njets.size(),theWeight);


        h_NBJetsMuMu->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);




    }


    //////DiElectron:
    if(isDiElectron )
    {
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
        if(nuSol.neutrino.pt()> 20 && nuSol.neutrinoBar.pt() > 20 )
        {
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
            //++n_afterTop;
            h_Nevents_top->Fill(1,theWeight);
            h_NeventsElEl_top->Fill(1,theWeight);
            RecoMTT.push_back(ttbar.M());

             t_antinu.push_back(LorentzVector(nuBar.Pt(),nuBar.Eta(),nuBar.Phi(),nuBar.M()));
            t_antinu_pt.push_back(nuBar.Pt());
            t_antinu_eta.push_back(nuBar.Eta());
            t_antinu_phi.push_back(nuBar.Phi());
            t_antinu_e.push_back(nuBar.E());
            t_nu.push_back(LorentzVector(nu.Pt(),nu.Eta(),nu.Phi(),nu.M()));
            t_nu_pt.push_back(nu.Pt());
            t_nu_eta.push_back(nu.Eta());
            t_nu_phi.push_back(nu.Phi());
            t_nu_e.push_back(nu.E());
            t_antiW.push_back(LorentzVector(W2.Pt(),W2.Eta(),W2.Phi(),W2.M()));
            t_antiW_pt.push_back(W2.Pt());
            t_antiW_eta.push_back(W2.Eta());
            t_antiW_phi.push_back(W2.Phi());
            t_antiW_e.push_back(W2.E());
            t_W.push_back(LorentzVector(W1.Pt(),W1.Eta(),W1.Phi(),W1.M()));
            t_W_pt.push_back(W1.Pt());
            t_W_eta.push_back(W1.Eta());
            t_W_phi.push_back(W1.Phi());
            t_W_e.push_back(W1.E());
            t_antitop.push_back(LorentzVector(t2.Pt(),t2.Eta(),t2.Phi(),t2.M()));
            t_antitop_pt.push_back(t2.Pt());
            t_antitop_eta.push_back(t2.Eta());
            t_antitop_phi.push_back(t2.Phi());
            t_antitop_e.push_back(t2.E());
            t_top.push_back(LorentzVector(t1.Pt(),t1.Eta(),t1.Phi(),t1.M()));
            t_top_pt.push_back(t1.Pt());
            t_top_eta.push_back(t1.Eta());
            t_top_phi.push_back(t1.Phi());
            t_top_e.push_back(t1.E());



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

            t_cos.push_back(TMath::Cos(theta1));
            t_cos.push_back(TMath::Cos(theta2));

            if(!isData)
            {
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
                    h_truthRecoElEl->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoCos->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoMTT->Fill(ttbar.M(),genttbar.M());
                    h_truthRecoElEl->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_truthRecoCos->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_RecoMinusTruthMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }
        //cout << "Am I here Now" <<endl;
        h_Nevents_ATS->Fill(1,theWeight);
        h_NeventsElEl_ATS->Fill(1,theWeight);
        h_ATS_ptLepElEl->Fill(PosLep.Pt(),theWeight);
        h_ATS_ptLepElEl->Fill(NegLep.Pt(),theWeight);
        h_ATS_etaLepElEl->Fill(PosLep.Eta(),theWeight);
        h_ATS_etaLepElEl->Fill(NegLep.Eta(),theWeight);
        h_ATS_mLepElEl->Fill(dilepton.M(),theWeight);
        h_ATS_METElEl->Fill(met.pt(),theWeight);
        h_ATS_NBJetsElEl->Fill(bjets.size(),theWeight);
        h_ATS_NJetsElEl->Fill(njets.size(),theWeight);


        // cout << "what 5" << endl;
        h_NBJetsElEl->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);
        // cout << "what 6" << endl;
    }
    //////ElectronMuon:
    if(isElMu )
    {
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
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 )
        {
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

             t_antinu.push_back(LorentzVector(nuBar.Pt(),nuBar.Eta(),nuBar.Phi(),nuBar.M()));
            t_antinu_pt.push_back(nuBar.Pt());
            t_antinu_eta.push_back(nuBar.Eta());
            t_antinu_phi.push_back(nuBar.Phi());
            t_antinu_e.push_back(nuBar.E());
            t_nu.push_back(LorentzVector(nu.Pt(),nu.Eta(),nu.Phi(),nu.M()));
            t_nu_pt.push_back(nu.Pt());
            t_nu_eta.push_back(nu.Eta());
            t_nu_phi.push_back(nu.Phi());
            t_nu_e.push_back(nu.E());
            t_antiW.push_back(LorentzVector(W2.Pt(),W2.Eta(),W2.Phi(),W2.M()));
            t_antiW_pt.push_back(W2.Pt());
            t_antiW_eta.push_back(W2.Eta());
            t_antiW_phi.push_back(W2.Phi());
            t_antiW_e.push_back(W2.E());
            t_W.push_back(LorentzVector(W1.Pt(),W1.Eta(),W1.Phi(),W1.M()));
            t_W_pt.push_back(W1.Pt());
            t_W_eta.push_back(W1.Eta());
            t_W_phi.push_back(W1.Phi());
            t_W_e.push_back(W1.E());
            t_antitop.push_back(LorentzVector(t2.Pt(),t2.Eta(),t2.Phi(),t2.M()));
            t_antitop_pt.push_back(t2.Pt());
            t_antitop_eta.push_back(t2.Eta());
            t_antitop_phi.push_back(t2.Phi());
            t_antitop_e.push_back(t2.E());
            t_top.push_back(LorentzVector(t1.Pt(),t1.Eta(),t1.Phi(),t1.M()));
            t_top_pt.push_back(t1.Pt());
            t_top_eta.push_back(t1.Eta());
            t_top_phi.push_back(t1.Phi());
            t_top_e.push_back(t1.E());

            //++n_afterTop;
            h_Nevents_top->Fill(1,theWeight);
            h_NeventsElMu_top->Fill(1,theWeight);
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


            t_cos.push_back(TMath::Cos(theta1));
            t_cos.push_back(TMath::Cos(theta2));

            if(!isData)
            {
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
                    h_truthRecoElMu->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoCos->Fill(TMath::Cos(theta1),TMath::Cos(gentheta1));
                    h_truthRecoMTT->Fill(ttbar.M(),genttbar.M());
                    h_truthRecoElMu->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_truthRecoCos->Fill(TMath::Cos(theta2),TMath::Cos(gentheta2));
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_RecoMinusTruthMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }
        //cout << "Am I here Now 2" <<endl;
        h_Nevents_ATS->Fill(1,theWeight);
        h_NeventsElMu_ATS->Fill(1,theWeight);
        h_ATS_ptLepElMu->Fill(PosLep.Pt(),theWeight);
        h_ATS_ptLepElMu->Fill(NegLep.Pt(),theWeight);
        h_ATS_etaLepElMu->Fill(PosLep.Eta(),theWeight);
        h_ATS_etaLepElMu->Fill(NegLep.Eta(),theWeight);
        h_ATS_mLepElMu->Fill(dilepton.M(),theWeight);
        h_ATS_METElMu->Fill(met.pt(),theWeight);
        h_ATS_NBJetsElMu->Fill(bjets.size(),theWeight);
        h_ATS_NJetsElMu->Fill(njets.size(),theWeight);

        h_NBJetsElMu->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);
    }

    if(isMuEl)
    {
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
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 )
        {
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

            //++n_afterTop;
            h_Nevents_top->Fill(1,theWeight);
            h_NeventsElMu_top->Fill(1,theWeight);
            RecoMTT.push_back(ttbar.M());
            t_antinu.push_back(LorentzVector(nuBar.Pt(),nuBar.Eta(),nuBar.Phi(),nuBar.M()));
            t_antinu_pt.push_back(nuBar.Pt());
            t_antinu_eta.push_back(nuBar.Eta());
            t_antinu_phi.push_back(nuBar.Phi());
            t_antinu_e.push_back(nuBar.E());
            t_nu.push_back(LorentzVector(nu.Pt(),nu.Eta(),nu.Phi(),nu.M()));
            t_nu_pt.push_back(nu.Pt());
            t_nu_eta.push_back(nu.Eta());
            t_nu_phi.push_back(nu.Phi());
            t_nu_e.push_back(nu.E());
            t_antiW.push_back(LorentzVector(W2.Pt(),W2.Eta(),W2.Phi(),W2.M()));
            t_antiW_pt.push_back(W2.Pt());
            t_antiW_eta.push_back(W2.Eta());
            t_antiW_phi.push_back(W2.Phi());
            t_antiW_e.push_back(W2.E());
            t_W.push_back(LorentzVector(W1.Pt(),W1.Eta(),W1.Phi(),W1.M()));
            t_W_pt.push_back(W1.Pt());
            t_W_eta.push_back(W1.Eta());
            t_W_phi.push_back(W1.Phi());
            t_W_e.push_back(W1.E());
            t_antitop.push_back(LorentzVector(t2.Pt(),t2.Eta(),t2.Phi(),t2.M()));
            t_antitop_pt.push_back(t2.Pt());
            t_antitop_eta.push_back(t2.Eta());
            t_antitop_phi.push_back(t2.Phi());
            t_antitop_e.push_back(t2.E());
            t_top.push_back(LorentzVector(t1.Pt(),t1.Eta(),t1.Phi(),t1.M()));
            t_top_pt.push_back(t1.Pt());
            t_top_eta.push_back(t1.Eta());
            t_top_phi.push_back(t1.Phi());
            t_top_e.push_back(t1.E());


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


            t_cos.push_back(TMath::Cos(theta1));
            t_cos.push_back(TMath::Cos(theta2));

            if(!isData)
            {
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
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta1)-TMath::Cos(gentheta1));
                    h_RecoMinusTruthCos->Fill(TMath::Cos(theta2)-TMath::Cos(gentheta2));
                    h_RecoMinusTruthMTT->Fill(ttbar.M() - genttbar.M());




                }
            }

        }
        //cout << "Am I here Now 3" <<endl;
        h_Nevents_ATS->Fill(1,theWeight);
        h_NeventsElMu_ATS->Fill(1,theWeight);
        h_ATS_ptLepElMu->Fill(PosLep.Pt(),theWeight);
        h_ATS_ptLepElMu->Fill(NegLep.Pt(),theWeight);
        h_ATS_etaLepElMu->Fill(PosLep.Eta(),theWeight);
        h_ATS_etaLepElMu->Fill(NegLep.Eta(),theWeight);
        h_ATS_mLepElMu->Fill(dilepton.M(),theWeight);
        h_ATS_METElMu->Fill(met.pt(),theWeight);
        h_ATS_NBJetsElMu->Fill(bjets.size(),theWeight);
        h_ATS_NJetsElMu->Fill(njets.size(),theWeight);

        h_NBJetsElMu->Fill(bjets.size(),theWeight);
        h_NBJetsDiLep->Fill(bjets.size(),theWeight);
    }

    }    
     t_outTree->Fill();
    
    
    
    //////////////////////GEN LEVEL////////////////
    if(!isData)
    {
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
        //        if(matchedPosMu.pt() != 0 )  cout << matchedPosMu.pt() << "muon match mother"<< posWMuMu.pdgId()<<"topMuMu Mother pdgId"<< topMuMu.pdgId() << endl;
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
            double w = sqrt(SF(t1.Pt())*SF(t2.Pt()));
            w_top = w_top*w;
            theWeight=theWeight*w;
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
        //        if(matchedPosEl.pt() != 0 )  cout << matchedPosEl.pt() << "Elon match mother"<< posWElEl.pdgId()<<"topElEl Mother pdgId"<< topElEl.pdgId() << endl;
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
            double w = sqrt(SF(t1.Pt())*SF(t2.Pt()));
            w_top = w_top*w;
            theWeight=theWeight*w;

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
        //        if(matchedPosEl.pt() != 0 )  cout << matchedPosEMEl.pt() << "Elon match mother"<< posWElMu.pdgId()<<"topElEl Mother pdgId"<< topElMu.pdgId() << endl;
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
            double w = sqrt(SF(t1.Pt())*SF(t2.Pt()));
            w_top = w_top*w;
            theWeight=theWeight*w;
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
        //        if(matchedPosEl.pt() != 0 )  cout << matchedPosEMMu.pt() << "Elon match mother"<< posWElMu2.pdgId()<<"topElEl Mother pdgId"<< topElMu2.pdgId() << endl;
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
            double w = sqrt(SF(t1.Pt())*SF(t2.Pt()));
            w_top = w_top*w;
            theWeight=theWeight*w;
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
    h_Weight->Fill(theWeight);



    if(DiMu)
    {
        if(!isDiMuon) return;
        h_NeventsMuMu_ALS->Fill(1,theWeight);
    }
    if(DiEl)
    {
        if(!isDiElectron) return;
        h_NeventsElEl_ALS->Fill(1,theWeight);
    }
    if(ElMu)
    {
        if(!isElMu && !isMuEl) return;
        h_NeventsElMu_ALS->Fill(1,theWeight);
    }

    h_Nevents_ALS->Fill(1,theWeight);

    /////////////////////2 Z Mass Veto //////////////////////////////////
    if(isDiMuon)
    {
        double w_ALS = w_posMu*w_negMu*w_mc*w_mm;
        PosLep.SetPtEtaPhiM(posMu.pt(),posMu.eta(),posMu.phi(),posMu.mass());
        NegLep.SetPtEtaPhiM(negMu.pt(),negMu.eta(),negMu.phi(),negMu.mass());
        DiLep = PosLep + NegLep;
        if(PosLep.Pt() > NegLep.Pt())
        {
            if(PosLep.Pt() < 25 ) return;
        }
        else
        {
            if(NegLep.Pt() < 25 ) return;
        }
        if(DiLep.M() < 20) return;
        h_Nevents_DiMu->Fill(1,w_ALS);;
        h_ALS_mLepNoVetoMuMu->Fill(DiLep.M(),w_ALS);
        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinMuMu_ALS->Fill(1,w_ALS);

        }
        if(!(DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NoutMuMu_ALS->Fill(1,w_ALS);;

            h_ALS_etaLMuMu->Fill(posMu.eta(),w_ALS);
            h_ALS_etaLMuMu->Fill(negMu.eta(),w_ALS);
            h_ALS_ptLMuMu->Fill(posMu.pt(),w_ALS);
            h_ALS_ptLMuMu->Fill(negMu.pt(),w_ALS);

            h_ALS_etaLDiLep->Fill(posMu.eta(),w_ALS);
            h_ALS_etaLDiLep->Fill(negMu.eta(),w_ALS);
            h_ALS_ptLDiLep->Fill(posMu.pt(),w_ALS);
            h_ALS_ptLDiLep->Fill(negMu.pt(),w_ALS);
        }

    }
    if(isDiElectron)
    {
        double w_ALS = w_posEl*w_negEl*w_mc*w_ee;
        PosLep.SetPtEtaPhiM(posEl.pt(),posEl.eta(),posEl.phi(),posEl.mass());
        NegLep.SetPtEtaPhiM(negEl.pt(),negEl.eta(),negEl.phi(),negEl.mass());
        DiLep = PosLep + NegLep;
        if(PosLep.Pt() > NegLep.Pt())
        {
            if(PosLep.Pt() < 25 ) return;
        }
        else
        {
            if(NegLep.Pt() < 25 ) return;
        }
        if(DiLep.M() < 20) return;
        h_Nevents_DiEl->Fill(1,w_ALS);;
        h_ALS_mLepNoVetoElEl->Fill(DiLep.M(),w_ALS);

        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinElEl_ALS->Fill(1,w_ALS);


        }
        if(!(DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NoutElEl_ALS->Fill(1,w_ALS);;

            h_ALS_etaLElEl->Fill(posEl.eta(),w_ALS);
            h_ALS_etaLElEl->Fill(negEl.eta(),w_ALS);
            h_ALS_ptLElEl->Fill(posEl.pt(),w_ALS);
            h_ALS_ptLElEl->Fill(negEl.pt(),w_ALS);

            h_ALS_etaLDiLep->Fill(posEl.eta(),w_ALS);
            h_ALS_etaLDiLep->Fill(negEl.eta(),w_ALS);
            h_ALS_ptLDiLep->Fill(posEl.pt(),w_ALS);
            h_ALS_ptLDiLep->Fill(negEl.pt(),w_ALS);
        }

    }
    if(isElMu)
    {
        double w_ALS = w_posEl*w_negMu*w_mc*w_em;
        PosLep.SetPtEtaPhiM(posEl.pt(),posEl.eta(),posEl.phi(),posEl.mass());
        NegLep.SetPtEtaPhiM(negMu.pt(),negMu.eta(),negMu.phi(),negMu.mass());
        DiLep = PosLep + NegLep;
        if(PosLep.Pt() > NegLep.Pt())
        {
            if(PosLep.Pt() < 25 ) return;
        }
        else
        {
            if(NegLep.Pt() < 25 ) return;
        }
        if(DiLep.M() < 20) return;
        h_Nevents_ElMu->Fill(1,w_ALS);;
        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinElMu_ALS->Fill(1,w_ALS);;


        }
        //        if((DiLep.M() < 106 && DiLep.M() > 76)) return;
        h_ALS_mLepNoVetoElMu->Fill(DiLep.M(),w_ALS);
        h_ALS_etaLElMu->Fill(posEl.eta(),w_ALS);
        h_ALS_etaLElMu->Fill(negMu.eta(),w_ALS);
        h_ALS_ptLElMu->Fill(posEl.pt(),w_ALS);
        h_ALS_ptLElMu->Fill(negMu.pt(),w_ALS);

        h_ALS_etaLDiLep->Fill(posEl.eta(),w_ALS);
        h_ALS_etaLDiLep->Fill(negMu.eta(),w_ALS);
        h_ALS_ptLDiLep->Fill(posEl.pt(),w_ALS);
        h_ALS_ptLDiLep->Fill(negMu.pt(),w_ALS);

    }
    if(isMuEl)
    {
        double w_ALS = w_posMu*w_negEl*w_mc*w_em;
        PosLep.SetPtEtaPhiM(posMu.pt(),posMu.eta(),posMu.phi(),posMu.mass());
        NegLep.SetPtEtaPhiM(negEl.pt(),negEl.eta(),negEl.phi(),negEl.mass());
        DiLep = PosLep + NegLep;
        if(PosLep.Pt() > NegLep.Pt())
        {
            if(PosLep.Pt() < 25 ) return;
        }
        else
        {
            if(NegLep.Pt() < 25 ) return;
        }
        if(DiLep.M() < 20) return;
        h_Nevents_ElMu->Fill(1,w_ALS);;

        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinElMu_ALS->Fill(1,w_ALS);;
        }
        //        if((DiLep.M() < 106 && DiLep.M() > 76)) return;
        h_ALS_mLepNoVetoElMu->Fill(DiLep.M(),w_ALS);

        h_ALS_etaLElMu->Fill(posMu.eta(),w_ALS);
        h_ALS_etaLElMu->Fill(negEl.eta(),w_ALS);
        h_ALS_ptLElMu->Fill(posMu.pt(),w_ALS);
        h_ALS_ptLElMu->Fill(negEl.pt(),w_ALS);

        h_ALS_etaLDiLep->Fill(posMu.eta(),w_ALS);
        h_ALS_etaLDiLep->Fill(negEl.eta(),w_ALS);
        h_ALS_ptLDiLep->Fill(posMu.pt(),w_ALS);
        h_ALS_ptLDiLep->Fill(negEl.pt(),w_ALS);

    }

    ////////////////3 Minimal Jet Multiplicity //////////////////
    if(njets.size() < 2) return;
    //    //++n_after2Jets;
    h_Nevents_AJS->Fill(1,theWeight);
    //After jet selection
    if(njets.size() > 0)
    {
        if(isDiMuon)
        {

            h_AJS_mLepNoVetoMuMu->Fill(DiLep.M(),theWeight);
            if(bjets.size()>=1)h_AoneBS_mLepNoVetoMuMu->Fill(DiLep.M(),theWeight);;

            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                h_NinMuMu_AJS->Fill(1,theWeight);;

                if(bjets.size()>=1)h_NinMuMu_NoMET_AoneBS->Fill(1,theWeight);;
                if(bjets.size()>=2)h_NinMuMu_NoMET_ABS->Fill(1,theWeight);;
                if(bjets.size()==2)h_NinMuMu_NoMET_ATS->Fill(1,theWeight);;

            }
            if(!(DiLep.M() < 106 && DiLep.M() > 76))
            {
                h_NoutMuMu_AJS->Fill(1,theWeight);;
                if(bjets.size()>=1)h_NoutMuMu_NoMET_AoneBS->Fill(1,theWeight);;
                if(bjets.size()>=2)h_NoutMuMu_NoMET_ABS->Fill(1,theWeight);;
                if(bjets.size()==2)h_NoutMuMu_NoMET_ATS->Fill(1,theWeight);;


                h_NeventsMuMu_AJS->Fill(1,theWeight);
                h_AJS_ptLepMuMu->Fill(posMu.pt(),theWeight);
                h_AJS_ptLepDiLep->Fill(posMu.pt(),theWeight);
                h_AJS_ptLepMuMu->Fill(negMu.pt(),theWeight);
                h_AJS_ptLepDiLep->Fill(negMu.pt(),theWeight);

                h_AJS_etaLepMuMu->Fill(posMu.eta(),theWeight);
                h_AJS_etaLepDiLep->Fill(posMu.eta(),theWeight);
                h_AJS_etaLepMuMu->Fill(negMu.eta(),theWeight);
                h_AJS_etaLepDiLep->Fill(negMu.eta(),theWeight);
            }
        }
        if(isDiElectron)
        {
            h_AJS_mLepNoVetoElEl->Fill(DiLep.M(),theWeight);
            if(bjets.size()>=1)h_AoneBS_mLepNoVetoElEl->Fill(DiLep.M(),theWeight);;
            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                h_NinElEl_AJS->Fill(1,theWeight);;

                if(bjets.size()>=1)h_NinElEl_NoMET_AoneBS->Fill(1,theWeight);;
                if(bjets.size()>=2)h_NinElEl_NoMET_ABS->Fill(1,theWeight);;
                if(bjets.size()==2)h_NinElEl_NoMET_ATS->Fill(1,theWeight);;


            }
            if(!(DiLep.M() < 106 && DiLep.M() > 76))
            {
                h_NoutElEl_AJS->Fill(1,theWeight);;
                if(bjets.size()>=1)h_NoutElEl_NoMET_AoneBS->Fill(1,theWeight);;
                if(bjets.size()>=2)h_NoutElEl_NoMET_ABS->Fill(1,theWeight);;
                if(bjets.size()==2)h_NoutElEl_NoMET_ATS->Fill(1,theWeight);;


                h_NeventsElEl_AJS->Fill(1,theWeight);
                h_AJS_ptLepElEl->Fill(posEl.pt(),theWeight);
                h_AJS_ptLepDiLep->Fill(posEl.pt(),theWeight);
                h_AJS_ptLepElEl->Fill(negEl.pt(),theWeight);
                h_AJS_ptLepDiLep->Fill(negEl.pt(),theWeight);

                h_AJS_etaLepElEl->Fill(posEl.eta(),theWeight);
                h_AJS_etaLepDiLep->Fill(posEl.eta(),theWeight);
                h_AJS_etaLepElEl->Fill(negEl.eta(),theWeight);
                h_AJS_etaLepDiLep->Fill(negEl.eta(),theWeight);
            }


        }
        if(isElMu)
        {
            h_AJS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);
            if(bjets.size()>=1)h_AoneBS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);;
            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                h_NinElMu_AJS->Fill(1,theWeight);;

                if(bjets.size()>=1)h_NinElMu_NoMET_AoneBS->Fill(1,theWeight);;
                if(bjets.size()>=2)h_NinElMu_NoMET_ABS->Fill(1,theWeight);;
                if(bjets.size()==2)h_NinElMu_NoMET_ATS->Fill(1,theWeight);;


            }
            h_NeventsElMu_AJS->Fill(1,theWeight);
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
            h_AJS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);
            if(bjets.size()>=1)h_AoneBS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);;
            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                h_NinElMu_AJS->Fill(1,theWeight);;
                if(bjets.size()>=1)h_NinElMu_NoMET_AoneBS->Fill(1,theWeight);;
                if(bjets.size()>=2)h_NinElMu_NoMET_ABS->Fill(1,theWeight);;
                if(bjets.size()==2)h_NinElMu_NoMET_ATS->Fill(1,theWeight);;

            }
            h_NeventsElMu_AJS->Fill(1,theWeight);
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


    ////////////////////4 METS PF1//////////////////////////

    if(isDiMuon && met.pt() < 40) return;
    if(isDiElectron && met.pt() < 40) return;
    //if((isElMu || isMuEl) && met.pt() <0) return;
    //++n_afterMet;
    h_Nevents_AMS->Fill(1,theWeight);;

    if(isDiMuon)
    {
        h_AMS_mLepNoVetoMuMu->Fill(DiLep.M(),theWeight);
        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinMuMu_AMS->Fill(1,theWeight);;

        }
        if(!(DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NoutMuMu_AMS->Fill(1,theWeight);;
            h_NeventsMuMu_AMS->Fill(1,theWeight);
            h_AMS_mLepMuMu->Fill(DiLep.M(),theWeight);
            h_AMS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_AMS_ptLepMuMu->Fill(posMu.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(posMu.pt(),theWeight);
            h_AMS_ptLepMuMu->Fill(negMu.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(negMu.pt(),theWeight);
            h_AMS_METMuMu->Fill(met.pt(),theWeight);
            h_AMS_METDiLep->Fill(met.pt(),theWeight);
        }
    }
    if(isDiElectron)
    {
        h_AMS_mLepNoVetoElEl->Fill(DiLep.M(),theWeight);
        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinElEl_AMS->Fill(1,theWeight);;

        }
        if(!(DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NoutElEl_AMS->Fill(1,theWeight);;
            h_NeventsElEl_AMS->Fill(1,theWeight);
            h_AMS_mLepElEl->Fill(DiLep.M(),theWeight);
            h_AMS_mLepDiLep->Fill(DiLep.M(),theWeight);
            h_AMS_ptLepElEl->Fill(posEl.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(posEl.pt(),theWeight);
            h_AMS_ptLepElEl->Fill(negEl.pt(),theWeight);
            h_AMS_ptLepDiLep->Fill(negEl.pt(),theWeight);
            h_AMS_METElEl->Fill(met.pt(),theWeight);
            h_AMS_METDiLep->Fill(met.pt(),theWeight);
        }


    }
    if(isElMu)
    {
        h_AMS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);
        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinElMu_AMS->Fill(1,theWeight);;

        }
        h_NeventsElMu_AMS->Fill(1,theWeight);
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
        h_AMS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);
        if((DiLep.M() < 106 && DiLep.M() > 76))
        {
            h_NinElMu_AMS->Fill(1,theWeight);;

        }
        h_NeventsElMu_AMS->Fill(1,theWeight);
        h_AMS_mLepElMu->Fill(DiLep.M(),theWeight);
        h_AMS_mLepDiLep->Fill(DiLep.M(),theWeight);
        h_AMS_ptLepElMu->Fill(posMu.pt(),theWeight);
        h_AMS_ptLepDiLep->Fill(posMu.pt(),theWeight);
        h_AMS_ptLepElMu->Fill(negEl.pt(),theWeight);
        h_AMS_ptLepDiLep->Fill(negEl.pt(),theWeight);
        h_AMS_METElMu->Fill(met.pt(),theWeight);
        h_AMS_METDiLep->Fill(met.pt(),theWeight);


    }

    /////////////////5 B-jets //////////////
    ////++n_after2BJets;
    if(bjets.size() > 0)
    {
        if(isDiMuon)
        {
            //if(bjets.size() >= 1) h_ABS_mLepNoVetoMuMu->Fill(DiLep.M(),theWeight);
            if(bjets.size() >= 1) h_ABS_mLepNoVetoDiLep->Fill(DiLep.M(),theWeight);
            if(bjets.size() >= 2) h_ABS_mLepNoVetoMuMu->Fill(DiLep.M(),theWeight);
            if(bjets.size() == 2) h_ATS_mLepNoVetoMuMu->Fill(DiLep.M(),theWeight);


            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                if(bjets.size() >= 1)h_NinMuMu_AoneBS->Fill(1,theWeight);;
                if(bjets.size() >= 2)h_NinMuMu_ABS->Fill(1,theWeight);
                if(bjets.size() == 2)h_NinMuMu_ATS->Fill(1,theWeight);

                return;
            }
            if(!(DiLep.M() < 106 && DiLep.M() > 76))
            {
                if(bjets.size() >= 1)h_NoutMuMu_AoneBS->Fill(1,theWeight);;
                if(bjets.size() >= 2)h_NoutMuMu_ABS->Fill(1,theWeight);;
                if(bjets.size() == 2) h_NoutMuMu_ATS->Fill(1,theWeight);
                if(bjets.size() >= 2)
                {
                    h_NeventsMuMu_ABS->Fill(1,theWeight);
                    h_ABS_NBJetsMuMu->Fill(bjets.size());
                    h_ABS_NJetsMuMu->Fill(njets.size());
                    h_ABS_NBJetsDiLep->Fill(bjets.size());
                    h_ABS_NJetsDiLep->Fill(njets.size());
                    h_ABS_etaLeadingJetMuMu->Fill(bjets.at(0).eta(),theWeight);
                    h_ABS_ptLeadingJetMuMu->Fill(bjets.at(0).pt(),theWeight);
                    h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
                    h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);


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
            }
        }
        if(isDiElectron)
        {
            if(bjets.size() >= 2) h_ABS_mLepNoVetoElEl->Fill(DiLep.M(),theWeight);
            if(bjets.size() >= 1) h_ABS_mLepNoVetoDiLep->Fill(DiLep.M(),theWeight);
            if(bjets.size() == 2) h_ATS_mLepNoVetoElEl->Fill(DiLep.M(),theWeight);

            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                if(bjets.size() >= 1)h_NinElEl_AoneBS->Fill(1,theWeight);;
                if(bjets.size() >= 2)h_NinElEl_ABS->Fill(1,theWeight);;
                if(bjets.size() == 2)h_NinElEl_ATS->Fill(1,theWeight);;

                return;
            }
            if(!(DiLep.M() < 106 && DiLep.M() > 76))
            {
                if(bjets.size() >= 1)h_NoutElEl_AoneBS->Fill(1,theWeight);;
                if(bjets.size() >= 2)h_NoutElEl_ABS->Fill(1,theWeight);;
                if(bjets.size() == 2)h_NoutElEl_ATS->Fill(1,theWeight);;

                if(bjets.size() >= 2)
                {
                    h_NeventsElEl_ABS->Fill(1,theWeight);
                    h_ABS_NBJetsElEl->Fill(bjets.size());
                    h_ABS_NJetsElEl->Fill(njets.size());
                    h_ABS_NBJetsDiLep->Fill(bjets.size());
                    h_ABS_NJetsDiLep->Fill(njets.size());
                    h_ABS_etaLeadingJetElEl->Fill(bjets.at(0).eta(),theWeight);
                    h_ABS_ptLeadingJetElEl->Fill(bjets.at(0).pt(),theWeight);
                    h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
                    h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);

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
            }


        }
        if(isElMu)
        {
            if(bjets.size() >= 2) h_ABS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);
            if(bjets.size() >= 1) h_ABS_mLepNoVetoDiLep->Fill(DiLep.M(),theWeight);
            if(bjets.size() == 2) h_ATS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);

            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                if(bjets.size() >= 1)h_NinElMu_AoneBS->Fill(1,theWeight);;
                if(bjets.size() >= 2)h_NinElMu_ABS->Fill(1,theWeight);;
                if(bjets.size() == 2)h_NinElMu_ATS->Fill(1,theWeight);;

            }
            if(bjets.size() >= 2)
            {
                h_NeventsElMu_ABS->Fill(1,theWeight);
                h_ABS_NBJetsElMu->Fill(bjets.size());
                h_ABS_NJetsElMu->Fill(njets.size());
                h_ABS_NBJetsDiLep->Fill(bjets.size());
                h_ABS_NJetsDiLep->Fill(njets.size());
                h_ABS_etaLeadingJetElMu->Fill(bjets.at(0).eta(),theWeight);
                h_ABS_ptLeadingJetElMu->Fill(bjets.at(0).pt(),theWeight);
                h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
                h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);

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

        }
        if(isMuEl)
        {
            if(bjets.size() >= 2) h_ABS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);
            if(bjets.size() >= 1) h_ABS_mLepNoVetoDiLep->Fill(DiLep.M(),theWeight);
            if(bjets.size() == 2) h_ATS_mLepNoVetoElMu->Fill(DiLep.M(),theWeight);

            if((DiLep.M() < 106 && DiLep.M() > 76))
            {
                if(bjets.size() >= 1)h_NinElMu_AoneBS->Fill(1,theWeight);;
                if(bjets.size() >= 2)h_NinElMu_ABS->Fill(1,theWeight);;
                if(bjets.size() == 2)h_NinElMu_ATS->Fill(1,theWeight);;

            }
            if(bjets.size() >= 2)
            {
                h_NeventsElMu_ABS->Fill(1,theWeight);
                h_ABS_NBJetsElMu->Fill(bjets.size());
                h_ABS_NJetsElMu->Fill(njets.size());
                h_ABS_NBJetsDiLep->Fill(bjets.size());
                h_ABS_NJetsDiLep->Fill(njets.size());
                h_ABS_etaLeadingJetElMu->Fill(bjets.at(0).eta(),theWeight);
                h_ABS_ptLeadingJetElMu->Fill(bjets.at(0).pt(),theWeight);
                h_ABS_etaLeadingJetDiLep->Fill(bjets.at(0).eta(),theWeight);
                h_ABS_ptLeadingJetDiLep->Fill(bjets.at(0).pt(),theWeight);

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
    }
    if(bjets.size() < 2) return;
    h_Nevents_ABS->Fill(1,theWeight);;





   



    //    cout << "isPythia: "<<isPythia << endl;



    t_lumi = 0;
    t_bunch = 0;
    t_Event = 0;
    t_Run = 0;


}

void MiniAnalyzer::beginJob()
{
    TH1::SetDefaultSumw2();

    h_cosMuMu = fs->make<TH1F>("h_cosMuMu",";cos(#theta);",22,-1.1,1.1);
    h_cosElEl = fs->make<TH1F>("h_cosElEl",";cos(#theta);",22,-1.1,1.1);
    h_cosElMu = fs->make<TH1F>("h_cosElMu",";cos(#theta);",22,-1.1,1.1);
    h_cosDiLep = fs->make<TH1F>("h_cosDiLep",";cos(#theta);",22,-1.1,1.1);
    h_cosGenElMu = fs->make<TH1F>("h_cosGenElMu",";cos(#theta);",22,-1.1,1.1);
    h_cosGenMuMu = fs->make<TH1F>("h_cosGenMuMu",";cos(#theta);",22,-1.1,1.1);
    h_cosGenElEl = fs->make<TH1F>("h_cosGenElEl",";cos(#theta);",22,-1.1,1.1);
    h_cosGen = fs->make<TH1F>("h_cosGen",";cos(#theta);",22,-1.1,1.1);
    h_GenTTbarM = fs->make<TH1F>("h_GenTTbarM",";M_{TTbar};",100,90,1300);
    h_etaMu = fs->make<TH1F>("h_EtaMu",";#eta_{l};",100,-3,3);
    h_TTbarM = fs->make<TH1F>("h_TTbarM",";M_{TTbar};",100,0,1300);
    h_PtMu = fs->make<TH1F>("h_PtMu",";Mu_{Pt};",100,0.,200.);
    h_NPV = fs->make<TH1F>("h_NPV",";NPV;",50,0,50);
    h_NBJetsMuMu = fs->make<TH1F>("h_NBJetsMuMu",";N_{BJets};",8,0,8);
    h_NBJetsElEl = fs->make<TH1F>("h_NBJetsElEl",";N_{BJets};",8,0,8);
    h_NBJetsElMu = fs->make<TH1F>("h_NBJetsElMu",";N_{BJets};",8,0,8);
    h_NBJetsDiLep = fs->make<TH1F>("h_NBJetsDiLep",";N_{BJets};",8,0,8);
    h_yTMuMu = fs->make<TH1F>("h_yTMuMu",";y^{Top};",35,-3,3);
    h_yTElEl = fs->make<TH1F>("h_yTElEl",";y^{Top};",35,-3,3);
    h_yTElMu = fs->make<TH1F>("h_yTElMu",";y^{Top};",35,-3,3);
    h_yTDiLep = fs->make<TH1F>("h_yTDiLep",";y^{Top};",35,-3,3);
    h_yWMuMu = fs->make<TH1F>("h_yWMuMu",";y^{W};",35,-3,3);
    h_yWElEl = fs->make<TH1F>("h_yWElEl",";y^{W};",35,-3,3);
    h_yWElMu = fs->make<TH1F>("h_yWElMu",";y^{W};",35,-3,3);
    h_yWDiLep = fs->make<TH1F>("h_yWDiLep",";y^{W};",35,-3,3);
    h_ptTMuMu= fs->make<TH1F>("h_ptTMuMu",";Pt^{Top};",35,0.,300.);
    h_ptTElEl= fs->make<TH1F>("h_ptTElEl",";Pt^{Top};",35,0.,300.);
    h_ptTElMu= fs->make<TH1F>("h_ptTElMu",";Pt^{Top};",35,0.,300.);
    h_ptTDiLep= fs->make<TH1F>("h_ptTDiLep",";Pt^{Top};",35,0.,500.);
    h_ptWMuMu= fs->make<TH1F>("h_ptWMuMu",";Pt^{W};",35,0.,1000.);
    h_ptWElEl= fs->make<TH1F>("h_ptWElEl",";Pt^{W};",35,0.,1000.);
    h_ptWElMu= fs->make<TH1F>("h_ptWElMu",";Pt^{W};",35,0.,1000.);
    h_ptWDiLep= fs->make<TH1F>("h_ptWDiLep",";Pt^{W};",35,0.,1000.);
    h_mTTbarMuMu =fs->make<TH1F>("h_mTTbarMuMu",";M^{tt};",35,0,1300);
    h_mTTbarElEl =fs->make<TH1F>("h_mTTbarElEl",";M^{tt};",35,0,1300);
    h_mTTbarElMu =fs->make<TH1F>("h_mTTbarElMu",";M^{tt};",35,0,1300);

    h_ATS_etaLepMuMu = fs->make<TH1F>("h_ATS_EtaLepMuMu",";#eta_{l};",35,-3,3);
    h_ATS_etaLepElEl = fs->make<TH1F>("h_ATS_EtaLepElEl",";#eta_{l};",35,-3,3);
    h_ATS_etaLepElMu = fs->make<TH1F>("h_ATS_EtaLepElMu",";#eta_{l};",35,-3,3);
    h_ATS_ptLepMuMu = fs->make<TH1F>("h_ATS_PtLepMuMu",";Pt_{l};",35,0.,300.);
    h_ATS_ptLepElMu = fs->make<TH1F>("h_ATS_PtLepElMu",";Pt_{l};",35,0.,300.);
    h_ATS_ptLepElEl = fs->make<TH1F>("h_ATS_PtLepElEl",";Pt_{l};",35,0.,300.);
    h_ATS_mLepMuMu = fs->make<TH1F>("h_ATS_mLepMuMu",";M_{ll};",35,0.,300.);
    h_ATS_mLepElEl = fs->make<TH1F>("h_ATS_mLepElEl",";M_{ll};",35,0.,300.);
    h_ATS_mLepElMu = fs->make<TH1F>("h_ATS_mLepElMu",";M_{ll};",35,0.,300.);

    h_ATS_METMuMu = fs->make<TH1F>("h_ATS_METMuMu",";MET;",35,0.,300.);
    h_ATS_METElEl = fs->make<TH1F>("h_ATS_METElEl",";MET;",35,0.,300.);
    h_ATS_METElMu = fs->make<TH1F>("h_ATS_METElMu",";MET;",35,0.,300.);
    h_ATS_NBJetsMuMu = fs->make<TH1F>("h_ATS_NBJetsMuMu",";N_{BJets};",8,0,8);
    h_ATS_NBJetsElEl = fs->make<TH1F>("h_ATS_NBJetsElEl",";N_{BJets};",8,0,8);
    h_ATS_NBJetsElMu = fs->make<TH1F>("h_ATS_NBJetsElMu",";N_{BJets};",8,0,8);
    h_ATS_NJetsMuMu = fs->make<TH1F>("h_ATS_NJetsMuMu",";N_{BJets};",8,0,8);
    h_ATS_NJetsElEl = fs->make<TH1F>("h_ATS_NJetsElEl",";N_{BJets};",8,0,8);
    h_ATS_NJetsElMu = fs->make<TH1F>("h_ATS_NJetsElMu",";N_{BJets};",8,0,8);
    h_NeventsElEl_ALS = fs->make<TH1F>("h_NeventsElEl_ALS","",10,-3,3);
    h_NeventsElMu_ALS = fs->make<TH1F>("h_NeventsElMu_ALS","",10,-3,3);
    h_NeventsMuMu_ALS = fs->make<TH1F>("h_NeventsMuMu_ALS","",10,-3,3);
    h_NeventsElEl_AJS = fs->make<TH1F>("h_NeventsElEl_AJS","",10,-3,3);
    h_NeventsElMu_AJS = fs->make<TH1F>("h_NeventsElMu_AJS","",10,-3,3);
    h_NeventsMuMu_AJS = fs->make<TH1F>("h_NeventsMuMu_AJS","",10,-3,3);
    h_NeventsElEl_AMS = fs->make<TH1F>("h_NeventsElEl_AMS","",10,-3,3);
    h_NeventsElMu_AMS = fs->make<TH1F>("h_NeventsElMu_AMS","",10,-3,3);
    h_NeventsMuMu_AMS = fs->make<TH1F>("h_NeventsMuMu_AMS","",10,-3,3);
    h_NeventsElEl_ABS = fs->make<TH1F>("h_NeventsElEl_ABS","",10,-3,3);
    h_NeventsElMu_ABS = fs->make<TH1F>("h_NeventsElMu_ABS","",10,-3,3);
    h_NeventsMuMu_ABS = fs->make<TH1F>("h_NeventsMuMu_ABS","",10,-3,3);
    h_NeventsElEl_ATS = fs->make<TH1F>("h_NeventsElEl_ATS","",10,-3,3);
    h_NeventsElMu_ATS = fs->make<TH1F>("h_NeventsElMu_ATS","",10,-3,3);
    h_NeventsMuMu_ATS = fs->make<TH1F>("h_NeventsMuMu_ATS","",10,-3,3);
    h_NeventsElEl_top = fs->make<TH1F>("h_NeventsElEl_top","",10,-3,3);
    h_NeventsElMu_top = fs->make<TH1F>("h_NeventsElMu_top","",10,-3,3);
    h_NeventsMuMu_top = fs->make<TH1F>("h_NeventsMuMu_top","",10,-3,3);
    h_Nevents_ATS = fs->make<TH1F>("h_Nevents_ATS","",10,-3,3);


    h_ALS_etaLMuMu = fs->make<TH1F>("h_ALS_EtaLepMuMu",";#eta_{l};",35,-3,3);
    h_ALS_etaLElEl = fs->make<TH1F>("h_ALS_EtaLepElEl",";#eta_{l};",35,-3,3);
    h_ALS_etaLElMu = fs->make<TH1F>("h_ALS_EtaLepElMu",";#eta_{l};",35,-3,3);
    h_ALS_etaLDiLep = fs->make<TH1F>("h_ALS_EtaLepDiLep",";#eta_{l};",35,-3,3);
    h_ALS_ptLMuMu = fs->make<TH1F>("h_ALS_PtLepMuMu",";Pt_{l};",35,0.,300.);
    h_ALS_ptLElMu = fs->make<TH1F>("h_ALS_PtLepElMu",";Pt_{l};",35,0.,300.);
    h_ALS_ptLElEl = fs->make<TH1F>("h_ALS_PtLepElEl",";Pt_{l};",35,0.,300.);
    h_ALS_ptLDiLep = fs->make<TH1F>("h_ALS_PtLepDiLep",";Pt_{l};",35,0.,300.);

    h_AJS_etaLepMuMu = fs->make<TH1F>("h_AJS_EtaLepMuMu",";#eta_{l};",35,-3,3);
    h_AJS_etaLepElEl = fs->make<TH1F>("h_AJS_EtaLepElEl",";#eta_{l};",35,-3,3);
    h_AJS_etaLepElMu = fs->make<TH1F>("h_AJS_EtaLepElMu",";#eta_{l};",35,-3,3);
    h_AJS_etaLepDiLep = fs->make<TH1F>("h_AJS_EtaLepDiLep",";#eta_{l};",35,-3,3);
    h_AJS_ptLepMuMu = fs->make<TH1F>("h_AJS_PtLepMuMu",";Pt_{l};",35,0.,300.);
    h_AJS_ptLepElMu = fs->make<TH1F>("h_AJS_PtLepElMu",";Pt_{l};",35,0.,300.);
    h_AJS_ptLepElEl = fs->make<TH1F>("h_AJS_PtLepElEl",";Pt_{l};",35,0.,300.);
    h_AJS_ptLepDiLep = fs->make<TH1F>("h_AJS_PtLepDiLep",";Pt_{l};",35,0.,300.);

    h_ABS_etaLepMuMu = fs->make<TH1F>("h_ABS_EtaLepMuMu",";#eta_{l};",35,-3,3);
    h_ABS_etaLepElEl = fs->make<TH1F>("h_ABS_EtaLepElEl",";#eta_{l};",35,-3,3);
    h_ABS_etaLepElMu = fs->make<TH1F>("h_ABS_EtaLepElMu",";#eta_{l};",35,-3,3);
    h_ABS_etaLepDiLep = fs->make<TH1F>("h_ABS_EtaLepDiLep",";#eta_{l};",35,-3,3);
    h_ABS_ptLepMuMu = fs->make<TH1F>("h_ABS_PtLepMuMu",";Pt_{l};",35,0.,300.);
    h_ABS_ptLepElMu = fs->make<TH1F>("h_ABS_PtLepElMu",";Pt_{l};",35,0.,300.);
    h_ABS_ptLepElEl = fs->make<TH1F>("h_ABS_PtLepElEl",";Pt_{l};",35,0.,300.);
    h_ABS_ptLepDiLep = fs->make<TH1F>("h_ABS_PtLepDiLep",";Pt_{l};",35,0.,300.);
    h_ABS_mLepMuMu = fs->make<TH1F>("h_ABS_mLepMuMu",";M_{ll};",35,0.,300.);
    h_ABS_mLepElEl = fs->make<TH1F>("h_ABS_mLepElEl",";M_{ll};",35,0.,300.);
    h_ABS_mLepElMu = fs->make<TH1F>("h_ABS_mLepElMu",";M_{ll};",35,0.,300.);
    h_ABS_mLepDiLep = fs->make<TH1F>("h_ABS_mLepDiLep",";M_{ll};",35,0.,300.);
    h_ABS_mLepNoVetoMuMu = fs->make<TH1F>("h_ABS_mLepNoVetoMuMu",";M_{ll};",300,0.,300.);
    h_ABS_mLepNoVetoElEl = fs->make<TH1F>("h_ABS_mLepNoVetoElEl",";M_{ll};",300,0.,300.);
    h_ABS_mLepNoVetoElMu = fs->make<TH1F>("h_ABS_mLepNoVetoElMu",";M_{ll};",300,0.,300.);
    h_ABS_mLepNoVetoDiLep = fs->make<TH1F>("h_ABS_mLepNoVetoDiLep",";M_{ll};",300,0.,300.);
    h_AoneBS_mLepNoVetoMuMu = fs->make<TH1F>("h_AoneBS_mLepNoVetoMuMu",";M_{ll};",300,0.,300.);
    h_AoneBS_mLepNoVetoElEl = fs->make<TH1F>("h_AoneBS_mLepNoVetoElEl",";M_{ll};",300,0.,300.);
    h_AoneBS_mLepNoVetoElMu = fs->make<TH1F>("h_AoneBS_mLepNoVetoElMu",";M_{ll};",300,0.,300.);
    h_ALS_mLepNoVetoMuMu = fs->make<TH1F>("h_ALS_mLepNoVetoMuMu",";M_{ll};",300,0.,300.);
    h_ALS_mLepNoVetoElEl = fs->make<TH1F>("h_ALS_mLepNoVetoElEl",";M_{ll};",300,0.,300.);
    h_ALS_mLepNoVetoElMu = fs->make<TH1F>("h_ALS_mLepNoVetoElMu",";M_{ll};",300,0.,300.);
    h_AJS_mLepNoVetoMuMu = fs->make<TH1F>("h_AJS_mLepNoVetoMuMu",";M_{ll};",300,0.,300.);
    h_AJS_mLepNoVetoElEl = fs->make<TH1F>("h_AJS_mLepNoVetoElEl",";M_{ll};",300,0.,300.);
    h_AJS_mLepNoVetoElMu = fs->make<TH1F>("h_AJS_mLepNoVetoElMu",";M_{ll};",300,0.,300.);
    h_AMS_mLepNoVetoMuMu = fs->make<TH1F>("h_AMS_mLepNoVetoMuMu",";M_{ll};",300,0.,300.);
    h_AMS_mLepNoVetoElEl = fs->make<TH1F>("h_AMS_mLepNoVetoElEl",";M_{ll};",300,0.,300.);
    h_AMS_mLepNoVetoElMu = fs->make<TH1F>("h_AMS_mLepNoVetoElMu",";M_{ll};",300,0.,300.);


    h_ABS_METMuMu = fs->make<TH1F>("h_ABS_METMuMu",";MET;",35,0.,300.);
    h_ABS_METElEl = fs->make<TH1F>("h_ABS_METElEl",";MET;",35,0.,300.);
    h_ABS_METElMu = fs->make<TH1F>("h_ABS_METElMu",";MET;",35,0.,300.);
    h_ABS_METDiLep = fs->make<TH1F>("h_ABS_METDiLep",";MET;",35,0.,300.);
    h_ABS_NBJetsMuMu = fs->make<TH1F>("h_ABS_NBJetsMuMu",";N_{BJets};",8,0,8);
    h_ABS_NBJetsElEl = fs->make<TH1F>("h_ABS_NBJetsElEl",";N_{BJets};",8,0,8);
    h_ABS_NBJetsElMu = fs->make<TH1F>("h_ABS_NBJetsElMu",";N_{BJets};",8,0,8);
    h_ABS_NBJetsDiLep = fs->make<TH1F>("h_ABS_NBJetsDiLep",";N_{BJets};",8,0,8);
    h_ABS_NJetsMuMu = fs->make<TH1F>("h_ABS_NJetsMuMu",";N_{BJets};",8,0,8);
    h_ABS_NJetsElEl = fs->make<TH1F>("h_ABS_NJetsElEl",";N_{BJets};",8,0,8);
    h_ABS_NJetsElMu = fs->make<TH1F>("h_ABS_NJetsElMu",";N_{BJets};",8,0,8);
    h_ABS_NJetsDiLep = fs->make<TH1F>("h_ABS_NJetsDiLep",";N_{BJets};",8,0,8);
    h_ABS_etaLeadingJetMuMu = fs->make<TH1F>("h_ABS_etaLeadingJetMuMu",";#eta_{l};",35,-3,3);
    h_ABS_etaLeadingJetElEl = fs->make<TH1F>("h_ABS_EtaLeadingJetElEl",";#eta_{l};",35,-3,3);
    h_ABS_etaLeadingJetElMu = fs->make<TH1F>("h_ABS_EtaLeadingJetElMu",";#eta_{l};",35,-3,3);
    h_ABS_etaLeadingJetDiLep = fs->make<TH1F>("h_ABS_EtaLeadingJetDiLep",";#eta_{l};",35,-3,3);
    h_ABS_ptLeadingJetMuMu = fs->make<TH1F>("h_ABS_PtLeadingJetMuMu",";Pt_{l};",35,0.,300.);
    h_ABS_ptLeadingJetElMu = fs->make<TH1F>("h_ABS_PtLeadingJetElMu",";Pt_{l};",35,0.,300.);
    h_ABS_ptLeadingJetElEl = fs->make<TH1F>("h_ABS_PtLeadingJetElEl",";Pt_{l};",35,0.,300.);
    h_ABS_ptLeadingJetDiLep = fs->make<TH1F>("h_ABS_ptLeadingJetDiLep",";Pt_{l};",35,0.,300.);

    h_AMS_ptLepMuMu = fs->make<TH1F>("h_AMS_PtLepMuMu",";Pt_{l};",35,0.,300.);
    h_AMS_ptLepElMu = fs->make<TH1F>("h_AMS_PtLepElMu",";Pt_{l};",35,0.,300.);
    h_AMS_ptLepElEl = fs->make<TH1F>("h_AMS_PtLepElEl",";Pt_{l};",35,0.,300.);
    h_AMS_ptLepDiLep = fs->make<TH1F>("h_AMS_PtLepDiLep",";Pt_{l};",35,0.,300.);
    h_AMS_mLepMuMu = fs->make<TH1F>("h_AMS_mLepMuMu",";M_{ll};",35,0.,300.);
    h_AMS_mLepElEl = fs->make<TH1F>("h_AMS_mLepElEl",";M_{ll};",35,0.,300.);
    h_AMS_mLepElMu = fs->make<TH1F>("h_AMS_mLepElMu",";M_{ll};",35,0.,300.);
    h_AMS_mLepDiLep = fs->make<TH1F>("h_AMS_mLepDiLep",";M_{ll};",35,0.,300.);
    h_AMS_METMuMu = fs->make<TH1F>("h_AMS_METMuMu",";MET;",35,0.,300.);
    h_AMS_METElEl = fs->make<TH1F>("h_AMS_METElEl",";MET;",35,0.,300.);
    h_AMS_METElMu = fs->make<TH1F>("h_AMS_METElMu",";MET;",35,0.,300.);
    h_AMS_METDiLep = fs->make<TH1F>("h_AMS_METDiLep",";MET;",35,0.,300.);

    h_truthRecoMuMu = fs->make<TH2F>("h_truthRecoMuMu","",22,-1.1,1.1,22,-1.1,1.1);
    h_truthRecoElEl = fs->make<TH2F>("h_truthRecoElEl","",22,-1.1,1.1,22,-1.1,1.1);
    h_truthRecoElMu = fs->make<TH2F>("h_truthRecoElMu","",22,-1.1,1.1,22,-1.1,1.1);
    h_truthRecoCos = fs->make<TH2F>("h_truthRecoCos","",22,-1.1,1.1,22,-1.1,1.1);
    h_truthRecoMTT = fs->make<TH2F>("h_truthRecoMTT","",50,0.,1000.,50,0.,1000.);
    h_RecoMinusTruthCos = fs->make<TH1F>("h_RecoMinusTruthCos","",100,-3,3);
    h_RecoMinusTruthMTT = fs->make<TH1F>("h_RecoMinusTruthMTT","",100,-400,400);
    h_Nevents = fs->make<TH1F>("h_Nevents","",10,-3,3);
    h_Nevents_weighted = fs->make<TH1F>("h_Nevents_weighted","",10,-3,3);
    h_Nevents_ABS = fs->make<TH1F>("h_Nevents_ABS","",10,-3,3);
    h_Nevents_ALS = fs->make<TH1F>("h_Nevents_ALS","",10,-3,3);
    h_Nevents_AMS = fs->make<TH1F>("h_Nevents_AMS","",10,-3,3);
    h_Nevents_AJS = fs->make<TH1F>("h_Nevents_AJS","",10,-3,3);
    h_Nevents_AVS = fs->make<TH1F>("h_Nevents_AVS","",10,-3,3);
    h_Nevents_AT = fs->make<TH1F>("h_Nevents_AT","",10,-3,3);
    h_Weight = fs->make<TH1F>("h_Weight","",50,-1.2,1.2);
    h_Nevents_DiEl = fs->make<TH1F>("h_Nevents_DiEl","",10,-3,3);
    h_Nevents_DiMu = fs->make<TH1F>("h_Nevents_DiMu","",10,-3,3);
    h_Nevents_ElMu = fs->make<TH1F>("h_Nevents_ElMu","",10,-3,3);
    h_Nevents_top = fs->make<TH1F>("h_Nevents_top","",10,-3,3);

    h_NinMuMu_ALS = fs->make<TH1F>("h_NinMuMu_ALS","",10,-3,3);
    h_NinElEl_ALS = fs->make<TH1F>("h_NinElEl_ALS","",10,-3,3);
    h_NinElMu_ALS = fs->make<TH1F>("h_NinElMu_ALS","",10,-3,3);
    h_NoutMuMu_ALS = fs->make<TH1F>("h_NoutMuMu_ALS","",10,-3,3);
    h_NoutElEl_ALS = fs->make<TH1F>("h_NoutElEl_ALS","",10,-3,3);
    h_NoutElMu_ALS = fs->make<TH1F>("h_NoutElMu_ALS","",10,-3,3);
    h_NinMuMu_AJS = fs->make<TH1F>("h_NinMuMu_AJS","",10,-3,3);
    h_NinElEl_AJS = fs->make<TH1F>("h_NinElEl_AJS","",10,-3,3);
    h_NinElMu_AJS = fs->make<TH1F>("h_NinElMu_AJS","",10,-3,3);
    h_NoutMuMu_AJS = fs->make<TH1F>("h_NoutMuMu_AJS","",10,-3,3);
    h_NoutElEl_AJS = fs->make<TH1F>("h_NoutElEl_AJS","",10,-3,3);
    h_NoutElMu_AJS = fs->make<TH1F>("h_NoutElMu_AJS","",10,-3,3);
    h_NinMuMu_AMS = fs->make<TH1F>("h_NinMuMu_AMS","",10,-3,3);
    h_NinElEl_AMS = fs->make<TH1F>("h_NinElEl_AMS","",10,-3,3);
    h_NinElMu_AMS = fs->make<TH1F>("h_NinElMu_AMS","",10,-3,3);
    h_NoutMuMu_AMS = fs->make<TH1F>("h_NoutMuMu_AMS","",10,-3,3);
    h_NoutElEl_AMS = fs->make<TH1F>("h_NoutElEl_AMS","",10,-3,3);
    h_NoutElMu_AMS = fs->make<TH1F>("h_NoutElMu_AMS","",10,-3,3);
    h_NinMuMu_ABS = fs->make<TH1F>("h_NinMuMu_ABS","",10,-3,3);
    h_NinElEl_ABS = fs->make<TH1F>("h_NinElEl_ABS","",10,-3,3);
    h_NinElMu_ABS = fs->make<TH1F>("h_NinElMu_ABS","",10,-3,3);
    h_NoutMuMu_ABS = fs->make<TH1F>("h_NoutMuMu_ABS","",10,-3,3);
    h_NoutElEl_ABS = fs->make<TH1F>("h_NoutElEl_ABS","",10,-3,3);
    h_NoutElMu_ABS = fs->make<TH1F>("h_NoutElMu_ABS","",10,-3,3);
    h_NinMuMu_NoMET_ABS = fs->make<TH1F>("h_NinMuMu_NoMET_ABS","",10,-3,3);
    h_NinElEl_NoMET_ABS = fs->make<TH1F>("h_NinElEl_NoMET_ABS","",10,-3,3);
    h_NinElMu_NoMET_ABS = fs->make<TH1F>("h_NinElMu_NoMET_ABS","",10,-3,3);
    h_NoutMuMu_NoMET_ABS = fs->make<TH1F>("h_NoutMuMu_NoMET_ABS","",10,-3,3);
    h_NoutElEl_NoMET_ABS = fs->make<TH1F>("h_NoutElEl_NoMET_ABS","",10,-3,3);
    h_NoutElMu_NoMET_ABS = fs->make<TH1F>("h_NoutElMu_NoMET_ABS","",10,-3,3);
    h_NinMuMu_AoneBS = fs->make<TH1F>("h_NinMuMu_AoneBS","",10,-3,3);
    h_NinElEl_AoneBS = fs->make<TH1F>("h_NinElEl_AoneBS","",10,-3,3);
    h_NinElMu_AoneBS = fs->make<TH1F>("h_NinElMu_AoneBS","",10,-3,3);
    h_NoutMuMu_AoneBS = fs->make<TH1F>("h_NoutMuMu_AoneBS","",10,-3,3);
    h_NoutElEl_AoneBS = fs->make<TH1F>("h_NoutElEl_AoneBS","",10,-3,3);
    h_NoutElMu_AoneBS = fs->make<TH1F>("h_NoutElMu_AoneBS","",10,-3,3);
    h_NinMuMu_NoMET_AoneBS = fs->make<TH1F>("h_NinMuMu_NoMET_AoneBS","",10,-3,3);
    h_NinElEl_NoMET_AoneBS = fs->make<TH1F>("h_NinElEl_NoMET_AoneBS","",10,-3,3);
    h_NinElMu_NoMET_AoneBS = fs->make<TH1F>("h_NinElMu_NoMET_AoneBS","",10,-3,3);
    h_NoutMuMu_NoMET_AoneBS = fs->make<TH1F>("h_NoutMuMu_NoMET_AoneBS","",10,-3,3);
    h_NoutElEl_NoMET_AoneBS = fs->make<TH1F>("h_NoutElEl_NoMET_AoneBS","",10,-3,3);
    h_NoutElMu_NoMET_AoneBS = fs->make<TH1F>("h_NoutElMu_NoMET_AoneBS","",10,-3,3);

    h_ATS_mLepNoVetoMuMu = fs->make<TH1F>("h_ATS_mLepNoVetoMuMu",";M_{ll};",300,0.,300.);
    h_ATS_mLepNoVetoElEl = fs->make<TH1F>("h_ATS_mLepNoVetoElEl",";M_{ll};",300,0.,300.);
    h_ATS_mLepNoVetoElMu = fs->make<TH1F>("h_ATS_mLepNoVetoElMu",";M_{ll};",300,0.,300.);
    h_ATS_mLepNoVetoDiLep = fs->make<TH1F>("h_ATS_mLepNoVetoDiLep",";M_{ll};",300,0.,300.);
    h_NinMuMu_ATS = fs->make<TH1F>("h_NinMuMu_ATS","",10,-3,3);
    h_NinElEl_ATS = fs->make<TH1F>("h_NinElEl_ATS","",10,-3,3);
    h_NinElMu_ATS = fs->make<TH1F>("h_NinElMu_ATS","",10,-3,3);
    h_NoutMuMu_ATS = fs->make<TH1F>("h_NoutMuMu_ATS","",10,-3,3);
    h_NoutElEl_ATS = fs->make<TH1F>("h_NoutElEl_ATS","",10,-3,3);
    h_NoutElMu_ATS = fs->make<TH1F>("h_NoutElMu_ATS","",10,-3,3);
    h_NinMuMu_NoMET_ATS = fs->make<TH1F>("h_NinMuMu_NoMET_ATS","",10,-3,3);
    h_NinElEl_NoMET_ATS = fs->make<TH1F>("h_NinElEl_NoMET_ATS","",10,-3,3);
    h_NinElMu_NoMET_ATS = fs->make<TH1F>("h_NinElMu_NoMET_ATS","",10,-3,3);
    h_NoutMuMu_NoMET_ATS = fs->make<TH1F>("h_NoutMuMu_NoMET_ATS","",10,-3,3);
    h_NoutElEl_NoMET_ATS = fs->make<TH1F>("h_NoutElEl_NoMET_ATS","",10,-3,3);
    h_NoutElMu_NoMET_ATS = fs->make<TH1F>("h_NoutElMu_NoMET_ATS","",10,-3,3);

    //initialize the tree
    f_outFile->cd();
    t_outTree =  new TTree("tree","tr");
    //    t_outTree->Branch("HLT_IsoMu24_v",&b_Mu1);
    //    t_outTree->Branch("HLT_IsoTkMu24_v",&b_Mu2);
    //    t_outTree->Branch("HLT_IsoMu22_eta2p1_v",&b_Mu3);
    //    t_outTree->Branch("HLT_IsoTkMu22_eta2p1_v",&b_Mu4);
    //    t_outTree->Branch("HLT_IsoTkMu22_eta2p1_v",&b_Mu3);
    //    t_outTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",&b_MuMu1);
    //    t_outTree->Branch("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",&b_MuMu2);
    //    t_outTree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",&b_MuMu3);
    //    t_outTree->Branch("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",&b_MuMu4);
    //    t_outTree->Branch("HLT_Ele32_eta2p1_WPTight_Gsf_v",&b_El1);
    //    t_outTree->Branch("HLT_Ele27_WPTight_Gsf_v",&b_El2);
    //    t_outTree->Branch("HLT_Ele25_eta2p1_WPTight_Gsf_v",&b_El3);
    //    t_outTree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&b_ElEl1);
    //    t_outTree->Branch("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v",&b_ElEl1);
    //    t_outTree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&b_MuMu1);
    //    t_outTree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&b_MuMu2);
    //    t_outTree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&b_MuMu3);
    //    t_outTree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&b_MuMu4);
    //    t_outTree->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&b_MuMu5);
    //    t_outTree->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&b_MuMu6);
    //    cout << "here I tree 0!" <<endl;

    // t_outTree->Branch("Leptons","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_Leptons);
    // t_outTree->Branch("bJets","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_bJets);
    //t_outTree->Branch("Jets","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_Jets);

    //    t_outTree->Branch("top","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_top);
    //   t_outTree->Branch("top_","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_antitop);
    //  t_outTree->Branch("W","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_W);
    // t_outTree->Branch("W_","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_antiW);
    //t_outTree->Branch("nu","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_nu);
    // t_outTree->Branch("nu_","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&t_antinu);
    t_outTree->Branch("cos","vector<double>",&t_cos);
    t_outTree->Branch("Run",&t_Run);
    t_outTree->Branch("Event",&t_Event);
    t_outTree->Branch("Lumi",&t_lumi);
    t_outTree->Branch("Bunch",&t_bunch);
    t_outTree->Branch("w_mc",&w_mc);
    t_outTree->Branch("w_trigger_ee",&w_ee);
    t_outTree->Branch("w_trigger_em",&w_em);
    t_outTree->Branch("w_trigger_mm",&w_mm);
    t_outTree->Branch("w_posEl",&w_posEl);
    t_outTree->Branch("w_negEl",&w_negEl);
    t_outTree->Branch("w_posMu",&w_posMu);
    t_outTree->Branch("w_negMu",&w_negMu);
    t_outTree->Branch("w_top",&w_top);
    t_outTree->Branch("w_bjets","vector<double>",&w_bjets);
    t_outTree->Branch("rho",&t_rho);
    t_outTree->Branch("PU_weight",&t_PU_weight);
    t_outTree->Branch("PU_weight_secondWay",&t_PU_weight_secondWay);

    t_outTree->Branch("pt_Leptons","vector<double>",&t_Leptons_pt);
    t_outTree->Branch("eta_Leptons","vector<double>",&t_Leptons_eta);
    t_outTree->Branch("phi_Leptons","vector<double>",&t_Leptons_phi);
    t_outTree->Branch("e_Leptons","vector<double>",&t_Leptons_e);
    t_outTree->Branch("charge_Leptons","vector<int>",&t_Leptons_charge);
    t_outTree->Branch("type_Leptons","vector<string>",&t_Leptons_type);
    t_outTree->Branch("SuperClusterEta_Leptons","vector<double>",&t_Leptons_etaSC);
    t_outTree->Branch("d0_Leptons","vector<double>",&t_Leptons_d0);
    t_outTree->Branch("dz_Leptons","vector<double>",&t_Leptons_dz);


    t_outTree->Branch("pt_gen_Leptons","vector<double>",&t_gen_Leptons_pt);
    t_outTree->Branch("eta_gen_Leptons","vector<double>",&t_gen_Leptons_eta);
    t_outTree->Branch("phi_gen_Leptons","vector<double>",&t_gen_Leptons_phi);
    t_outTree->Branch("e_gen_Leptons","vector<double>",&t_gen_Leptons_e);
    t_outTree->Branch("charge_gen_Leptons","vector<int>",&t_gen_Leptons_charge);
    t_outTree->Branch("id_gen_Leptons","vector<int>",&t_gen_Leptons_id);

    t_outTree->Branch("pt_gen_daughters","vector<double>",&t_gen_daughters_pt);
    t_outTree->Branch("eta_gen_daughters","vector<double>",&t_gen_daughters_eta);
    t_outTree->Branch("phi_gen_daughters","vector<double>",&t_gen_daughters_phi);
    t_outTree->Branch("e_gen_daughters","vector<double>",&t_gen_daughters_e);
    t_outTree->Branch("id_gen_daughters","vector<int>",&t_gen_daughters_id);

    t_outTree->Branch("pt_gen_mothers","vector<double>",&t_gen_mothers_pt);
    t_outTree->Branch("eta_gen_mothers","vector<double>",&t_gen_mothers_eta);
    t_outTree->Branch("phi_gen_mothers","vector<double>",&t_gen_mothers_phi);
    t_outTree->Branch("e_gen_mothers","vector<double>",&t_gen_mothers_e);
    t_outTree->Branch("id_gen_mothers","vector<int>",&t_gen_mothers_id);

    t_outTree->Branch("pt_Jets","vector<double>",&t_Jets_pt);
    t_outTree->Branch("eta_Jets","vector<double>",&t_Jets_eta);
    t_outTree->Branch("phi_Jets","vector<double>",&t_Jets_phi);
    t_outTree->Branch("e_Jets","vector<double>",&t_Jets_e);
    t_outTree->Branch("hadflav_Jets","vector<int>",&t_Jets_hadflav);
    t_outTree->Branch("ptresolution_Jets","vector<double>",&t_Jets_pt_resolution);
    t_outTree->Branch("phiresolution_Jets","vector<double>",&t_Jets_phi_resolution);
    t_outTree->Branch("sf_nominal_Jets","vector<double>",&t_Jets_sf_nominal);
    t_outTree->Branch("sf_up_Jets","vector<double>",&t_Jets_sf_up);
    t_outTree->Branch("sf_down_Jets","vector<double>",&t_Jets_sf_down);
    t_outTree->Branch("bdis_Jets","vector<double>",&t_Jets_dis);


    t_outTree->Branch("Loose_Bdiscriminator","vector<bool>",&t_Jets_LBdiscriminator);
    t_outTree->Branch("Medium_Bdiscriminator","vector<bool>",&t_Jets_MBdiscriminator);
    t_outTree->Branch("Tight_Bdiscriminator","vector<bool>",&t_Jets_TBdiscriminator);

    t_outTree->Branch("pt_gen_Jets","vector<double>",&t_gen_Jets_pt);
    t_outTree->Branch("eta_gen_Jets","vector<double>",&t_gen_Jets_eta);
    t_outTree->Branch("phi_gen_Jets","vector<double>",&t_gen_Jets_phi);
    t_outTree->Branch("e_gen_Jets","vector<double>",&t_gen_Jets_e);
    t_outTree->Branch("hadFlavour_gen_Jets","vector<double>",&t_gen_Jets_hadflav);
    t_outTree->Branch("partonFlavour_gen_Jets","vector<double>",&t_gen_Jets_prtnflav);

    t_outTree->Branch("pt_bJets","vector<double>",&t_bJets_pt);
    t_outTree->Branch("eta_bJets","vector<double>",&t_bJets_eta);
    t_outTree->Branch("phi_bJets","vector<double>",&t_bJets_phi);
    t_outTree->Branch("e_bJets","vector<double>",&t_bJets_e);

    t_outTree->Branch("pt_mets","vector<double>",&t_Met_pt);
    t_outTree->Branch("px_mets","vector<double>",&t_Met_px);
    t_outTree->Branch("py_mets","vector<double>",&t_Met_py);
    t_outTree->Branch("sum_mets","vector<double>",&t_Met_sum);



    t_outTree->Branch("num_PU_vertices",&t_num_PU_vertices);
    t_outTree->Branch("PU_BunchCrossing",&t_PU_BunchCrossing);
    t_outTree->Branch("num_PU_gen_vertices",&t_num_PU_gen_vertices);
    t_outTree->Branch("num_PV",&t_num_PV);
    t_outTree->Branch("size_Leptons",&t_Leptons_size);
    t_outTree->Branch("size_Electrons",&t_Electrons_size);
    t_outTree->Branch("size_Muons",&t_Muons_size);
    
    t_outTree->Branch("pt_top","vector<double>",&t_top_pt);
    t_outTree->Branch("eta_top","vector<double>",&t_top_eta);
    t_outTree->Branch("phi_top","vector<double>",&t_top_phi);
    t_outTree->Branch("e_top","vector<double>",&t_top_e);
    
     t_outTree->Branch("pt_antitop","vector<double>",&t_antitop_pt);
    t_outTree->Branch("eta_antitop","vector<double>",&t_antitop_eta);
    t_outTree->Branch("phi_antitop","vector<double>",&t_antitop_phi);
    t_outTree->Branch("e_antitop","vector<double>",&t_antitop_e);
    
     t_outTree->Branch("pt_W","vector<double>",&t_W_pt);
    t_outTree->Branch("eta_W","vector<double>",&t_W_eta);
    t_outTree->Branch("phi_W","vector<double>",&t_W_phi);
    t_outTree->Branch("e_W","vector<double>",&t_W_e);
    
     t_outTree->Branch("pt_antiW","vector<double>",&t_antiW_pt);
    t_outTree->Branch("eta_antiW","vector<double>",&t_antiW_eta);
    t_outTree->Branch("phi_antiW","vector<double>",&t_antiW_phi);
    t_outTree->Branch("e_antiW","vector<double>",&t_antiW_e);
    
     t_outTree->Branch("pt_nu","vector<double>",&t_nu_pt);
    t_outTree->Branch("eta_nu","vector<double>",&t_nu_eta);
    t_outTree->Branch("phi_nu","vector<double>",&t_nu_phi);
    t_outTree->Branch("e_nu","vector<double>",&t_nu_e);
    
     t_outTree->Branch("pt_antinu","vector<double>",&t_antinu_pt);
    t_outTree->Branch("eta_antinu","vector<double>",&t_antinu_eta);
    t_outTree->Branch("phi_antinu","vector<double>",&t_antinu_phi);
    t_outTree->Branch("e_antinu","vector<double>",&t_antinu_e);
    //    cout << "here I tree 1!" <<endl;



}
void MiniAnalyzer::endJob()
{
    //    cout << "here I start1!" <<endl;

    f_outFile->cd();
    h_Nevents->Write();
    h_NPV->Write();
    t_outTree->Write();
    f_outFile->Close();
    amwtSolver->writeOut();
    //    cout << "here I start2!" <<endl;

}


//____________________________________________________________________________
bool MiniAnalyzer::isMediumMuon(const reco::Muon & recoMu)
{
    bool goodGlob = recoMu.isGlobalMuon() &&
            recoMu.globalTrack()->normalizedChi2() < 3 &&
            recoMu.combinedQuality().chi2LocalPosition < 12 &&
            recoMu.combinedQuality().trkKink < 20;
    bool isMedium = muon::isLooseMuon(recoMu) &&
            recoMu.innerTrack()->validFraction() > 0.49 &&
            muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451);
    return isMedium;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//____________________________________________________________________________
void MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
void MiniAnalyzer::beginRun(const Run & iRun, const EventSetup &)
{


}

void MiniAnalyzer::endRun(const Run & iRun, const EventSetup &)
{
    //cout << "Im in End Run" <<endl;
    if(isData)return;
    if(isPythia) return;
    //cout << "isPythia: "<<isPythia << endl;
    //cout << "Why here" <<endl;

    edm::Handle<LHERunInfoProduct> run;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    iRun.getByToken(lheInfo_,run);
    //    iRun.getByLabel( "externalLHEProducer", run );
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());

    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++)
    {
        std::cout << iter->tag() << std::endl;
        std::vector<std::string> lines = iter->lines();
        for (unsigned int iLine = 0; iLine<lines.size(); iLine++)
        {
            std::cout << lines.at(iLine);
        }
    }

}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
