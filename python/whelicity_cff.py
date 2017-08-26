import FWCore.ParameterSet.Config as cms


Whelicity = cms.EDAnalyzer("MiniAnalyzer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
	packed = cms.InputTag("packedGenParticles"),
	pruned = cms.InputTag("prunedGenParticles"),
	genEvtInfo = cms.InputTag("generator"),
	lheEvtInfo = cms.InputTag("externalLHEProducer"),
	rho = cms.InputTag("fixedGridRhoFastjetAll"),
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt"),
	ttgen = cms.InputTag("genEvt"),
	isData = cms.bool(False),
	ptRes = cms.string("Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"),
	phiRes = cms.string("Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt"),
	sfRes = cms.string("Spring16_25nsV10_MC_SF_AK4PFchs.txt"),
	outFileName = cms.string("outFile.root"),
	triggerResults = cms.InputTag("TriggerResults","","HLT"),
	externalLHEProducer = cms.InputTag("externalLHEProducer"),
	isPythia = cms.bool(False),
	triggerFilters = cms.InputTag("TriggerResults","","PAT"),
        DiMu = cms.bool(False),
        DiEl = cms.bool(False),
        ElMu = cms.bool(False),
		egammaSF = cms.string("egammaEffi.txt_EGM2D.root"),
		isRunGH = cms.bool(False),
		muonISOSF = cms.string("ISOEfficienciesAndSF_BCDEF.root"),
		muonIDSF = cms.string("IDEfficienciesAndSF_BCDEF.root"),
		egammaTkSF = cms.string("TkegammaEffi.txt_EGM2D.root"),
		btagSf = cms.string("CSVv2_Moriond17_B_H.csv"),
		muonTkSF = cms.string("Tracking_EfficienciesAndSF_BCDEFGH.root")

	
)


