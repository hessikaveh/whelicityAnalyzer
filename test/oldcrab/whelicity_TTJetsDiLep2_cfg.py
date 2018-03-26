import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
## configure process options
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True)
)


##____________________________________________________________________________||
process.noscraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )

##____________________________________________________________________________||
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)

##____________________________________________________________________________||
#process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

##____________________________________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

##____________________________________________________________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.debug = cms.bool(True)

##____________________________________________________________________________||
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')

process.load('RecoMET.METFilters.jetIDFailureFilter_cfi')
process.jetIDFailure.MinJetPt  = cms.double(30.0)
process.jetIDFailure.MaxJetEta = cms.double(999.0)

process.noscraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )



## configure geometry & conditions
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles)

process.load("whelicity1.MiniAnalyzer.whelicity_cff")
process.Whelicity.isData = cms.bool(False)

process.Whelicity.outFileName = cms.string("tree.root")
process.TFileService = cms.Service("TFileService",
 fileName = cms.string("histos.root")
)
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJets')
#process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3Residual')

process.load('RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi')

process.load('RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi')
process.load('RecoMET.METFilters.eeNoiseFilter_cfi')

process.load('RecoMET/METAnalyzers/CSCHaloFilter_cfi')
process.rejectRecov = cms.EDFilter(
  "RecovRecHitFilter",
  EERecHitSource = cms.InputTag("reducedEcalRecHitsEE"),
  MinRecovE = cms.double(30),
  TaggingMode=cms.bool(False)
)

process.Vertex = cms.Path(process.goodOfflinePrimaryVertices*getattr(process,"patPF2PATSequence"+postfix)*~process.primaryVertexFilter)
process.Scraping = cms.Path(~process.noscraping)
process.HBHENoise = cms.Path(~process.HBHENoiseFilter)
process.CSCTightHalo = cms.Path(~process.CSCTightHaloFilter)
process.RecovRecHit = cms.Path(~process.rejectRecov)

process.HCALLaser = cms.Path(~process.hcalLaserEventFilter)
process.ECALDeadCellTP = cms.Path(~process.EcalDeadCellTriggerPrimitiveFilter)
process.ECALDeadCellBE = cms.Path(~process.EcalDeadCellBoundaryEnergyFilter)
process.jetID = cms.Path(~process.jetIDFailure)
process.trackingFailure = cms.Path(process.goodVertices*~process.trackingFailureFilter)
process.inconsistentMuon = cms.Path(~process.inconsistentMuonPFCandidateFilter)
process.greedyMuon = cms.Path(~process.greedyMuonPFCandidateFilter)
process.eeNoise = cms.Path(~process.eeNoiseFilter)
#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
useAOD = False
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.egmGsfElectronIDSequence * process.Whelicity)


