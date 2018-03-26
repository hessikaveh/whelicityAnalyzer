import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(75000) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
## configure process options
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True)
)
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
## configure geometry & conditions
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( [
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/02F0F9B7-EFBA-E611-9AE6-842B2B76670F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/123E406C-ECBA-E611-99CF-0CC47A7EEE80.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/209A94FB-E4BA-E611-842B-D4AE526A0461.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/36D6C860-26BB-E611-8BA1-001E674FB207.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3C36B95F-26BB-E611-A9BB-0CC47AB0B704.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/406146E0-E4BA-E611-BD30-002590D9D88C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/40BC2BA4-26BB-E611-A441-001E67DDC2CC.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/685DA755-26BB-E611-A3B1-5065F381B271.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/709F00A3-AEBA-E611-88E4-D4AE526A0B03.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/70DD10EC-B2BA-E611-B84C-842B2B7682B8.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8065AC88-B5BA-E611-917A-D4AE526A109A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/9C8926F5-B1BA-E611-A10B-0CC47A57D164.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B2692000-E9BA-E611-B465-0CC47A0AD3BC.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BC6EC879-26BB-E611-8A26-0025905A612C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CCBA6A52-EBBA-E611-89C5-001EC9AF1FBF.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D68DEE5B-26BB-E611-B0BC-90B11C267182.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/DA4FBC3E-BABA-E611-9A96-0CC47A7EEDB0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/DE900E81-26BB-E611-807F-0CC47A7EEC70.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/EE877568-26BB-E611-BA63-A4BADB1E6031.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F69FBE86-26BB-E611-8C23-0CC47A78A3EC.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FA06C7A1-EDBA-E611-A02C-0CC47A7EEE76.root' ] );


secFiles.extend( [
               ] )
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

#genParticleCollection = 'prunedGenParticles'
#process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
#process.initSubset.src = genParticleCollection
#process.decaySubset.src = genParticleCollection
#process.decaySubset.runMode = "Run2"
process.load("whelicity1.MiniAnalyzer.whelicity_cff")
process.Whelicity.isData = cms.bool(False)
process.Whelicity.outFileName = cms.string("TTtree4March.root")
process.TFileService = cms.Service("TFileService",
 fileName = cms.string("TTJetsDiLep4March.root")
)


# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.egmGsfElectronIDSequence * process.Whelicity)


