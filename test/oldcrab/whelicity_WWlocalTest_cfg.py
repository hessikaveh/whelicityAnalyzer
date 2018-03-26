import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(81000) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
## configure process options
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True)
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
readFiles.extend( [
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/08E155A9-FAB6-E611-92BF-00259073E45E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0CB16E4E-F0B6-E611-9D13-0090FAA58294.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/1AA5E3E6-F4B6-E611-8A54-0090FAA575E0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/265C6533-CDB6-E611-9E71-0090FAA58B94.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/287ADB8D-E1B6-E611-8332-00259073E4EA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/309E843E-CFB6-E611-9E6F-00259073E496.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/34023683-C1B6-E611-ACA4-00259074AE94.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/3A344242-D7B6-E611-8A3D-0CC47A4DEE54.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/3AF4E14E-C5B6-E611-9A6B-00259074AE98.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/4462CF47-BCB6-E611-A504-00259073E482.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/50B4039E-DDB6-E611-89EC-00259073E488.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/72CD4EDC-C8B6-E611-A1EC-002590D0B0BE.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/76742A18-D9B6-E611-A98F-00259073E384.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/82EEAB2E-E3B6-E611-A906-002590D0AFC2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/84BD3443-EDB6-E611-A863-20CF305B0501.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/8894C5DA-D5B6-E611-BE17-0090FAA587C4.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/968C1EFA-DEB6-E611-91AF-20CF30561726.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A2BFEAB5-EAB6-E611-8842-20CF305B04D2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/AE805894-E5B6-E611-B84F-00259073E4E2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B415524C-CBB6-E611-AD28-00259073E496.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B880FA55-FFB6-E611-BBD1-002590D0B056.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/BA714726-D1B6-E611-B6CF-00259073E3D2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/BC259EE6-DBB6-E611-91DB-002590747DD8.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DE039D90-D3B6-E611-BBD4-002590747E1C.root' ] );


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


