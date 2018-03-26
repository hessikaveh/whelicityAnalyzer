import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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
readFiles.extend( ['file:/cmsdata3/kaveh/workarea/16E7A136-1FBE-E611-B8DA-0025905C3DD0.root'
]);
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


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
    ,outputCommands = cms.untracked.vstring('drop *',
      "keep *_offlineSlimmedPrimaryVertices_*_*",
      "keep *_slimmedMuons_*_*",
       "keep *_slimmedElectrons_*_*",
       "keep *_slimmedJets_*_*",
       "keep *_slimmedMETs_*_*",
       "keep *_packedGenParticles_*_*",
       "keep *_prunedGenParticles_*_*",
       "keep *_packedGenParticles_*_*",
       "keep *_generator_*_*",
       "keep *_externalLHEProducer_*_*",
       "keep *_fixedGridRhoFastjetAll_*_*",
       "keep *_reducedEgamma_*_*",
       "keep *_TriggerResults_*_*")

)





# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.egmGsfElectronIDSequence * process.Whelicity)

