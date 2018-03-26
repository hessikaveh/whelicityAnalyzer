import FWCore.ParameterSet.Config as cms



process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles)
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
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles)

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles)

readFiles.extend( [
'root://xrootd-cms.infn.it//store/data/Run2016E/DoubleEG/MINIAOD/03Feb2017-v1/110000/EA7C2D56-A1EA-E611-86B2-0CC47A13CC7A.root']);
# Bad Charged Hadron and Bad Muon Filters from MiniAOD
#process.load('RecoMET.METFilters.metFilters_cff')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
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

process.load("whelicity1.MiniAnalyzer.whelicity_cff")
process.Whelicity.isData = cms.bool(True)
process.Whelicity.isPythia = cms.bool(False)
process.Whelicity.isRunH = cms.bool(True)
process.Whelicity.DiEl = cms.bool(True)
process.Whelicity.outFileName = cms.string("tree.root")
process.TFileService = cms.Service("TFileService",
 fileName = cms.string("histos.root")
)


# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.BadPFMuonFilter           *
					 process.BadChargedCandidateFilter * 
					 process.egmGsfElectronIDSequence  * 
					 process.Whelicity                   )




