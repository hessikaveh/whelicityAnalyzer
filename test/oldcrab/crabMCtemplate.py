from CRABClient.UserUtilities import config, getUsernameFromSiteDB
#import FWCore.ParameterSet.VarParsing as VarParsing

# setup 'analysis'  options
#options = VarParsing.VarParsing ('analysis')
#options.register ('name',
#                  'MCtemplate', 
#                  VarParsing.VarParsing.multiplicity.singleton, 
#                  VarParsing.VarParsing.varType.string,          
#                  "request name")

#options.inputFiles= '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'


#options.parseArguments()

config = config()


config.General.requestName   = 'MCTemplate'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'test/whelicity_template_cfg.py'
config.JobType.inputFiles  = ['Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt','Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt','Spring16_25nsV10_MC_SF_AK4PFchs.txt']
config.JobType.outputFiles = ['tree.root','sanityCheckHistos.root']

config.Data.inputDataset = 'DataSet'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality  = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = True
config.Data.outputDatasetTag = 'MC'

config.Site.storageSite = 'T3_IR_IPM'
