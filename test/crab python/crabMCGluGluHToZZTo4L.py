from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName   = 'STtW27Feb'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/whelicity_ST_cfg.py'
config.JobType.inputFiles  = ['Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt','Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt','Spring16_25nsV10_MC_SF_AK4PFchs.txt']
config.JobType.outputFiles = ['tree.root','sanityCheckHistos.root']
config.Data.inputDataset = '/GluGluHToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/RunIISpring15MiniAODv2-AsymptFlat10to50bx25Raw_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality  = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = True
config.Data.outputDatasetTag = 'MC'

config.Site.storageSite = 'T3_IR_IPM'
