from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName   = 'DYJetsToLL'
#config.General.requestName   = 'TTJetsDiLep4'
#config.General.requestName   = 'DataMuEG'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/whelicity_DYJetsToLL_cfg.py'
#config.JobType.psetName = 'test/whelicity_TTJetsDiLep2_cfg.py'
config.JobType.inputFiles  = ['Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt','Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt']
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-herwigpp_30M/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'

#config.Data.inputDataset = '/MuonEG/Run2016B-23Sep2016-v3/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.ignoreLocality  = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = True
config.Data.outputDatasetTag = 'MC'

config.Site.storageSite = 'T3_IR_IPM'
