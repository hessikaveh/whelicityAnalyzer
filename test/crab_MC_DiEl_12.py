from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.General.requestName   = 'crabMC'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True 
config.General.transferLogs = False 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'test/whelicity_madgraph_DiEl_cfg.py' 
config.JobType.inputFiles  = ['Spring16_25nsV6_DATA_SF_AK4PFchs.txt','Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt','Spring16_25nsV6_DATA_PhiResolution_AK4PFchs.txt','pileup_distribution_moriond17.root','pileup_distribution_data16.root','Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt','Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt','Spring16_25nsV10_MC_SF_AK4PFchs.txt','egammaEffi.txt_EGM2D.root','ISOEfficienciesAndSF_BCDEF.root','IDEfficienciesAndSF_BCDEF.root','TkegammaEffi.txt_EGM2D.root','Tracking_EfficienciesAndSF_BCDEFGH.root','CSVv2_Moriond17_B_H.csv','ISOEfficienciesAndSF_GH.root','IDEfficienciesAndSF_GH.root','triggerSummary_ee.root','triggerSummary_emu.root','triggerSummary_mumu.root'] 
config.JobType.outputFiles = ['tree.root','sanityCheckHistos.root']
config.Data.inputDataset = '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outputDatasetTag = 'MC'
config.Site.blacklist = ['T2_US_*']
config.Site.storageSite = 'T3_IR_IPM'
