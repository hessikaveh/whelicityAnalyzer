from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName   = 'DataDoubleEG01'
config.General.requestName   = 'DataMuEG25Feb17'
#config.General.requestName   = 'DataDoubleMu01'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/whelicity_DataElMu_SM_H_GH_cfg.py'
config.JobType.inputFiles  = ['Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt','Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt','Spring16_25nsV10_MC_SF_AK4PFchs.txt','egammaEffi.txt_EGM2D.root','ISOEfficienciesAndSF_GH.root','IDEfficienciesAndSF_GH.root','TkegammaEffi.txt_EGM2D.root','Tracking_EfficienciesAndSF_BCDEFGH.root','CSVv2_Moriond17_B_H.csv']
config.JobType.outputFiles = ['tree.root','sanityCheckHistos.root']
#config.Data.inputDataset = '/DoubleEG/Run2016B-23Sep2016-v3/MINIAOD'
config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2016B-23Sep2016-v3/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.ignoreLocality  = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = True
config.Data.outputDatasetTag = 'March'

config.Site.storageSite = 'T3_IR_IPM'
