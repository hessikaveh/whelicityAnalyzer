#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs

#Data DoubleMuon Single
crab submit -c test/crab_Data_DiMu_11.py General.requestName=SingleMuonMuMuRunEJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_9.py General.requestName=SingleMuonMuMuRunCJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_10.py General.requestName=SingleMuonMuMuRunDJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_12.py General.requestName=SingleMuonMuMuRunFJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_13.py General.requestName=SingleMuonMuMuRunGJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_8.py General.requestName=SingleMuonMuMuRunBJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_14.py General.requestName=SingleMuonMuMuRunH1Jan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_15.py General.requestName=SingleMuonMuMuRunH2Jan13 Data.outputDatasetTag=Jan13MuMuDATA
#Data DoubleMuon
crab submit -c test/crab_Data_DiMu_0.py General.requestName=DoubleMuonRunBJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_1.py General.requestName=DoubleMuonRunCJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_2.py General.requestName=DoubleMuonRunDJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_3.py General.requestName=DoubleMuonRunEJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_4.py General.requestName=DoubleMuonRunFJan13 Data.outputDatasetTag=Jan13MuMuDATA

crab submit -c test/crab_Data_DiMu_5.py General.requestName=DoubleMuonRunGJan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_6.py General.requestName=DoubleMuonRunH1Jan13 Data.outputDatasetTag=Jan13MuMuDATA
crab submit -c test/crab_Data_DiMu_7.py General.requestName=DoubleMuonRunH2Jan13 Data.outputDatasetTag=Jan13MuMuDATA

#TT
crab submit -c test/crab_MC_DiMu_0.py General.requestName=TTJan13MuMu Data.outputDatasetTag=TTJan13MuMuMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiMu_1.py General.requestName=tbarWJan13MuMu Data.outputDatasetTag=tbarWJan13MuMuMC
crab submit -c test/crab_MC_DiMu_2.py General.requestName=tWJan13MuMu Data.outputDatasetTag=tWJan13MuMuMC

#WW
crab submit -c test/crab_MC_DiMu_3.py General.requestName=WWJan13MuMu Data.outputDatasetTag=WWJan13MuMuMC
#WZ
crab submit -c test/crab_MC_DiMu_4.py General.requestName=WZJan13MuMu Data.outputDatasetTag=WZJan13MuMuMC
#ZZ
crab submit -c test/crab_MC_DiMu_5.py General.requestName=ZZJan13MuMu Data.outputDatasetTag=ZZJan13MuMuMC
#W+Jets
crab submit -c test/crab_MC_DiMu_6.py General.requestName=WJetsJan13MuMu Data.outputDatasetTag=WJetsJan13MuMuMC
#DY
crab submit -c test/crab_MC_DiMu_7.py General.requestName=DY10Jan13MuMu Data.outputDatasetTag=DY10Jan13MuMuMC
crab submit -c test/crab_MC_DiMu_8.py General.requestName=DY50Jan13MuMu Data.outputDatasetTag=DY50Jan13MuMuMC

#TTWJets
crab submit -c test/crab_MC_DiMu_9.py General.requestName=TTWJetsToLNuJan13_2MuMu Data.outputDatasetTag=TTWJetsToLNuJan13_2MuMuMC
crab submit -c test/crab_MC_DiMu_10.py General.requestName=TTWJetsToQQJan13MuMu Data.outputDatasetTag=TTWJetsToQQJan13MuMuMC
#TTZ
crab submit -c test/crab_MC_DiMu_12.py General.requestName=TTZToLLNuNuJan13_2MuMu Data.outputDatasetTag=TTZToLLNuNuJan13_2MuMuMC
crab submit -c test/crab_MC_DiMu_11.py General.requestName=TTZToQQJan13MuMu Data.outputDatasetTag=TTZToQQJan13MuMuMC


