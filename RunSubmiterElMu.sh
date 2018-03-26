#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs

#Data MuonEG Single
crab submit -c test/crab_Data_ElMu_11.py General.requestName=SingleElectronElMuRunEJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_9.py General.requestName=SingleElectronElMuRunCJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_10.py General.requestName=SingleElectronElMuRunDJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_12.py General.requestName=SingleElectronElMuRunFJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_13.py General.requestName=SingleElectronElMuRunGJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_8.py General.requestName=SingleElectronElMuRunBJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_14.py General.requestName=SingleElectronElMuRunH1Jan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_15.py General.requestName=SingleElectronElMuRunH2Jan07 Data.outputDatasetTag=Jan07ElMuDATA

crab submit -c test/crab_Data_ElMu_16.py General.requestName=SingleMuonElMuRunBJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_17.py General.requestName=SingleMuonElMuRunCJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_18.py General.requestName=SingleMuonElMuRunDJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_19.py General.requestName=SingleMuonElMuRunEJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_20.py General.requestName=SingleMuonElMuRunFJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_21.py General.requestName=SingleMuonElMuRunGJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_22.py General.requestName=SingleMuonElMuRunH1Jan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_23.py General.requestName=SingleMuonElMuRunH2Jan07 Data.outputDatasetTag=Jan07ElMuDATA

#Data MuonEG
crab submit -c test/crab_Data_ElMu_0.py General.requestName=MuonEGRunBJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_1.py General.requestName=MuonEGRunCJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_2.py General.requestName=MuonEGRunDJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_3.py General.requestName=MuonEGRunEJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_4.py General.requestName=MuonEGRunFJan07 Data.outputDatasetTag=Jan07ElMuDATA

crab submit -c test/crab_Data_ElMu_5.py General.requestName=MuonEGRunGJan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_6.py General.requestName=MuonEGRunH1Jan07 Data.outputDatasetTag=Jan07ElMuDATA
crab submit -c test/crab_Data_ElMu_7.py General.requestName=MuonEGRunH2Jan07 Data.outputDatasetTag=Jan07ElMuDATA

#TT
crab submit -c test/crab_MC_ElMu_0.py General.requestName=TTJan07_3ElMu Data.outputDatasetTag=TTJan07_3ElMuMC
#ST#tW tbarW
crab submit -c test/crab_MC_ElMu_1.py General.requestName=tbarWJan07_2ElMu Data.outputDatasetTag=tbarWJan07_2ElMuMC
crab submit -c test/crab_MC_ElMu_2.py General.requestName=tWJan07ElMu Data.outputDatasetTag=tWJan07ElMuMC

#WW
crab submit -c test/crab_MC_ElMu_3.py General.requestName=WWJan07ElMu Data.outputDatasetTag=WWJan07ElMuMC
#WZ
crab submit -c test/crab_MC_ElMu_4.py General.requestName=WZJan07ElMu Data.outputDatasetTag=WZJan07ElMuMC
#ZZ
crab submit -c test/crab_MC_ElMu_5.py General.requestName=ZZJan07ElMu Data.outputDatasetTag=ZZJan07ElMuMC
#W+Jets
crab submit -c test/crab_MC_ElMu_6.py General.requestName=WJetsJan07ElMu Data.outputDatasetTag=WJetsJan07ElMuMC
#DY
crab submit -c test/crab_MC_ElMu_7.py General.requestName=DY10Jan07ElMu Data.outputDatasetTag=DY10Jan07ElMuMC
crab submit -c test/crab_MC_ElMu_8.py General.requestName=DY50Jan07ElMu Data.outputDatasetTag=DY50Jan07ElMuMC

#TTWJets
crab submit -c test/crab_MC_ElMu_9.py General.requestName=TTWJetsToLNuJan07_3ElMu Data.outputDatasetTag=TTWJetsToLNuJan07_3ElMuMC
crab submit -c test/crab_MC_ElMu_10.py General.requestName=TTWJetsToQQJan07_2ElMu Data.outputDatasetTag=TTWJetsToQQJan07_2ElMuMC
#TTZ
crab submit -c test/crab_MC_ElMu_12.py General.requestName=TTZToLLNuNuJan07_3ElMu Data.outputDatasetTag=TTZToLLNuNuJan07_3ElMuMC
crab submit -c test/crab_MC_ElMu_11.py General.requestName=TTZToQQJan07_2ElMu Data.outputDatasetTag=TTZToQQJan07_2ElMuMC



