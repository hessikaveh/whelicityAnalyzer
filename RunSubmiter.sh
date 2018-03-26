#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs
#Data DoubleEG Single
crab submit -c test/crab_Data_DiEl_11.py General.requestName=SingleElectronElElRunENov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_9.py General.requestName=SingleElectronElElRunCNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_10.py General.requestName=SingleElectronElElRunDNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_12.py General.requestName=SingleElectronElElRunFNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_13.py General.requestName=SingleElectronElElRunGNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_8.py General.requestName=SingleElectronElElRunBNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_14.py General.requestName=SingleElectronElElRunH1Nov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_15.py General.requestName=SingleElectronElElRunH2Nov15 Data.outputDatasetTag=Nov15ElElDATA
#Data DoubleEG
crab submit -c test/crab_Data_DiEl_0.py General.requestName=DoubleEGRunBNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_1.py General.requestName=DoubleEGRunCNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_2.py General.requestName=DoubleEGRunDNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_3.py General.requestName=DoubleEGRunENov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_4.py General.requestName=DoubleEGRunFNov15 Data.outputDatasetTag=Nov15ElElDATA

crab submit -c test/crab_Data_DiEl_5.py General.requestName=DoubleEGRunGNov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_6.py General.requestName=DoubleEGRunH1Nov15 Data.outputDatasetTag=Nov15ElElDATA
crab submit -c test/crab_Data_DiEl_7.py General.requestName=DoubleEGRunH2Nov15 Data.outputDatasetTag=Nov15ElElDATA

#TT
crab submit -c test/crab_MC_DiEl_0.py General.requestName=TTNov15ElEl Data.outputDatasetTag=TTNov15ElElMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiEl_1.py General.requestName=tbarWNov15ElEl Data.outputDatasetTag=tbarWNov15ElElMC
crab submit -c test/crab_MC_DiEl_2.py General.requestName=tWNov15ElEl Data.outputDatasetTag=tWNov15ElElMC

#WW
crab submit -c test/crab_MC_DiEl_3.py General.requestName=WWNov15ElEl Data.outputDatasetTag=WWNov15ElElMC
#WZ
crab submit -c test/crab_MC_DiEl_4.py General.requestName=WZNov15ElEl Data.outputDatasetTag=WZNov15ElElMC
#ZZ
crab submit -c test/crab_MC_DiEl_5.py General.requestName=ZZNov15ElEl Data.outputDatasetTag=ZZNov15ElElMC
#W+Jets
crab submit -c test/crab_MC_DiEl_6.py General.requestName=WJetsNov15ElEl Data.outputDatasetTag=WJetsNov15ElElMC
#DY
crab submit -c test/crab_MC_DiEl_7.py General.requestName=DY10Nov15ElEl Data.outputDatasetTag=DY10Nov15ElElMC
crab submit -c test/crab_MC_DiEl_8.py General.requestName=DY50Nov15ElEl Data.outputDatasetTag=DY50Nov15ElElMC

#TTWJets
crab submit -c test/crab_MC_DiEl_9.py General.requestName=TTWJetsToLNuNov15ElEl Data.outputDatasetTag=TTWJetsToLNuNov15ElElMC
crab submit -c test/crab_MC_DiEl_10.py General.requestName=TTWJetsToQQNov15ElEl Data.outputDatasetTag=TTWJetsToQQNov15ElElMC
#TTZ
crab submit -c test/crab_MC_DiEl_12.py General.requestName=TTZToLLNuNuNov15ElEl Data.outputDatasetTag=TTZToLLNuNuNov15ElElMC
crab submit -c test/crab_MC_DiEl_11.py General.requestName=TTZToQQNov15ElEl Data.outputDatasetTag=TTZToQQNov15ElElMC





#--------------------------------------
#Data DoubleMuon Single
crab submit -c test/crab_Data_DiMu_11.py General.requestName=SingleMuonMuMuRunENov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_9.py General.requestName=SingleMuonMuMuRunCNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_10.py General.requestName=SingleMuonMuMuRunDNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_12.py General.requestName=SingleMuonMuMuRunFNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_13.py General.requestName=SingleMuonMuMuRunGNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_8.py General.requestName=SingleMuonMuMuRunBNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_14.py General.requestName=SingleMuonMuMuRunH1Nov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_15.py General.requestName=SingleMuonMuMuRunH2Nov15 Data.outputDatasetTag=Nov15MuMuDATA
#Data DoubleMuon
crab submit -c test/crab_Data_DiMu_0.py General.requestName=DoubleMuonRunBNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_1.py General.requestName=DoubleMuonRunCNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_2.py General.requestName=DoubleMuonRunDNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_3.py General.requestName=DoubleMuonRunENov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_4.py General.requestName=DoubleMuonRunFNov15 Data.outputDatasetTag=Nov15MuMuDATA

crab submit -c test/crab_Data_DiMu_5.py General.requestName=DoubleMuonRunGNov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_6.py General.requestName=DoubleMuonRunH1Nov15 Data.outputDatasetTag=Nov15MuMuDATA
crab submit -c test/crab_Data_DiMu_7.py General.requestName=DoubleMuonRunH2Nov15 Data.outputDatasetTag=Nov15MuMuDATA

#TT
crab submit -c test/crab_MC_DiMu_0.py General.requestName=TTNov15MuMu Data.outputDatasetTag=TTNov15MuMuMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiMu_1.py General.requestName=tbarWNov15MuMu Data.outputDatasetTag=tbarWNov15MuMuMC
crab submit -c test/crab_MC_DiMu_2.py General.requestName=tWNov15MuMu Data.outputDatasetTag=tWNov15MuMuMC

#WW
crab submit -c test/crab_MC_DiMu_3.py General.requestName=WWNov15MuMu Data.outputDatasetTag=WWNov15MuMuMC
#WZ
crab submit -c test/crab_MC_DiMu_4.py General.requestName=WZNov15MuMu Data.outputDatasetTag=WZNov15MuMuMC
#ZZ
crab submit -c test/crab_MC_DiMu_5.py General.requestName=ZZNov15MuMu Data.outputDatasetTag=ZZNov15MuMuMC
#W+Jets
crab submit -c test/crab_MC_DiMu_6.py General.requestName=WJetsNov15MuMu Data.outputDatasetTag=WJetsNov15MuMuMC
#DY
crab submit -c test/crab_MC_DiMu_7.py General.requestName=DY10Nov15MuMu Data.outputDatasetTag=DY10Nov15MuMuMC
crab submit -c test/crab_MC_DiMu_8.py General.requestName=DY50Nov15MuMu Data.outputDatasetTag=DY50Nov15MuMuMC

#TTWJets
crab submit -c test/crab_MC_DiMu_9.py General.requestName=TTWJetsToLNuNov15MuMu Data.outputDatasetTag=TTWJetsToLNuNov15MuMuMC
crab submit -c test/crab_MC_DiMu_10.py General.requestName=TTWJetsToQQNov15MuMu Data.outputDatasetTag=TTWJetsToQQNov15MuMuMC
#TTZ
crab submit -c test/crab_MC_DiMu_12.py General.requestName=TTZToLLNuNuNov15MuMu Data.outputDatasetTag=TTZToLLNuNuNov15MuMuMC
crab submit -c test/crab_MC_DiMu_11.py General.requestName=TTZToQQNov15MuMu Data.outputDatasetTag=TTZToQQNov15MuMuMC



#-------------------------------------------------------------------------------------------------
#Data MuonEG Single
crab submit -c test/crab_Data_ElMu_11.py General.requestName=SingleElectronElMuRunENov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_9.py General.requestName=SingleElectronElMuRunCNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_10.py General.requestName=SingleElectronElMuRunDNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_12.py General.requestName=SingleElectronElMuRunFNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_13.py General.requestName=SingleElectronElMuRunGNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_8.py General.requestName=SingleElectronElMuRunBNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_14.py General.requestName=SingleElectronElMuRunH1Nov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_15.py General.requestName=SingleElectronElMuRunH2Nov15 Data.outputDatasetTag=Nov15ElMuDATA

crab submit -c test/crab_Data_ElMu_16.py General.requestName=SingleMuonElMuRunBNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_17.py General.requestName=SingleMuonElMuRunCNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_18.py General.requestName=SingleMuonElMuRunDNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_19.py General.requestName=SingleMuonElMuRunENov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_20.py General.requestName=SingleMuonElMuRunFNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_21.py General.requestName=SingleMuonElMuRunGNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_22.py General.requestName=SingleMuonElMuRunH1Nov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_23.py General.requestName=SingleMuonElMuRunH2Nov15 Data.outputDatasetTag=Nov15ElMuDATA

#Data MuonEG
crab submit -c test/crab_Data_ElMu_0.py General.requestName=MuonEGRunBNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_1.py General.requestName=MuonEGRunCNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_2.py General.requestName=MuonEGRunDNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_3.py General.requestName=MuonEGRunENov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_4.py General.requestName=MuonEGRunFNov15 Data.outputDatasetTag=Nov15ElMuDATA

crab submit -c test/crab_Data_ElMu_5.py General.requestName=MuonEGRunGNov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_6.py General.requestName=MuonEGRunH1Nov15 Data.outputDatasetTag=Nov15ElMuDATA
crab submit -c test/crab_Data_ElMu_7.py General.requestName=MuonEGRunH2Nov15 Data.outputDatasetTag=Nov15ElMuDATA

#TT
crab submit -c test/crab_MC_ElMu_0.py General.requestName=TTNov15ElMu Data.outputDatasetTag=TTNov15ElMuMC
#ST#tW tbarW
crab submit -c test/crab_MC_ElMu_1.py General.requestName=tbarWNov15ElMu Data.outputDatasetTag=tbarWNov15ElMuMC
crab submit -c test/crab_MC_ElMu_2.py General.requestName=tWNov15ElMu Data.outputDatasetTag=tWNov15ElMuMC

#WW
crab submit -c test/crab_MC_ElMu_3.py General.requestName=WWNov15ElMu Data.outputDatasetTag=WWNov15ElMuMC
#WZ
crab submit -c test/crab_MC_ElMu_4.py General.requestName=WZNov15ElMu Data.outputDatasetTag=WZNov15ElMuMC
#ZZ
crab submit -c test/crab_MC_ElMu_5.py General.requestName=ZZNov15ElMu Data.outputDatasetTag=ZZNov15ElMuMC
#W+Jets
crab submit -c test/crab_MC_ElMu_6.py General.requestName=WJetsNov15ElMu Data.outputDatasetTag=WJetsNov15ElMuMC
#DY
crab submit -c test/crab_MC_ElMu_7.py General.requestName=DY10Nov15ElMu Data.outputDatasetTag=DY10Nov15ElMuMC
crab submit -c test/crab_MC_ElMu_8.py General.requestName=DY50Nov15ElMu Data.outputDatasetTag=DY50Nov15ElMuMC

#TTWJets
crab submit -c test/crab_MC_ElMu_9.py General.requestName=TTWJetsToLNuNov15ElMu Data.outputDatasetTag=TTWJetsToLNuNov15ElMuMC
crab submit -c test/crab_MC_ElMu_10.py General.requestName=TTWJetsToQQNov15ElMu Data.outputDatasetTag=TTWJetsToQQNov15ElMuMC
#TTZ
crab submit -c test/crab_MC_ElMu_12.py General.requestName=TTZToLLNuNuNov15ElMu Data.outputDatasetTag=TTZToLLNuNuNov15ElMuMC
crab submit -c test/crab_MC_ElMu_11.py General.requestName=TTZToQQNov15ElMu Data.outputDatasetTag=TTZToQQNov15ElMuMC





