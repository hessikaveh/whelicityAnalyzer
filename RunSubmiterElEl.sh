#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs
#Data DoubleEG Single
crab submit -c test/crab_Data_DiEl_11.py General.requestName=SingleElectronElElRunEJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_9.py General.requestName=SingleElectronElElRunCJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_10.py General.requestName=SingleElectronElElRunDJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_12.py General.requestName=SingleElectronElElRunFJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_13.py General.requestName=SingleElectronElElRunGJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_8.py General.requestName=SingleElectronElElRunBJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_14.py General.requestName=SingleElectronElElRunH1Jan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_15.py General.requestName=SingleElectronElElRunH2Jan14 Data.outputDatasetTag=Jan14ElElDATA
#Data DoubleEG
crab submit -c test/crab_Data_DiEl_0.py General.requestName=DoubleEGRunBJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_1.py General.requestName=DoubleEGRunCJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_2.py General.requestName=DoubleEGRunDJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_3.py General.requestName=DoubleEGRunEJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_4.py General.requestName=DoubleEGRunFJan14 Data.outputDatasetTag=Jan14ElElDATA

crab submit -c test/crab_Data_DiEl_5.py General.requestName=DoubleEGRunGJan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_6.py General.requestName=DoubleEGRunH1Jan14 Data.outputDatasetTag=Jan14ElElDATA
crab submit -c test/crab_Data_DiEl_7.py General.requestName=DoubleEGRunH2Jan14 Data.outputDatasetTag=Jan14ElElDATA

#TT
crab submit -c test/crab_MC_DiEl_0.py General.requestName=TTJan14ElEl Data.outputDatasetTag=TTJan14ElElMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiEl_1.py General.requestName=tbarWJan14ElEl Data.outputDatasetTag=tbarWJan14ElElMC
crab submit -c test/crab_MC_DiEl_2.py General.requestName=tWJan14ElEl Data.outputDatasetTag=tWJan14ElElMC

#WW
crab submit -c test/crab_MC_DiEl_3.py General.requestName=WWJan14ElEl Data.outputDatasetTag=WWJan14ElElMC
#WZ
crab submit -c test/crab_MC_DiEl_4.py General.requestName=WZJan14ElEl Data.outputDatasetTag=WZJan14ElElMC
#ZZ
crab submit -c test/crab_MC_DiEl_5.py General.requestName=ZZJan14ElEl Data.outputDatasetTag=ZZJan14ElElMC
#W+Jets
crab submit -c test/crab_MC_DiEl_6.py General.requestName=WJetsJan14ElEl Data.outputDatasetTag=WJetsJan14ElElMC
#DY
crab submit -c test/crab_MC_DiEl_7.py General.requestName=DY10Jan14ElEl Data.outputDatasetTag=DY10Jan14ElElMC
crab submit -c test/crab_MC_DiEl_8.py General.requestName=DY50Jan14ElEl Data.outputDatasetTag=DY50Jan14ElElMC

#TTWJets
crab submit -c test/crab_MC_DiEl_9.py General.requestName=TTWJetsToLNuJan14ElEl Data.outputDatasetTag=TTWJetsToLNuJan14ElElMC
crab submit -c test/crab_MC_DiEl_10.py General.requestName=TTWJetsToQQJan14ElEl Data.outputDatasetTag=TTWJetsToQQJan14ElElMC
#TTZ
crab submit -c test/crab_MC_DiEl_12.py General.requestName=TTZToLLNuNuJan14ElEl Data.outputDatasetTag=TTZToLLNuNuJan14ElElMC
crab submit -c test/crab_MC_DiEl_11.py General.requestName=TTZToQQJan14ElEl Data.outputDatasetTag=TTZToQQJan14ElElMC


