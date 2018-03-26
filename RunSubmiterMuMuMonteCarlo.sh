#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs


#TT
crab submit -c test/crab_MC_DiMu_0.py General.requestName=TTNov01MuMu Data.outputDatasetTag=TTNov01MuMuMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiMu_1.py General.requestName=tbarWNov01MuMu Data.outputDatasetTag=tbarWNov01MuMuMC
crab submit -c test/crab_MC_DiMu_2.py General.requestName=tWNov01MuMu Data.outputDatasetTag=tWNov01MuMuMC

#WW
crab submit -c test/crab_MC_DiMu_3.py General.requestName=WWNov01MuMu Data.outputDatasetTag=WWNov01MuMuMC
#WZ
crab submit -c test/crab_MC_DiMu_4.py General.requestName=WZNov01MuMu Data.outputDatasetTag=WZNov01MuMuMC
#ZZ
crab submit -c test/crab_MC_DiMu_5.py General.requestName=ZZNov01MuMu Data.outputDatasetTag=ZZNov01MuMuMC
#W+Jets
crab submit -c test/crab_MC_DiMu_6.py General.requestName=WJetsNov01MuMu Data.outputDatasetTag=WJetsNov01MuMuMC
#DY
crab submit -c test/crab_MC_DiMu_7.py General.requestName=DY10Nov01MuMu Data.outputDatasetTag=DY10Nov01MuMuMC
crab submit -c test/crab_MC_DiMu_8.py General.requestName=DY50Nov01MuMu Data.outputDatasetTag=DY50Nov01MuMuMC

#TTWJets
crab submit -c test/crab_MC_DiMu_9.py General.requestName=TTWJetsToLNuNov01MuMu Data.outputDatasetTag=TTWJetsToLNuNov01MuMuMC
crab submit -c test/crab_MC_DiMu_10.py General.requestName=TTWJetsToQQNov01MuMu Data.outputDatasetTag=TTWJetsToQQNov01MuMuMC
#TTZ
crab submit -c test/crab_MC_DiMu_12.py General.requestName=TTZToLLNuNuNov01MuMu Data.outputDatasetTag=TTZToLLNuNuNov01MuMuMC
crab submit -c test/crab_MC_DiMu_11.py General.requestName=TTZToQQNov01MuMu Data.outputDatasetTag=TTZToQQNov01MuMuMC


#TT
crab submit -c test/crab_MC_DiMu_GH_0.py General.requestName=TTNov01MuMuGH Data.outputDatasetTag=TTNov01MuMuGHMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiMu_GH_1.py General.requestName=tbarWNov01MuMuGH Data.outputDatasetTag=tbarWNov01MuMuGHMC
crab submit -c test/crab_MC_DiMu_GH_2.py General.requestName=tWNov01MuMuGH Data.outputDatasetTag=tWNov01MuMuGHMC

#WW
crab submit -c test/crab_MC_DiMu_GH_3.py General.requestName=WWNov01MuMuGH Data.outputDatasetTag=WWNov01MuMuGHMC
#WZ
crab submit -c test/crab_MC_DiMu_GH_4.py General.requestName=WZNov01MuMuGH Data.outputDatasetTag=WZNov01MuMuGHMC
#ZZ
crab submit -c test/crab_MC_DiMu_GH_5.py General.requestName=ZZNov01MuMuGH Data.outputDatasetTag=ZZNov01MuMuGHMC
#W+Jets
crab submit -c test/crab_MC_DiMu_GH_6.py General.requestName=WJetsNov01MuMuGH Data.outputDatasetTag=WJetsNov01MuMuGHMC
#DY
crab submit -c test/crab_MC_DiMu_GH_7.py General.requestName=DY10Nov01MuMuGH Data.outputDatasetTag=DY10Nov01MuMuGHMC
crab submit -c test/crab_MC_DiMu_GH_8.py General.requestName=DY50Nov01MuMuGH Data.outputDatasetTag=DY50Nov01MuMuGHMC

#TTWJets
crab submit -c test/crab_MC_DiMu_GH_9.py General.requestName=TTWJetsToLNuNov01MuMuGH Data.outputDatasetTag=TTWJetsToLNuNov01MuMuGHMC
crab submit -c test/crab_MC_DiMu_GH_10.py General.requestName=TTWJetsToQQNov01MuMuGH Data.outputDatasetTag=TTWJetsToQQNov01MuMuGHMC
#TTZ
crab submit -c test/crab_MC_DiMu_GH_12.py General.requestName=TTZToLLNuNuNov01MuMuGH Data.outputDatasetTag=TTZToLLNuNuNov01MuMuGHMC
crab submit -c test/crab_MC_DiMu_GH_11.py General.requestName=TTZToQQNov01MuMuGH Data.outputDatasetTag=TTZToQQNov01MuMuGHMC


