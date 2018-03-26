#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs

#TT
crab submit -c test/crab_MC_DiEl_0.py General.requestName=TTOct30ElEl Data.outputDatasetTag=TTOct30ElElMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiEl_1.py General.requestName=tbarWOct30ElEl Data.outputDatasetTag=tbarWOct30ElElMC
crab submit -c test/crab_MC_DiEl_2.py General.requestName=tWOct30ElEl Data.outputDatasetTag=tWOct30ElElMC

#WW
crab submit -c test/crab_MC_DiEl_3.py General.requestName=WWOct30ElEl Data.outputDatasetTag=WWOct30ElElMC
#WZ
crab submit -c test/crab_MC_DiEl_4.py General.requestName=WZOct30ElEl Data.outputDatasetTag=WZOct30ElElMC
#ZZ
crab submit -c test/crab_MC_DiEl_5.py General.requestName=ZZOct30ElEl Data.outputDatasetTag=ZZOct30ElElMC
#W+Jets
crab submit -c test/crab_MC_DiEl_6.py General.requestName=WJetsOct30ElEl Data.outputDatasetTag=WJetsOct30ElElMC
#DY
crab submit -c test/crab_MC_DiEl_7.py General.requestName=DY10Oct30ElEl Data.outputDatasetTag=DY10Oct30ElElMC
crab submit -c test/crab_MC_DiEl_8.py General.requestName=DY50Oct30ElEl Data.outputDatasetTag=DY50Oct30ElElMC

#TTWJets
crab submit -c test/crab_MC_DiEl_9.py General.requestName=TTWJetsToLNuOct30ElEl Data.outputDatasetTag=TTWJetsToLNuOct30ElElMC
crab submit -c test/crab_MC_DiEl_10.py General.requestName=TTWJetsToQQOct30ElEl Data.outputDatasetTag=TTWJetsToQQOct30ElElMC
#TTZ
crab submit -c test/crab_MC_DiEl_12.py General.requestName=TTZToLLNuNuOct30ElEl Data.outputDatasetTag=TTZToLLNuNuOct30ElElMC
crab submit -c test/crab_MC_DiEl_11.py General.requestName=TTZToQQOct30ElEl Data.outputDatasetTag=TTZToQQOct30ElElMC


#TT
crab submit -c test/crab_MC_DiEl_GH_0.py General.requestName=TTOct30ElElGH Data.outputDatasetTag=TTOct30ElElGHMC
#ST#tW tbarW
crab submit -c test/crab_MC_DiEl_GH_1.py General.requestName=tbarWOct30ElElGH Data.outputDatasetTag=tbarWOct30ElElGHMC
crab submit -c test/crab_MC_DiEl_GH_2.py General.requestName=tWOct30ElElGH Data.outputDatasetTag=tWOct30ElElGHMC

#WW
crab submit -c test/crab_MC_DiEl_GH_3.py General.requestName=WWOct30ElElGH Data.outputDatasetTag=WWOct30ElElGHMC
#WZ
crab submit -c test/crab_MC_DiEl_GH_4.py General.requestName=WZOct30ElElGH Data.outputDatasetTag=WZOct30ElElGHMC
#ZZ
crab submit -c test/crab_MC_DiEl_GH_5.py General.requestName=ZZOct30ElElGH Data.outputDatasetTag=ZZOct30ElElGHMC
#W+Jets
crab submit -c test/crab_MC_DiEl_GH_6.py General.requestName=WJetsOct30ElElGH Data.outputDatasetTag=WJetsOct30ElElGHMC
#DY
crab submit -c test/crab_MC_DiEl_GH_7.py General.requestName=DY10Oct30ElElGH Data.outputDatasetTag=DY10Oct30ElElGHMC
crab submit -c test/crab_MC_DiEl_GH_8.py General.requestName=DY50Oct30ElElGH Data.outputDatasetTag=DY50Oct30ElElGHMC

#TTWJets
crab submit -c test/crab_MC_DiEl_GH_9.py General.requestName=TTWJetsToLNuOct30ElElGH Data.outputDatasetTag=TTWJetsToLNuOct30ElElGHMC
crab submit -c test/crab_MC_DiEl_GH_10.py General.requestName=TTWJetsToQQOct30ElElGH Data.outputDatasetTag=TTWJetsToQQOct30ElElGHMC
#TTZ
crab submit -c test/crab_MC_DiEl_GH_12.py General.requestName=TTZToLLNuNuOct30ElElGH Data.outputDatasetTag=TTZToLLNuNuOct30ElElGHMC
crab submit -c test/crab_MC_DiEl_GH_11.py General.requestName=TTZToQQOct30ElElGH Data.outputDatasetTag=TTZToQQOct30ElElGHMC


