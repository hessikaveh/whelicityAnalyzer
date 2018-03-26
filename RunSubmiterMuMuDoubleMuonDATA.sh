#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs


#Data DoubleMuon
crab submit -c test/crab_Data_DiMu_0.py General.requestName=DoubleMuonRunBNov06 Data.outputDatasetTag=Nov06MuMuDATA
crab submit -c test/crab_Data_DiMu_1.py General.requestName=DoubleMuonRunCNov06 Data.outputDatasetTag=Nov06MuMuDATA
crab submit -c test/crab_Data_DiMu_2.py General.requestName=DoubleMuonRunDNov06 Data.outputDatasetTag=Nov06MuMuDATA
crab submit -c test/crab_Data_DiMu_3.py General.requestName=DoubleMuonRunENov06 Data.outputDatasetTag=Nov06MuMuDATA
crab submit -c test/crab_Data_DiMu_4.py General.requestName=DoubleMuonRunFNov06 Data.outputDatasetTag=Nov06MuMuDATA

crab submit -c test/crab_Data_DiMu_5.py General.requestName=DoubleMuonRunGNov06 Data.outputDatasetTag=Nov06MuMuDATA
crab submit -c test/crab_Data_DiMu_6.py General.requestName=DoubleMuonRunH1Nov06 Data.outputDatasetTag=Nov06MuMuDATA
crab submit -c test/crab_Data_DiMu_7.py General.requestName=DoubleMuonRunH2Nov06 Data.outputDatasetTag=Nov06MuMuDATA


