#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs

#--------------------------------------
#Data DoubleMuon Single
crab submit -c test/crab_Data_DiMu_11.py General.requestName=SingleMuonMuMuRunENov01 Data.outputDatasetTag=Nov01MuMuDATA
crab submit -c test/crab_Data_DiMu_9.py General.requestName=SingleMuonMuMuRunCNov01 Data.outputDatasetTag=Nov01MuMuDATA
crab submit -c test/crab_Data_DiMu_10.py General.requestName=SingleMuonMuMuRunDNov01 Data.outputDatasetTag=Nov01MuMuDATA
crab submit -c test/crab_Data_DiMu_12.py General.requestName=SingleMuonMuMuRunFNov01 Data.outputDatasetTag=Nov01MuMuDATA
crab submit -c test/crab_Data_DiMu_13.py General.requestName=SingleMuonMuMuRunGNov01 Data.outputDatasetTag=Nov01MuMuDATA
crab submit -c test/crab_Data_DiMu_8.py General.requestName=SingleMuonMuMuRunBNov01 Data.outputDatasetTag=Nov01MuMuDATA
crab submit -c test/crab_Data_DiMu_14.py General.requestName=SingleMuonMuMuRunH1Nov01 Data.outputDatasetTag=Nov01MuMuDATA
crab submit -c test/crab_Data_DiMu_15.py General.requestName=SingleMuonMuMuRunH2Nov01 Data.outputDatasetTag=Nov01MuMuDATA

