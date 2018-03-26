#!/bin/bash
echo This is the WHelicity analysis
echo First Lets Compile

scram b -j4

echo Now Lets Add Runs
#Data DoubleEG Single
crab submit -c test/crab_Data_DiEl_11.py General.requestName=SingleElectronElElRunEOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_9.py General.requestName=SingleElectronElElRunCOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_10.py General.requestName=SingleElectronElElRunDOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_12.py General.requestName=SingleElectronElElRunFOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_13.py General.requestName=SingleElectronElElRunGOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_8.py General.requestName=SingleElectronElElRunBOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_14.py General.requestName=SingleElectronElElRunH1Oct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_15.py General.requestName=SingleElectronElElRunH2Oct30 Data.outputDatasetTag=Oct30ElElDATA
#Data DoubleEG
crab submit -c test/crab_Data_DiEl_0.py General.requestName=DoubleEGRunBOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_1.py General.requestName=DoubleEGRunCOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_2.py General.requestName=DoubleEGRunDOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_3.py General.requestName=DoubleEGRunEOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_4.py General.requestName=DoubleEGRunFOct30 Data.outputDatasetTag=Oct30ElElDATA

crab submit -c test/crab_Data_DiEl_5.py General.requestName=DoubleEGRunGOct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_6.py General.requestName=DoubleEGRunH1Oct30 Data.outputDatasetTag=Oct30ElElDATA
crab submit -c test/crab_Data_DiEl_7.py General.requestName=DoubleEGRunH2Oct30 Data.outputDatasetTag=Oct30ElElDATA

