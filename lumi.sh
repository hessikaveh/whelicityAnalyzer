#!/bin/bash
echo This is the WHelicity analysis lumi calculation

 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunBJan13/results/processedLumis.json  > MuMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunCJan13/results/processedLumis.json  >> MuMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunDJan13/results/processedLumis.json  >> MuMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunEJan13/results/processedLumis.json  >> MuMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunFJan13/results/processedLumis.json  >> MuMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunGJan13/results/processedLumis.json  >> MuMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunH1Jan13/results/processedLumis.json  >> MuMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuonRunH2Jan13/results/processedLumis.json  >> MuMuLumi.txt
# ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleMuRunDApr3/results/processedLumis.json --hltpath HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v* >> MuMuLumi.txt

echo Now Lets calculate the DoubleEG part
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunBJan14/results/processedLumis.json    > ElElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunCJan14/results/processedLumis.json  >> ElElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunDJan14/results/processedLumis.json  >> ElElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunEJan14/results/processedLumis.json  >> ElElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunFJan14/results/processedLumis.json  >> ElElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunGJan14/results/processedLumis.json   >> ElElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunH1Jan14/results/processedLumis.json   >> ElElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_DoubleEGRunH2Jan14/results/processedLumis.json   >> ElElLumi.txt

 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunBJan07/results/processedLumis.json   > SingleElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunCJan07/results/processedLumis.json   >> SingleElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunDJan07/results/processedLumis.json   >> SingleElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunEJan07/results/processedLumis.json   >> SingleElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunFJan07/results/processedLumis.json   >> SingleElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunGJan07/results/processedLumis.json   >> SingleElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunH1Jan07/results/processedLumis.json   >> SingleElLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_SingleElectronElElRunH2Jan07/results/processedLumis.json   >> SingleElLumi.txt
echo Now Lets calculate the MuonEG part
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunBJan07/results/processedLumis.json  > ElMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunCJan07/results/processedLumis.json  >> ElMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunDJan07/results/processedLumis.json  >> ElMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunEJan07/results/processedLumis.json  >> ElMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunFJan07/results/processedLumis.json  >> ElMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunGJan07/results/processedLumis.json  >> ElMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunH1Jan07/results/processedLumis.json   >> ElMuLumi.txt
 ~/brilconda/bin/brilcalc lumi -c web -b "STABLE BEAMS" --normtag normtag_DATACERT.json -u /pb -i crab_projects/crab_MuonEGRunH2Jan07/results/processedLumis.json  >> ElMuLumi.txt

echo Done!
