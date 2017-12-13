# whelicityAnalyzer

cmsrel CMSSW_8_0_26_patch2
cd CMSS_8_0_26_patch2/src
mkdir whelicity1
git clone https://github.com/hessikaveh/whelicityAnalyzer.git
mv whelicityAnalyzer MiniAnalyzer
cmsRun test/whelicity_madgraph_DiEl_cfg.py
