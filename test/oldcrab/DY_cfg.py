import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
## configure process options
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/125D0D47-6F30-E611-8C9B-848F69FD29D0.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/12FDA0DB-E730-E611-90BE-001E67E5E8B6.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/16FFE3C3-7F30-E611-9CC4-001D09FDD7C8.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/24BC8143-6F30-E611-8465-848F69FD2958.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/3655A415-E830-E611-A245-00259074AE2E.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/40C90AE1-E730-E611-8109-848F69FD4E98.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/5EC9BCD9-E730-E611-9226-44A842CFC98B.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/6A5956EF-E730-E611-9A44-001E67396BB7.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/72DAA04A-8C30-E611-BF00-44A84225CFF0.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/7C5E8756-E830-E611-8803-2C600CAFEF72.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/823BF0D3-9730-E611-92C6-901B0E542962.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/84B7344D-7F30-E611-A463-848F69FD4667.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/84B73DC2-7E30-E611-9230-002590A831B4.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/9017F570-5C30-E611-999E-0025905D1D7A.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/96100D0C-9B30-E611-87CB-00259090846E.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/A450FF50-7F30-E611-97DE-001D09FDD7C8.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/A4E946FD-E730-E611-9653-0CC47A4DED94.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/A8AC3402-E830-E611-99AF-D4AE526A0488.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/D456C801-E830-E611-9112-549F358EB721.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/D89D9691-7B30-E611-B10D-0CC47A4D76AC.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/E6671938-E830-E611-888A-0025905C95F8.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/F40E4002-E830-E611-A9E8-0025907BAD4A.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/FEC2EEE3-7630-E611-B0CB-002590D0B00E.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/025FD935-9332-E611-98CE-FA163E472E2A.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/06507482-E530-E611-B82F-7845C4FC37AF.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/0EE3D05D-E530-E611-B1B2-7845C4FC3A4C.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/10CE32EA-F030-E611-B73B-008CFA002FF4.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/16910939-6F32-E611-88B1-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/1803F72B-7032-E611-ACDA-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/24E4DCF4-6F32-E611-8E62-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/26D10939-6F32-E611-9528-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/2ACBC294-9332-E611-AEBB-02163E014982.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/2CDC2754-9332-E611-A091-02163E01325B.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/2CFCC83C-9432-E611-9518-02163E014919.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/2E08FEE8-E631-E611-BE45-848F69FD28AA.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/322F4BE9-E631-E611-9856-180373FF8CD4.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/44732346-6F32-E611-9CF0-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/44D72346-6F32-E611-BF74-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/4CCE0939-6F32-E611-87ED-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/4E91B898-7832-E611-89CC-F04DA275C2CE.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/50475D9B-7832-E611-87CB-00266CF9B878.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/5AB05499-9232-E611-9D69-02163E00F7B3.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/5E8A2446-6F32-E611-B186-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/64591533-EB30-E611-BCAC-848F69FD29AF.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/669D32F7-E630-E611-99B4-7845C4FC3A4C.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/6C900939-6F32-E611-A039-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/6E590633-EB30-E611-912C-848F69FD29AF.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/74AC2346-6F32-E611-9CB3-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/7C1CFC7B-7032-E611-961C-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/7CBE8D9D-EA30-E611-9670-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/8A900939-6F32-E611-A4F3-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/8E0C3D1E-BE31-E611-9D17-7845C4FC3A40.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/90082446-6F32-E611-822C-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/9483C67E-9232-E611-B134-02163E014BF0.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/96CC2FD6-EC30-E611-8F39-848F69FD45CB.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/98960939-6F32-E611-941C-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/A4EE8B59-F130-E611-B995-180373FF8D6A.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/A88B2446-6F32-E611-A37E-7845C4F92F87.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/B48482FE-E630-E611-BB35-7845C4FC37AF.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/B616C81A-BE31-E611-974F-00266CF9B878.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/BE7D87EB-C331-E611-A264-7845C4F9321B.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/E2B00939-6F32-E611-BB6C-7845C4F91450.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/EA701333-EB30-E611-9E58-848F69FD29AF.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/EC259B2B-EB30-E611-AF24-F04DA275C2CE.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/EEA871E7-C331-E611-9426-848F69FD2961.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/F666FB12-9332-E611-AC32-02163E0177E5.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/FCA94C53-9032-E611-B26A-848F69FD292B.root',
'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/FCFB2346-6F32-E611-BFF0-7845C4F92F87.root'
    )
)
genParticleCollection = 'prunedGenParticles'
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
process.initSubset.src = genParticleCollection
process.decaySubset.src = genParticleCollection
process.decaySubset.runMode = "Run2"
process.load("whelicity1.MiniAnalyzer.whelicity_cff")

process.TFileService = cms.Service("TFileService",
 fileName = cms.string("DY2Lhistos.root")
)


process.p = cms.Path(process.Whelicity)
