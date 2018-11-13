###############################
####### Parameters ############
#1;95;0c##############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('python')

options.register('outFilename',
    'Maps/test_withMuonsVeto.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'output file name'
)
options.register('reportEvery',
    1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Report every N events (default is N=1000)'
)
options.register('wantSummary',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 150000)

options.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process('USER')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        "/store/mc/RunIISummer16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00DB3FA3-02B6-E611-8A56-9CB65404FC30.root",
        "/store/mc/RunIISummer16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/02FBE5E0-1AB6-E611-85B1-0CC47A7DFF82.root",
        "/store/mc/RunIISummer16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/04CD4E12-0CB6-E611-8AB6-E41D2D08DDF0.root",
        "/store/mc/RunIISummer16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/06D23025-06B6-E611-804A-001E6739AC59.root",
        "/store/mc/RunIISummer16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0EBEA982-09B6-E611-A978-A0000220FE80.root",
        #"/store/mc/RunIISummer16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00DB3FA3-02B6-E611-8A56-9CB65404FC30.root",
        #"/store/mc/RunIISummer16MiniAODv2/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0809E9C0-23BB-E611-A10E-6C3BE5B5B078.root",
        #"/store/mc/RunIISummer16MiniAODv2/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00559533-06B7-E611-8594-0CC47AD9914C.root",
        #"/store/mc/RunIISummer16MiniAODv2/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00D17FD4-8EBD-E611-B17D-002590D0AFC2.root",
        #"/store/mc/RunIISummer16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/000316AF-9FBE-E611-9761-0CC47A7C35F8.root",
        #"/store/mc/RunIISummer16MiniAODv2/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/029929BB-ACBD-E611-96C1-002590DE6C9A.root",
        #"/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root",

        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/2670E1F5-8E39-E611-ACB3-002590D0B05E.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/2E241F23-C63A-E611-A3E6-0CC47AB35EF0.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/6A41327D-F639-E611-9249-00266CFFCD14.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/8A621648-653A-E611-A85B-001E67F8F9FC.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/20802904-1F3A-E611-AE7F-0CC47A6C17FC.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/C23993D9-A239-E611-8AA8-B083FED07198.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v2/20000/12FD7419-4F3F-E611-8F04-002590AC4B5C.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/180C6D56-143A-E611-9A3B-F04DA275407C.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/081AFD03-BB39-E611-B4DD-0025904B8934.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/1A51CAF8-D739-E611-9317-0025902BD8CE.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/94BF057F-873A-E611-97A3-0CC47AB35C66.root",
        #"/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/0C55FF4B-B939-E611-A3F9-001E67F3332A.root",

        #"/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/004A0552-3929-E611-BD44-0025905A48F0.root",
        #"/store/backfill/1/RunIISpring16MiniAODv2/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/610000/00424948-1827-E611-8598-0CC47A1DF1A6.root",
        #"/store/mc/RunIISpring16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/001FED44-BA3E-E611-8DAF-0025905C5474.root",
        #"/store/mc/RunIISpring16MiniAODv2/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/00223778-771A-E611-B535-3417EBE5280A.root",
        #"/store/mc/RunIISpring16MiniAODv2/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/08540B3F-5F1C-E611-ADC7-00266CF3E128.root",
        #"/store/mc/RunIISpring16MiniAODv2/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/163EF19E-021D-E611-A0D4-002590E7DFFE.root",
        #"/store/mc/RunIISpring16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/000843D6-AC1C-E611-AB18-0025901A9EFC.root"

        #Bprime
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082629/0000/B2GEDMNtuple_1.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082705/0000/B2GEDMNtuple_2.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082646/0000/B2GEDMNtuple_2.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082629/0000/B2GEDMNtuple_1.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082612/0000/B2GEDMNtuple_2.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082554/0000/B2GEDMNtuple_1.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082537/0000/B2GEDMNtuple_1.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/B2GAnaFW_80X_V2p1/BprimeBToHB_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_80X_V2p1/161014_064618/0000/B2GEDMNtuple_1.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_082504/0000/B2GEDMNtuple_2.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_111936/0000/B2GEDMNtuple_2.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_111914/0000/B2GEDMNtuple_2.root',
        #'/store/user/grauco/B2GAnaFwk80x_v2p1/B2GAnaFW/v80x_v2p1/BprimeBToHB_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8/v80x_v2p1/161013_111851/0000/B2GEDMNtuple_1.root',
        #QCD
        #'/store/group/phys_b2g/B2GAnaFW_80X_V2p1/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1_B2GAnaFW_80X_V2p1/161018_194944/0000/B2GEDMNtuple_10.root"
        #'/store/group/phys_b2g/B2GAnaFW_80X_V2p1/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1/161018_195000/0000/B2GEDMNtuple_10.root',
        #'/store/user/grauco/B2GAnaFW/B2GAnaFW_80X_V2p1/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p1/161018_070023/0000/B2GEDMNtuple_1.root',
        #'/store/user/grauco/B2GAnaFW/B2GAnaFW_80X_V2p1/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p1/161018_070036/0000/B2GEDMNtuple_101.root',
        #'/store/user/grauco/B2GAnaFW/B2GAnaFW_80X_V2p1/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p1/161018_070104/0000/B2GEDMNtuple_1.root',
        #'/store/mc/RunIISpring16MiniAODv2/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/00000/02525F30-F827-E611-AF32-0CC47A74527A.root',


        #'/store/mc/RunIISpring16MiniAODv2/BprimeBToHB_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/2670E1F5-8E39-E611-ACB3-002590D0B05E.root',
        #'file:B2GEDMNtuple.root',
    )
)

process.TFileService = cms.Service('TFileService',
   fileName = cms.string(options.outFilename)
)

### Selected leptons and jets
process.skimmedPatMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 10.0 && abs(eta) < 2.4 && isLooseMuon")
    )

process.bTaggingEffAnalyzerAK4PF = cms.EDAnalyzer('bTaggingEffAnalyzer',
    JetsTag            = cms.InputTag('slimmedJets'),
    MuonsTag           = cms.InputTag("slimmedMuons"),
    DiscriminatorTag   = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    DiscriminatorValue = cms.double(0.8484),
    PtNBins            = cms.int32(100),
    PtMin              = cms.double(0.),
    PtMax              = cms.double(1000.),
    EtaNBins           = cms.int32(60),
    EtaMin             = cms.double(-2.4),
    EtaMax             = cms.double(2.4)
)

#process.bTaggingEffAnalyzerAK4PF = cms.EDAnalyzer('bTaggingEffAnalyzer',
#    JetsTag            = cms.InputTag('slimmedJetsAK8PFCHSSoftDropPacked:SubJets'),
#    DiscriminatorTag   = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
#    DiscriminatorValue = cms.double(0.8484),
#    PtNBins            = cms.int32(100),
#    PtMin              = cms.double(0.),
#    PtMax              = cms.double(1000.),
#    EtaNBins           = cms.int32(60),
#    EtaMin             = cms.double(-2.4),
#    EtaMax             = cms.double(2.4)
#)

process.p = cms.Path(process.bTaggingEffAnalyzerAK4PF)
