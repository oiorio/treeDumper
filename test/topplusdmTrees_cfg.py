import commands, os
### *****************************************************************************************
### Usage:
###
### cmsRun topplusdmanaEDMntuples_cfg.py maxEvts=N sample="mySample/sample.root" version="1"7 outputLabel="myoutput"
###
### Default values for the options are set:
### maxEvts     = -1
### sample      = 'file:/scratch/decosa/ttDM/testSample/tlbsm_53x_v3_mc_10_1_qPV.root'
### outputLabel = 'analysisTTDM.root'
### *****************************************************************************************
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')

options.register('maxEvts',
                 -1,# default value: process all events
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.int,
                 'Number of events to process')

SingleElTriggers = []
SingleMuTriggers = []
hadronTriggers = []
HadronPFHT800Triggers = []
HadronPFHT900Triggers = []
HadronPFJet450Triggers = []

chan = "MET_Prompt"

chan = "TTbarDMJets_scalar_Mchi-50_Mphi-50"

chan = "DY"
chan = "WJ"
filedir= "/tmp/oiorio/"
cmd = "ls "+filedir+"/"+chan+"/"

status,ls_la = commands.getstatusoutput( cmd )
listFiles = ls_la.split(os.linesep)
files = []
#files = ["file:re-MiniAOD_17Jul/"+l for l in listFiles]
#files = ["file:"+filedir+"MET_Prompt/"+l for l in listFiles]
files = ["file:"+filedir+"/"+chan+"/"+l for l in listFiles]
options.register('sample',
                 #'/store/user/decosa/B2GAnaFW_80X_v2p4/BprimeBToHB_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/BprimeBToHB_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_v2p4/170213_085503/0000/B2GEDMNtuple_2.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/170104_182124/0000/B2GEDMNtuple_104.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016C-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_155142/0000/B2GEDMNtuple_10.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016H-PromptReco-v2_B2GAnaFW_80X_v2p4/161221_180745/0000/B2GEDMNtuple_10.root',
                 #'/store/user/grauco/B2GAnaFW_80X_V2p4/BprimeBToHB_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/BprimeBToHB_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161227_111236/0000/B2GEDMNtuple_2.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p3/JetHT/Run2016B/JetHT/Run2016B-23Sep2016-v3_B2GAnaFW_80X_V2p3/161216_214635/0000/B2GEDMNtuple_1.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1_B2GAnaFW_80X_V2p4/170124_202648/0000/B2GEDMNtuple_1.root',

#Recent files:
               #  '/store/user/decosa/EPS17/B2GAnaFW_80X_V3p1/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_103545/0000/B2GEDMNtuple_1.root',
                 'file:B2GEDMNtuple_1.root',#l'ho preso dal mio TT
                 #'file:/eos/user/o/oiorio/samples/synch/mc_jecv4/B2GSynchMC.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')

options.register('version',
                 #'53',
                 '71',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'ntuple version (53 or 71)')

options.register('outputLabel',
                 'treesTTJets.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

options.register('EraLabel',
                 'BCD',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Data Era Label')

options.register('isData',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Is data?')

options.register('applyRes',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'ApplyResiduals?')

options.register('addPartonInfo',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Add parton info??')

options.register('changeJECs',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Apply new JECs?')

options.register('LHE',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Keep LHEProducts')

options.register('lhesource',
                 'externalLHEProducer',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'LHEProducts source')

options.register('channel',
                 "",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'channel for weight evaluation'
                 )

options.parseArguments()

if(options.isData):options.LHE = False

if(not options.isData): options.applyRes = False

l = ["singleTrigger"+str(s) for s in xrange(15)]
l = l + ["trigger2"+str(s) for s in xrange(15)]

SingleElTriggers = ["HLT_Ele27_eta2p1_WP75_Gsf"]
SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WP75_Gsf_v"+str(s) for s in xrange(15)]

PhotonTriggers = ["HLT_AK8PFJet360_TrimMass30"]
PhotonTriggers = PhotonTriggers + ["HLT_AK8PFJet360_TrimMass30_v"+str(s) for s in xrange(15)]

SingleMuTriggers = ["HLT_Mu50","HLT_TkMu50"]
SingleMuTriggers = SingleMuTriggers + ["HLT_Mu50_v"+str(s) for s in xrange(15)]
SingleMuTriggers = SingleMuTriggers + ["HLT_TkMu50_v"+str(s) for s in xrange(15)]

HadronPFHT900Triggers = ["HLT_PFHT900"]
HadronPFHT900Triggers = HadronPFHT900Triggers + ["HLT_PFHT900_v"+str(s) for s in xrange(15)]

HadronPFHT800Triggers = ["HLT_PFHT800"]
HadronPFHT800Triggers =HadronPFHT800Triggers + ["HLT_PFHT800_v"+str(s) for s in xrange(15)]

HadronPFJet450Triggers = ["HLT_PFJet450"]
HadronPFJet450Triggers =HadronPFJet450Triggers + ["HLT_PFJet450_v"+str(s) for s in xrange(15)]
#HadronPFJet450Triggers =HadronPFJet450Triggers + ["HLT_PFJet450_v"+str(s) for s in xrange(15)]

#hadronTriggers = ["HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v7","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v7", "HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v3", "HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v2", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v2", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v2", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v2", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v2",  "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight","HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v7", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v7", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v7"]
hadronTriggers = ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"]

#hadronTriggers = hadronTriggers + ["HLT_PFMET90_PFMHT190_IDTight_v"+str(s) for s in xrange(15)]
#hadronTriggers = hadronTriggers + ["HLT_PFMET100_PFMHT100_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]       
#hadronTriggers = hadronTriggers + ["HLT_PFMET120_NoiseCleaned_PFMHT120_IDTight_v"+str(s) for s in xrange(15)]
#hadronTriggers = hadronTriggers + ["HLT_PFMET120_JetIdCleaned_PFMHT120_IDTight_v"+str(s) for s in xrange(15)]
#hadronTriggers = hadronTriggers + ["HLT_PFMET90_NoiseCleaned_PFMHT90_IDTight_v"+str(s) for s in xrange(15)]
#hadronTriggers = hadronTriggers + ["HLT_PFMET90_JetIdCleaned_PFMHT90_IDTight_v"+str(s) for s in xrange(15)]

if(options.isData):

    SingleElTriggers = ["HLT_Ele27_eta2p1_WPLoose_Gsf", "HLT_Ele23_WPLoose_Gsf"]
    SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WPLoose_Gsf_v"+str(s) for s in xrange(15)]
    SingleElTriggers = SingleElTriggers + ["HLT_Ele23_WPLoose_Gsf_v"+str(s) for s in xrange(15)]

    SingleMuTriggers = ["HLT_Mu50","HLT_TkMu50"]
    SingleMuTriggers = SingleMuTriggers + ["HLT_Mu50_v"+str(s) for s in xrange(15)]
    SingleMuTriggers = SingleMuTriggers + ["HLT_TkMu50_v"+str(s) for s in xrange(15)]

    HadronPFHT900Triggers = ["HLT_PFHT900"]
    HadronPFHT900Triggers = HadronPFHT900Triggers + ["HLT_PFHT900_v"+str(s) for s in xrange(15)]
    
    HadronPFHT800Triggers = ["HLT_PFHT800"]
    HadronPFHT800Triggers =HadronPFHT800Triggers + ["HLT_PFHT800_v"+str(s) for s in xrange(15)]
    
    HadronPFJet450Triggers = ["HLT_PFJet450"]
    HadronPFJet450Triggers =HadronPFJet450Triggers + ["HLT_PFJet450_v"+str(s) for s in xrange(15)]
    #HadronPFJet450Triggers =HadronPFJet450Triggers + ["HLT_PFJet450_v"+str(s) for s in xrange(15)]

    #hadronTriggers = ["HLT_PFMET90_PFMHT90_IDTight_v1", "HLT_PFMET100_PFMHT100_IDTight_v1", "HLT_PFMET100_PFMHT100_IDTight_v1", "HLT_PFMET120_PFMHT120_IDTight_v1", "HLT_PFMET120_NoiseCleaned_PFMHT120_IDTight","HLT_PFMET120_JetIdCleaned_PFMHT120_IDTight","HLT_PFMET90_JetIdCleaned_PFMHT90_IDTight","HLT_PFMET90_NoiseCleaned_PFMHT90_IDTight","HLT_PFMET90_NoiseCleaned_PFMHT90_NoID", "HLT_PFMET90_PFMHT90_IDTight", "HLT_PFMET90_PFMHT90_IDTight", "HLT_PFMET100_PFMHT100_IDTight_v7", "HLT_PFMET110_PFMHT110_IDTight_v7", "HLT_PFMET120_PFMHT120_IDTight_v7", "HLT_PFMET90_PFMHT90_IDTight_v2", "HLT_PFMET100_PFMHT100_IDTight_v2", "HLT_PFMET110_PFMHT110_IDTight_v2", "HLT_PFMET120_PFMHT120_IDTight_v2", "HLT_AK8PFJet200_v1", "HLT_AK8PFJet260_v1"]

    #hadronTriggers = hadronTriggers+ ["HLT_PFMET120_NoiseCleaned_PFMHT120_IDTight_v"+str(s) for s in xrange(15)]
    #hadronTriggers = hadronTriggers+ ["HLT_PFMET120_JetIdCleaned_PFMHT120_IDTight_v"+str(s) for s in xrange(15)]
    #hadronTriggers = hadronTriggers+ ["HLT_PFMET90_JetIdCleaned_PFMHT90_IDTight_v"+str(s) for s in xrange(15)]
    #hadronTriggers = hadronTriggers+ ["HLT_PFMET90_NoiseCleaned_PFMHT90_IDTight_v"+str(s) for s in xrange(15)]
    #hadronTriggers = hadronTriggers+ ["HLT_PFMET90_NoiseCleaned_PFMHT90_NoID_v"+str(s) for s in xrange(15)]
    #hadronTriggers = hadronTriggers+ ["HLT_PFMET90_PFMHT90_IDTight_v"+str(s) for s in xrange(15)]

    hadronTriggers = ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"]

    #hadronTriggers = hadronTriggers + ["HLT_PFMET90_PFMHT190_IDTight_v"+str(s) for s in xrange(15)]
    #hadronTriggers = hadronTriggers + ["HLT_PFMET100_PFMHT100_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]

process = cms.Process("ttDManalysisTrees")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('HLTrigReport')
### Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
### Number of maximum events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts) )
### Source file
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#        "root://cms-xrd-global.cern.ch//store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_1.root",
#        "/store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_1.root",
#        "/store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_2.root",
#        "/store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_3.root",
#        "/store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_4.root",
#        "/store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_5.root",
#        "/store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_6.root",
#        "/store/user/ggiannin/B2GAnaFW_80X_V3p1/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/TprimeBToTZ_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_125043/0000/B2GEDMNtuple_7.root",
#        "file:/tmp/oiorio/samples/B2GEDMNtuple_1.root",
##        "file:/tmp/oiorio/samples/B2GEDMNtuple_2.root",
#       "file:/tmp/oiorio/samples/B2GEDMNtuple_3.root",
#        "file:/tmp/oiorio/samples/B2GEDMNtuple_4.root",
#        "file:/tmp/oiorio/samples/B2GEDMNtuple_5.root",
#        "file:/tmp/oiorio/samples/B2GEDMNtuple_6.root",
#        "file:/tmp/oiorio/samples/B2GEDMNtuple_7.root",
        options.sample
#        options.sample
        )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag

process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2')  
#process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

### Rootplizer

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputLabel))

process.load("ttDM.treeDumper.topplusdmedmRootTreeMaker_cff")

process.DMTreesDumper.channelInfo.SingleElTriggers=cms.vstring(SingleElTriggers)
process.DMTreesDumper.channelInfo.SingleMuTriggers=cms.vstring(SingleMuTriggers)
process.DMTreesDumper.channelInfo.hadronicTriggers=cms.vstring(hadronTriggers)
process.DMTreesDumper.channelInfo.HadronPFHT800Triggers=cms.vstring(HadronPFHT800Triggers)
process.DMTreesDumper.channelInfo.HadronPFHT900Triggers=cms.vstring(HadronPFHT900Triggers)
process.DMTreesDumper.channelInfo.HadronPFJet450Triggers=cms.vstring(HadronPFJet450Triggers)

if options.addPartonInfo:
#    testMC = False
    testMC = True
    if not options.isData:                                                                                                                            
        process.DMTreesDumper.channelInfo.getPartonW=cms.untracked.bool(testMC)
        process.DMTreesDumper.channelInfo.getPartonTop=cms.untracked.bool(testMC)
        process.DMTreesDumper.channelInfo.doWReweighting=cms.untracked.bool(testMC)
        process.DMTreesDumper.channelInfo.doTopReweighting=cms.untracked.bool(testMC)

process.DMTreesDumper.lhes =cms.InputTag(options.lhesource)
process.DMTreesDumper.changeJECs = cms.untracked.bool(options.changeJECs)
process.DMTreesDumper.isData = cms.untracked.bool(options.isData)#This adds the L2L3Residuals
process.DMTreesDumper.applyRes = cms.untracked.bool(options.applyRes)#This adds the L2L3Residuals
process.DMTreesDumper.EraLabel = cms.untracked.string(options.EraLabel)
process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(True)
process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(True)

   
if options.channel == "ttbar":
    process.DMTreesDumper.getPartonTop  = cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.getParticleWZ  = cms.untracked.bool(True) 
if options.channel == "wzjets":
    print "channel is " + options.channel 
    process.DMTreesDumper.channelInfo.getPartonW  = cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.getParticleWZ  = cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.getWZFlavour  = cms.untracked.bool(True)
if options.isData:
    process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(False)
    process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(False)

if not options.isData:
    process.DMTreesDumper.getPartonTop  = cms.untracked.bool(True)
    process.DMTreesDumper.getPartonW  = cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(True) 
    process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(True)

if options.channel == "wprime":
    process.DMTreesDumper.lhes =cms.InputTag("source")
    process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(True)
    #    process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(False)
    process.DMTreesDumper.channelInfo.getPartonW=cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.getPartonTop=cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.doWReweighting=cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.doTopReweighting=cms.untracked.bool(True)


process.analysisPath = cms.Path(
    process.DMTreesDumper
    )

for p in process.DMTreesDumper.physicsObjects:
    if(p.prefix == cms.string("el")):
        p.variablesF.append(cms.InputTag("electrons","elvidVeto"))
        p.variablesF.append(cms.InputTag("electrons","elvidLoose"))
        p.variablesF.append(cms.InputTag("electrons","elvidMedium"))
        p.variablesF.append(cms.InputTag("electrons","elvidTight"))
        p.toSave.append("elvidVeto")
        p.toSave.append( "elvidLoose")
        p.toSave.append("elvidMedium")
        p.toSave.append("elvidTight")
            
