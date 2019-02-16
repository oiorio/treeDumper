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

options.register('sample',
                 #'/store/user/decosa/B2GAnaFW_80X_v2p4/BprimeBToHB_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/BprimeBToHB_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_v2p4/170213_085503/0000/B2GEDMNtuple_2.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/170104_182124/0000/B2GEDMNtuple_104.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016C-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_155142/0000/B2GEDMNtuple_10.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016H-PromptReco-v2_B2GAnaFW_80X_v2p4/161221_180745/0000/B2GEDMNtuple_10.root',
                 #'/store/user/grauco/B2GAnaFW_80X_V2p4/BprimeBToHB_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2/BprimeBToHB_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161227_111236/0000/B2GEDMNtuple_2.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p3/JetHT/Run2016B/JetHT/Run2016B-23Sep2016-v3_B2GAnaFW_80X_V2p3/161216_214635/0000/B2GEDMNtuple_1.root',
                 #'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1_B2GAnaFW_80X_V2p4/170124_202648/0000/B2GEDMNtuple_1.root',

#Recent files:
                 #2016 signal t'
                 #"/TprimeBToTZ_M-1400_LH_TuneCP5_13TeV-madgraph-pythia8/oiorio-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_B2GAnaFW_94X_V0_Sep_14-b20ae75985a288300defff97749d07bd/USER"
                 #'/store/user/decosa/EPS17/B2GAnaFW_80X_V3p1/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p1/170510_103545/0000/B2GEDMNtuple_1.root',
                 #2017 ntuples
                 #"/store/user/oiorio/ttDM/samples/2018/Sep/14Sep/B2GAnaFW_94X_V0_Sep_14/TprimeBToTZ_M-1400_LH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2/TprimeBToTZ_M-1400_LH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_B2GAnaFW_94X_V0_Sep_14/181003_184239/0000/B2GEDMNtuple_9.root",
                 "/store/user/lvigilan/Tprime/samples/94X_18Jan19/TprimeBToTZ_M-1400_Width-30p_LH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2/TprimeBToTZ_M-1400_Width-30p_LH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_94X_18Jan19/190118_170039/0000/B2GEDMNtuple_1.root",
                 #'file:B2GEDMNtuple_1.root',#l'ho preso dal mio TT
                 #'file:/eos/user/o/oiorio/samples/synch/mc_jecv4/B2GSynchMC.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')

options.register('version',
                 #'94X','94X_2016','80X','10X'
                 '94X',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'ntuple version (94X, 80X, 10X)')

options.register('outputLabel',
                 'treesTTJets.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

options.register('EraLabel',
                 'B',
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

options.register('mode',
                 'local',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Execution mode: local or crab')

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
print "version is: ", options.version, " is it data? ", options.isData

version = options.version #version: alias for commodity

if(options.isData):options.LHE = False

if(not options.isData): options.applyRes = False

#Trigger and filters setup:

if(version=="80X" or version=="94X" or version == "94X_2016"):
    l = ["singleTrigger"+str(s) for s in xrange(15)]
    l = l + ["trigger2"+str(s) for s in xrange(15)]

    SingleElTriggers = ["HLT_Ele27_eta2p1_WP75_Gsf"]
    SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WP75_Gsf_v"+str(s) for s in xrange(15)]
    
    SingleElControlTriggers = ["HLT_Ele23_WPLoose_Gsf"]##
    SingleElControlTriggers = SingleElControlTriggers + ["HLT_Ele23_WPLoose_Gsf_v"+str(s) for s in xrange(15)]

    SingleElHighPtTriggers =  ["HLT_Ele115_CaloIdVT_GsfTrkIdT"]
    SingleElHighPtTriggers = SingleElHighPtTriggers + ["HLT_Ele115_CaloIdVT_GsfTrkIdT_v"+str(s) for s in xrange(15)]

    if (version=="94X"):
        SingleElControlTriggers = ["HLT_Ele20_WPLoose_Gsf"]
        SingleElControlTriggers = SingleElControlTriggers + ["HLT_Ele20_WPLoose_Gsf_v"+str(s) for s in xrange(15)]
   
        SingleElTriggers = ["HLT_Ele27_WPTight_Gsf"]
        SingleElTriggers = SingleElTriggers + ["HLT_Ele27_WPTight_Gsf_v"+str(s) for s in xrange(15)]

        SingleElHighPtTriggers =  ["HLT_Ele115_CaloIdVT_GsfTrkIdT"]###
        SingleElHighPtTriggers = SingleElHighPtTriggers + ["HLT_Ele115_CaloIdVT_GsfTrkIdT_v"+str(s) for s in xrange(15)]

    SingleMuTriggers = ["HLT_IsoMu22_v"+str(s) for s in range(1,5)]
    SingleMuTriggers = SingleMuTriggers + ["HLT_IsoMu24_v"+str(s) for s in range(1,5)]
    SingleMuTriggers = SingleMuTriggers + ["HLT_IsoTkMu24_v"+str(s) for s in range(1,5)]
    
    SingleMuHighPtTriggers = ["HLT_Mu50","HLT_TkMu50"]
    SingleMuHighPtTriggers = SingleMuHighPtTriggers + ["HLT_Mu50_v"+str(s) for s in xrange(15)]
    SingleMuHighPtTriggers = SingleMuHighPtTriggers + ["HLT_TkMu50_v"+str(s) for s in xrange(15)]

    SingleMuControlTriggers = ["HLT_Mu50","HLT_TkMu50"]
    SingleMuControlTriggers = SingleMuControlTriggers + ["HLT_Mu50_v"+str(s) for s in xrange(15)]
    

    PhotonTriggers = [""]
    
    HadronPFHT900Triggers = ["HLT_PFHT900"]
    HadronPFHT900Triggers = HadronPFHT900Triggers + ["HLT_PFHT900_v"+str(s) for s in xrange(15)]

    HadronPFHT800Triggers = ["HLT_PFHT800"]
    HadronPFHT800Triggers =HadronPFHT800Triggers + ["HLT_PFHT800_v"+str(s) for s in xrange(15)]

    HadronPFJet450Triggers = ["HLT_PFJet450"]
    HadronPFJet450Triggers =HadronPFJet450Triggers + ["HLT_PFJet450_v"+str(s) for s in xrange(15)]
    metTriggers = ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"] #metTriggers = metTriggers + ["HLT_PFMET90_PFMHT190_IDTight_v"+str(s) for s in xrange(15)] 
    #metTriggers = metTriggers + ["HLT_PFMET100_PFMHT100_IDTight_v"+str(s) for s in xrange(15)]
    metTriggers = metTriggers + ["HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v"+str(s) for s in xrange(15)]
    metTriggers = metTriggers + ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]       

    if(options.isData):

        SingleElTriggers = ["HLT_Ele27_eta2p1_WPLoose_Gsf"]
        SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WPLoose_Gsf_v"+str(s) for s in xrange(15)]

        SingleElControlTriggers = ["HLT_Ele23_WPLoose_Gsf"]
        SingleElControlTriggers = SingleElControlTriggers + ["HLT_Ele23_WPLoose_Gsf_v"+str(s) for s in xrange(15)]
                                    
        SingleElHighPtTriggers =  ["HLT_Ele115_CaloIdVT_GsfTrkIdT"]
        SingleElHighPtTriggers = SingleElHighPtTriggers + ["HLT_Ele115_CaloIdVT_GsfTrkIdT_v"+str(s) for s in xrange(15)]

        if (version=="94X"):
            SingleElTriggers = ["HLT_Ele27_WPTight_Gsf"]
            SingleElTriggers = SingleElTriggers + ["HLT_Ele27_WPTight_Gsf_v"+str(s) for s in xrange(15)]
            
            SingleElControlTriggers = ["HLT_Ele20_WPLoose_Gsf"]
            SingleElControlTriggers = SingleElControlTriggers + ["HLT_Ele20_WPLoose_Gsf_v"+str(s) for s in xrange(15)]
                        
            SingleElHighPtTriggers =  ["HLT_Ele115_CaloIdVT_GsfTrkIdT"]
            SingleElHighPtTriggers = SingleElHighPtTriggers + ["HLT_Ele115_CaloIdVT_GsfTrkIdT_v"+str(s) for s in xrange(15)]
        
        SingleMuTriggers = ["HLT_Mu50","HLT_TkMu50"]
        SingleMuTriggers = SingleMuTriggers + ["HLT_Mu50_v"+str(s) for s in xrange(15)]
        SingleMuTriggers = SingleMuTriggers + ["HLT_TkMu50_v"+str(s) for s in xrange(15)]

        SingleMuHighPtTriggers = ["HLT_Mu50","HLT_TkMu50"]
        SingleMuHighPtTriggers = SingleMuHighPtTriggers + ["HLT_Mu50_v"+str(s) for s in xrange(15)]
        SingleMuHighPtTriggers = SingleMuHighPtTriggers + ["HLT_TkMu50_v"+str(s) for s in xrange(15)]
        
        HadronPFHT900Triggers = ["HLT_PFHT900"]
        HadronPFHT900Triggers = HadronPFHT900Triggers + ["HLT_PFHT900_v"+str(s) for s in xrange(15)]
    
        HadronPFHT800Triggers = ["HLT_PFHT800"]
        HadronPFHT800Triggers =HadronPFHT800Triggers + ["HLT_PFHT800_v"+str(s) for s in xrange(15)]
    
        HadronPFJet450Triggers = ["HLT_PFJet450"]
        HadronPFJet450Triggers =HadronPFJet450Triggers + ["HLT_PFJet450_v"+str(s) for s in xrange(15)]

        metTriggers = ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"]

        metTriggers = metTriggers + ["HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v"+str(s) for s in xrange(15)]
        metTriggers = metTriggers + ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]

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

process.load("TreeFWK.treeDumper.topplusdmedmRootTreeMaker_cff")

process.DMTreesDumper.channelInfo.SingleElTriggers=cms.vstring(SingleElTriggers)
process.DMTreesDumper.channelInfo.SingleMuTriggers=cms.vstring(SingleMuTriggers)
process.DMTreesDumper.channelInfo.SingleElControlTriggers=cms.vstring(SingleElControlTriggers)
process.DMTreesDumper.channelInfo.SingleMuControlTriggers=cms.vstring(SingleMuControlTriggers)
process.DMTreesDumper.channelInfo.SingleElHighPtTriggers=cms.vstring(SingleElHighPtTriggers)
process.DMTreesDumper.channelInfo.SingleMuHighPtTriggers=cms.vstring(SingleMuHighPtTriggers)
process.DMTreesDumper.channelInfo.MetTriggers=cms.vstring(metTriggers)
process.DMTreesDumper.channelInfo.MetControlTriggers=cms.vstring(metTriggers)
process.DMTreesDumper.channelInfo.HadronicHHTTriggers=cms.vstring(HadronPFHT800Triggers)
process.DMTreesDumper.channelInfo.HadronicHHTControlTriggers=cms.vstring(HadronPFHT900Triggers)
process.DMTreesDumper.channelInfo.SingleJetTriggers=cms.vstring(HadronPFJet450Triggers)
process.DMTreesDumper.channelInfo.SingleJetControlTriggers=cms.vstring(HadronPFJet450Triggers)
process.DMTreesDumper.channelInfo.SingleJetSubstructureTriggers=cms.vstring(HadronPFJet450Triggers)
process.DMTreesDumper.channelInfo.SingleJetSubstructureControlTriggers=cms.vstring(HadronPFJet450Triggers)

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


#if addAK8CHS:
#    DMTreesDumper.physicsObjects.append(jet8CHS)
#    DMTreesDumper.physicsObjects.append(subjetCHS)
#if addAK8PUPPI:
#    DMTreesDumper.physicsObjects.append(jet8Puppi)
#    DMTreesDumper.physicsObjects.append(subjetPuppi)


#Setup categories and systematics:
#SingleTop:

catMuST = ["Tight","TightAntiIso","Loose"]
catElST = ["Tight","TightAntiIso","Veto","Antiveto"]
catJetST = ["Tight"]
catMetST = ["CorrT1"]

scanMuST = ["Iso04_0p06_LE","Iso04_0p15_LE","Iso04_0p06_GE","Iso04_0p15_GE"]
scanEl = []
scanJet = ["CorrPt_20","CorrPt_40"]
sysMu = []
sysEl = []
sysJet = ["JESUp","JESDown","JERUp","JERDown"]

#Setup JECs and Jet collections:

sj = "subjetsAK8CHS"
sjpref = "subjetAK8CHS"
sjpup = "subjetsAK8Puppi"
sjpuppref = "subjetAK8Puppi"
jetak8label = cms.string("jetsAK8CHS")
subjetak8label = cms.string("subjetsAK8CHS")
jetak8puplabel = cms.string("jetsAK8Puppi")
subjetak8puplabel = cms.string("subjetsAK8Puppi")


addAK8CHS=True
addAK8Puppi=True

jecFolder="JECs/"
if(options.mode=="crab"):
    jecFolder="./"

if(version=='80X'):
    process.DMTreesDumper.era = cms.untracked.string("2016_80X")
    process.DMTreesDumper.prefixLabelData = cms.untracked.string("Summer16_23Sep2016V4")
    process.DMTreesDumper.prefixLabelMC = cms.untracked.string("Summer16_23Sep2016V4")
    process.DMTreesDumper.postfixLabelData = cms.untracked.string("V4_DATA")
    process.DMTreesDumper.postfixLabelMC = cms.untracked.string("_MC")
    process.DMTreesDumper.jetType = cms.untracked.string("AK4PFchs")
    process.DMTreesDumper.jetType8 = cms.untracked.string("AK8PFchs")
    process.DMTreesDumper.boostedTopsLabel = jetak8label
    process.DMTreesDumper.boostedTopsSubjetsLabel = subjetak8label    

if(version=='94X'):
    process.DMTreesDumper.era = cms.untracked.string("2017_94X")
    process.DMTreesDumper.prefixLabelData = cms.untracked.string(jecFolder+"/Fall17_17Nov2017")
    process.DMTreesDumper.prefixLabelMC = cms.untracked.string(jecFolder+"/Fall17_17Nov2017")
    process.DMTreesDumper.postfixLabelData = cms.untracked.string("_V32_DATA")
    process.DMTreesDumper.postfixLabelMC = cms.untracked.string("_V32_MC")
    process.DMTreesDumper.jetType = cms.untracked.string("AK4PFchs")
    process.DMTreesDumper.jetType8 = cms.untracked.string("AK8PFPuppi")
    process.DMTreesDumper.boostedTopsLabel = jetak8puplabel
    process.DMTreesDumper.boostedTopsSubjetsLabel = subjetak8puplabel    
    process.DMTreesDumper.prefiringWeight = cms.InputTag("prefiringweight","NonPrefiringProb")
    process.DMTreesDumper.prefiringWeightUp = cms.InputTag("prefiringweight","NonPrefiringProbUp")
    process.DMTreesDumper.prefiringWeightDown = cms.InputTag("prefiringweight","NonPrefiringProbDown")

if(version=='94X_2016'):
    process.DMTreesDumper.era = cms.untracked.string("2016_94X")
    process.DMTreesDumper.prefixLabelData = cms.untracked.string(jecFolder+"/Summer16_07Aug2017")
    process.DMTreesDumper.prefixLabelMC = cms.untracked.string(jecFolder+"/Summer16_23Sep2016V4")
    process.DMTreesDumper.postfixLabelData = cms.untracked.string("V4_DATA")
    process.DMTreesDumper.postfixLabelMC = cms.untracked.string("_MC")
    process.DMTreesDumper.jetType = cms.untracked.string("AK4PFchs")
    process.DMTreesDumper.jetType8 = cms.untracked.string("AK8PFPuppi")
    process.DMTreesDumper.boostedTopsLabel = jetak8puplabel
    process.DMTreesDumper.boostedTopsSubjetsLabel = subjetak8puplabel    
    process.DMTreesDumper.prefiringWeight = cms.InputTag("prefiringweight","NonPrefiringProb")
    process.DMTreesDumper.prefiringWeightUp = cms.InputTag("prefiringweight","NonPrefiringProbUp")
    process.DMTreesDumper.prefiringWeightDown = cms.InputTag("prefiringweight","NonPrefiringProbDown")

    
if(addAK8CHS):
   process.DMTreesDumper.physicsObjects.append(process.jet8CHS)
   process.DMTreesDumper.physicsObjects.append(process.subjetCHS)

if(addAK8Puppi):
   process.DMTreesDumper.physicsObjects.append(process.jet8Puppi)
   process.DMTreesDumper.physicsObjects.append(process.subjetPuppi)


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
    process.DMTreesDumper.channelInfo.getPartonW=cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.getPartonTop=cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.doWReweighting=cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.doTopReweighting=cms.untracked.bool(True)


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
            
