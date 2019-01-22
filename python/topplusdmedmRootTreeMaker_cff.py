import FWCore.ParameterSet.Config as cms
import copy
#from ttDM.treeDumper.topplusdmTrees_cfg import isData

leptonssize = cms.untracked.int32(6)
jetssize = cms.untracked.int32(20)
genpartsize = cms.untracked.int32(100)

pholabel = cms.string("photons")
mulabel = cms.string("muons")
elelabel = cms.string("electrons")
jetlabel = cms.string("jetsAK4CHS")
jetak8label = cms.string("jetsAK8CHS")
subjetak8label = cms.string("subjetsAK8CHS")

jetak8puplabel = cms.string("jetsAK8Puppi")
subjetak8puplabel = cms.string("subjetsAK8Puppi")

metlabel = cms.string("metFull")
eventlabel = cms.string("eventShapePFVars")
eventJetlabel = cms.string("eventShapePFJetVars")
genlabel = cms.string("genPart")



#Systematics:
#systsToSave = ["noSyst","jes__up","jes__down"]
#systsToSave = ["noSyst","jer__up","jer__down"]
#systsToSave = ["noSyst","jer__up","jer__down","jes__up","jes__down"]
systsToSave = ["noSyst"]


#metFilters = ["Flag_CSCTightHaloFilter","Flag_goodVertices", "Flag_eeBadScFilter"]
metFilters=["Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_globalTightHalo2016Filter"]

catMu = ["Tight","Loose"]
catEl = ["Tight","Veto"]
catJet = ["Tight"]
catMet = ["CorrT1"]

scanMu = []
scanEl = []
scanJet = []
#scanJet = ["CorrPt_30"]

sysMu = []
sysEl = []
#sysJet = [""]
sysJet = ["JESUp","JESDown","JERUp","JERDown"]
#sysJet = ["JESUp","JESDown"]

#SingleTop
#catMu = ["Tight","TightAntiIso","Loose"]
#catEl = ["Tight","TightAntiIso","Veto","Antiveto"]
#catJet = ["Tight"]
#catMet = ["CorrT1"]
#catMet = []
#scanMu = ["Iso04_0p06_LE","Iso04_0p15_LE","Iso04_0p06_GE","Iso04_0p15_GE"]
#scanEl = []
#scanJet = ["CorrPt_20","CorrPt_40"]
#sysMu = [""]
#sysEl = [""]
#sysJet = ["JESUp","JESDown","JERUp","JERDown"]


#Setup categories and systematics:
#SingleTop:
SingleTopSetup=False
if SingleTopSetup:
    catMu = ["Tight","TightAntiIso","Loose"]
    catEl = ["Tight","TightAntiIso","Veto","Antiveto"]
    catJet = ["Tight"]
    catMet = ["CorrT1"]
    
    scanMu = ["Iso04_0p06_LE","Iso04_0p15_LE","Iso04_0p06_GE","Iso04_0p15_GE"]
    scanEl = []
    scanJet = ["CorrPt_20","CorrPt_40"]
    
    sysMu = [""]
    sysEl = [""]
    sysJet = ["JESUp","JESDown","JERUp","JERDown"]


cutOnTriggers = False
doPreselectionCuts = False

#What to use for jets/other variables
saveBase = cms.untracked.bool(False)
saveNoCategory = True
j= "jetsAK4CHS"
jpref= "jetAK4CHS"

sj = "subjetsAK8CHS"
sjpref = "subjetAK8CHS"

sjpup = "subjetsAK8Puppi"
sjpuppref = "subjetAK8Puppi"


#sj = "subjetsCmsTopTag"
#sjpref = "subjetsCmsTopTag"
subjetak8label = cms.string(sj)

#Initializing the analyzer
DMTreesDumper = cms.EDAnalyzer(
    'DMAnalysisTreeMaker',
    lhes = cms.InputTag('source'),
    genprod = cms.InputTag('generator'),
    muLabel = mulabel,
    eleLabel = elelabel,
    jetsLabel = jetlabel,
    photonLabel = pholabel,
    boostedTopsLabel = jetak8label,
    boostedTopsSubjetsLabel = subjetak8label,
    metLabel = metlabel,
    genLabel = genlabel,
    eventLabel = eventlabel,
    #resolutionsFile = cms.string('Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt'),
    #scaleFactorsFile = cms.string('Spring16_25nsV10_MC_SF_AK8PFchs.txt'),
    physicsObjects = cms.VPSet(
        cms.PSet(
            label = metlabel,
            prefix = cms.string("metFull"),
            maxInstances =  cms.untracked.int32(1),
            saveBaseVariables = cms.untracked.bool(True),
            categories = cms.vstring(catMet),
            scanCuts = cms.vstring(),
            systCats = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(
                cms.InputTag("metFull","metFulluncorPt"),
                cms.InputTag("metFull","metFulluncorPhi"),
                cms.InputTag("metFull","metFullPt"),
                cms.InputTag("metFull","metFullPhi"),
                cms.InputTag("metFull","metFullPx"),
                cms.InputTag("metFull","metFullPy"),
                ),
            variablesI = cms.VInputTag(),
            singleD = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            toSave = cms.vstring(),
            )
        ),
    version = cms.untracked.string("2017_94X"),
    doPreselection = cms.untracked.bool(doPreselectionCuts), 
  
    ak8jetvSubjetIndex0 = cms.InputTag( "jetsAK8", "jetAK8vSubjetIndex0"),
    ak8jetvSubjetIndex1 = cms.InputTag( "jetsAK8", "jetAK8vSubjetIndex1"),

    #trigger part:
    useTriggers = cms.untracked.bool(True),
    cutOnTriggers = cms.untracked.bool(cutOnTriggers),
    triggerBits = cms.InputTag("TriggerUserData","triggerBitTree"),
    triggerNames = cms.InputTag("TriggerUserData","triggerNameTree"),
    triggerPrescales = cms.InputTag("TriggerUserData","triggerPrescaleTree"),
    #met filters
    metBits = cms.InputTag("METUserData","triggerBitTree"),
    metNames = cms.InputTag("METUserData","triggerNameTree"),
    #lumi,run,number
    lumiBlock = cms.InputTag("eventInfo","evtInfoLumiBlock"),
    runNumber = cms.InputTag("eventInfo","evtInfoRunNumber"),
    eventNumber = cms.InputTag("eventInfo","evtInfoEventNumber"),
    #HBHE
    HBHEIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),

    #Filters
    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter", ""),
    BadPFMuonFilter =cms.InputTag("BadPFMuonFilter", ""),
    
    jetKeysAK4CHS = cms.InputTag("jetKeysAK4CHS", ""),
    muonKeys = cms.InputTag("muonKeys", ""),

    #vertex
    vertexZ =  cms.InputTag("vertexInfo","z"),
    vertexChi2 =  cms.InputTag("vertexInfo","chi"),
    vertexNdof =  cms.InputTag("vertexInfo","ndof"),
    vertexRho =  cms.InputTag("vertexInfo","rho"),

    #resolved top part:
    doResolvedTopHad=cms.untracked.bool(False),
    doResolvedTopSemiLep=cms.untracked.bool(False),
    #cuts for the jet scan
    jetScanCuts=cms.vdouble(30), #Note: the order is important, as the jet collection with the first cut is used for the definition of mt2w.

    #corrections
    prefixLabelData = cms.untracked.string("Summer16_23Sep2016V4"),
    prefixLabelMC = cms.untracked.string("Summer16_23Sep2016V4"),
    postfixLabelData = cms.untracked.string("V4_DATA"),
    postfixLabelMC = cms.untracked.string("_MV"),
    jetType = cms.untracked.string("AK4PFchs"),
    jetType8 = cms.untracked.string("AK8PFchs"),

    #Systematics trees to produce. Include:  
    systematics = cms.vstring(systsToSave), 

    channelInfo = cms.PSet(
        channel = cms.string("ttDM"),#Name of the channel, to use in the trees
        crossSection = cms.double(1),#Cross section in pb
        originalEvents = cms.double(1),#Number of events in the MC
        SingleElTriggers = cms.vstring(""),
        SingleMuTriggers = cms.vstring(""),
        SingleElHighPtTriggers = cms.vstring(""),
        SingleMuHighPtTriggers = cms.vstring(""),
        SingleElControlTriggers = cms.vstring(""),
        SingleMuControlTriggers = cms.vstring(""),
        PhotonTriggers = cms.vstring(""),
        HadronicHTTriggers = cms.vstring(""),
        HadronicHTControlTriggers = cms.vstring(""),
        MetTriggers = cms.vstring(""),
        MetControlTriggers = cms.vstring(""),
        SingleJetTriggers = cms.vstring(""),
        SingleJetControlTriggers = cms.vstring(""),
        SingleJetSubstructureTriggers = cms.vstring(""),
        SingleJetSubstructureControlTriggers = cms.vstring(""),
        metFilters = cms.vstring(metFilters),
        useLHE = cms.untracked.bool(False),#Whether one uses the weights from the LHE in order to get scale uncertainties
        useLHEWeights = cms.untracked.bool(False),#Whether one uses the weights from the LHE in order to get scale uncertaintiesxb
        addLHAPDFWeights = cms.untracked.bool(False), #Whether to add the PDF for uncertainty evaluation (time consuming)
        maxWeights = cms.untracked.int32(111),
        )
    )


#Now taking the other input objects:
DMTreesDumper.physicsObjects.append(  
    cms.PSet(
        label = mulabel,
        prefix = cms.string("mu"),
        maxInstances = leptonssize,
        saveBaseVariables = saveBase,
        saveNoCat = cms.untracked.bool(saveNoCategory),
        #        categories = cms.vstring(),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag("muons","muE"),
            cms.InputTag("muons","muPt"),
            #cms.InputTag("muons","muMass"),
            cms.InputTag("muons","muEta"),
            cms.InputTag("muons","muPhi"),
            cms.InputTag("muons","muCharge"),
            cms.InputTag("muons","muIsLooseMuon"),
            cms.InputTag("muons","muIsSoftMuon"),
            cms.InputTag("muons","muIsMediumMuon"),
            cms.InputTag("muons","muIsTightMuon"),
            cms.InputTag("muons","muIsHighPtMuon"),
            cms.InputTag("muons","muDB"),
            cms.InputTag("muons","muDBerr"),
            cms.InputTag("muons","muDz"),
            cms.InputTag("muons","muDxy"),
            cms.InputTag("muons","muGenMuonCharge"),
            cms.InputTag("muons","muGenMuonEta"),
            cms.InputTag("muons","muGenMuonPt"),
            cms.InputTag("muons","muGenMuonE"),
            cms.InputTag("muons","muGenMuonPhi"),
            cms.InputTag("muons","muGlbTrkNormChi2"),
            cms.InputTag("muons","muInTrkNormChi2"),
            cms.InputTag("muons","muIsGlobalMuon"),
            cms.InputTag("muons","muIsPFMuon"),
            cms.InputTag("muons","muIsTrackerMuon"),
            cms.InputTag("muons","muIso04"),
            cms.InputTag("muons","muNumberMatchedStations"),
            cms.InputTag("muons","muNumberOfPixelLayers"),
            cms.InputTag("muons","muNumberOfValidTrackerHits"),
            cms.InputTag("muons","muNumberTrackerLayers"),
            cms.InputTag("muons","muNumberValidMuonHits"),
            cms.InputTag("muons","muNumberValidPixelHits"),
            cms.InputTag("muons","muSumChargedHadronPt"),
            cms.InputTag("muons","muSumNeutralHadronPt"),
            cms.InputTag("muons","muSumPUPt"),
            cms.InputTag("muons","muSumPhotonPt"),            
            ),
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        scanCuts = cms.vstring(scanMu),
        categories = cms.vstring(catMu),
        systCats = cms.vstring(sysMu),
        #        categories = cms.vstring("Tight","Loose"),
        toSave = cms.vstring("muE","muPt","muEta","muPhi","muIso04","muCharge","muIsMediumMuon","muIsTightMuon","muIsLooseMuon","muIsHighPtMuon","allExtra"),
        )
    )

addPho = False
if addPho:
    DMTreesDumper.physicsObjects.append( 
    cms.PSet(
        label = pholabel,
        prefix = cms.string("pho"),
        maxInstances =  leptonssize,
        saveBaseVariables = cms.untracked.bool(True),
        categories = cms.vstring(),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag("photons","phoEta"),
            cms.InputTag("photons","phoPt"),
            cms.InputTag("photons","phoHoverE"),
            cms.InputTag("photons","phoSigmaIEtaIEta"),
            cms.InputTag("photons","phoNeutralHadronIso"),
            cms.InputTag("photons","phoChargedHadronIso"),
            cms.InputTag("photons","phoPhotonIso"),
            cms.InputTag("photons","phoHasPixelSeed"),
            cms.InputTag("photons","phoPassLooseID"),
            cms.InputTag("photons","phoPassMediumID"),
            cms.InputTag("photons","phoPassTightID"),
            ),
        variablesI = cms.VInputTag(
            ),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        scanCuts = cms.vstring(),
        systCats = cms.vstring(),
        toSave = cms.vstring("phoPt","phoHoverE","phoSigmaIEtaIEta","phoNeutralHadronIso","phoNeutralHadronIsoEAcorrected","phoChargedHadronIso","phoChargedHadronIsoEAcorrected","phoPhotonIso","phoPhotonIsoEAcorrected","phoHasPixelSeed","phoPassLooseID","phoPassMediumID","phoPassTightID","allExtra"),
        )
    )
DMTreesDumper.physicsObjects.append( 
    cms.PSet(
        label = elelabel,
        prefix = cms.string("el"),
        maxInstances = leptonssize,
        saveBaseVariables = saveBase,
        saveNoCat = cms.untracked.bool(saveNoCategory),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag("electrons","elE"),
            cms.InputTag("electrons","elPt"),
            cms.InputTag("electrons","elEta"),
            cms.InputTag("electrons","elPhi"),
            cms.InputTag("electrons","elCharge"),
            cms.InputTag("electrons","elDz"),
            cms.InputTag("electrons","elEta"),
            cms.InputTag("electrons","elSCEta"),
            cms.InputTag("electrons","elHoE"),
            cms.InputTag("electrons","elIso03"),
            cms.InputTag("electrons","eldEtaIn"),
            cms.InputTag("electrons","eldPhiIn"),
            cms.InputTag("electrons","elmissHits"),
            cms.InputTag("electrons","elfull5x5siee"),
            cms.InputTag("electrons","elooEmooP"),
            cms.InputTag("electrons","elhasMatchedConVeto"),
            cms.InputTag("electrons","elvidVeto"),
            cms.InputTag("electrons","elvidTight"),
            cms.InputTag("electrons","elvidMedium"),
        ),
        variablesI = cms.VInputTag(
            ),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        categories = cms.vstring(catEl),
        scanCuts = cms.vstring(scanEl),
        systCats = cms.vstring(sysEl),
        toSave = cms.vstring("elE","elPt","elEta","elPhi","elIso03","elisTight","elisMedium","elisLoose","elisVeto","elscEta","allExtra"),
        )
    ) 



DMTreesDumper.physicsObjects.append( 
    cms.PSet(
        label = eventJetlabel,
        prefix = cms.string(""),
        maxInstances =  cms.untracked.int32(1),
        saveBaseVariables = cms.untracked.bool(True),
        categories = cms.vstring(),
        scanCuts = cms.vstring(),
        systCats = cms.vstring(),
        variablesF = cms.VInputTag(),
        variablesD = cms.VInputTag(),
        variablesI = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        singleD = cms.VInputTag(
            #cms.InputTag("eventShapePFJetVars","isotropy"),
            #cms.InputTag("eventShapePFJetVars","C"),
            #cms.InputTag("eventShapePFJetVars","D"),
            #cms.InputTag("eventShapePFJetVars","aplanarity"),
            #cms.InputTag("eventShapePFJetVars","circularity"),
            #cms.InputTag("eventShapePFJetVars","sphericity"),
            #cms.InputTag("eventShapePFJetVars","thrust"),
            #cms.InputTag("eventShapePFJetVars","thrust"),
            ),
        toSave = cms.vstring(),
        )
    )

DMTreesDumper.physicsObjects.append( 
    cms.PSet(
        label = jetlabel,
        prefix = cms.string(jpref),
        maxInstances = jetssize,
        saveBaseVariables = saveBase,
        saveNoCat = cms.untracked.bool(saveNoCategory),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag(j,jpref+"E"),
            cms.InputTag(j,jpref+"Pt"),
            #cms.InputTag(j,jpref+"CorrPt"),
            #cms.InputTag(j,jpref+"Mass"),
            cms.InputTag(j,jpref+"Eta"),
            cms.InputTag(j,jpref+"Phi"),
            cms.InputTag(j,jpref+"QGL"),
            #cms.InputTag(j,jpref+"mult"),
            #cms.InputTag(j,jpref+"axis2"),
            #cms.InputTag(j,jpref+"ptD"),
            cms.InputTag(j,jpref+"PartonFlavour"),
            cms.InputTag(j,jpref+"Phi"),
            cms.InputTag(j,jpref+"CSVv2"),
            cms.InputTag(j,jpref+"CMVAv2"),
            cms.InputTag(j,jpref+"DeepCSV"),
            cms.InputTag(j,jpref+"DeepCSVProbb"),
            cms.InputTag(j,jpref+"DeepCSVBvsAll"),
            cms.InputTag(j,jpref+"Charge"),
            cms.InputTag(j,jpref+"ElectronEnergy"),
            cms.InputTag(j,jpref+"GenJetCharge"),
            cms.InputTag(j,jpref+"GenJetE"),
            cms.InputTag(j,jpref+"GenJetEta"),
            cms.InputTag(j,jpref+"GenJetPhi"),
            cms.InputTag(j,jpref+"GenJetPt"),
            cms.InputTag(j,jpref+"GenPartonCharge"),
            cms.InputTag(j,jpref+"GenPartonE"),
            cms.InputTag(j,jpref+"GenPartonEta"),
            cms.InputTag(j,jpref+"GenPartonPhi"),
            cms.InputTag(j,jpref+"GenPartonPt"),
            cms.InputTag(j,jpref+"HadronFlavour"),
            #cms.InputTag(j,jpref+"SmearedE"),
            #cms.InputTag(j,jpref+"SmearedPt"),
            #cms.InputTag(j,jpref+"SmearedPEta"),
            #cms.InputTag(j,jpref+"SmearedPhi"),
            #cms.InputTag(j,jpref+"Y"),
            cms.InputTag(j,jpref+"chargedEmEnergyFrac"),
            cms.InputTag(j,jpref+"chargedHadronEnergyFrac"),
            cms.InputTag(j,jpref+"neutralEmEnergyFrac"),
            cms.InputTag(j,jpref+"neutralHadronEnergyFrac"),
            cms.InputTag(j,jpref+"electronMultiplicity"),
            cms.InputTag(j,jpref+"muonMultiplicity"),
            cms.InputTag(j,jpref+"neutralMultiplicity"),
            cms.InputTag(j,jpref+"photonMultiplicity"),
            cms.InputTag(j,jpref+"chargedMultiplicity"),
            cms.InputTag(j,jpref+"neutralEmEnergy"),
            cms.InputTag(j,jpref+"jecFactor0"),
            cms.InputTag(j,jpref+"jetArea"),
            ),
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        categories = cms.vstring(catJet),
        scanCuts = cms.vstring(scanJet),
        systCats = cms.vstring(sysJet),
        #categories = cms.vstring("Tight"),
        toSave = cms.vstring(jpref+"E",jpref+"Pt",jpref+"Eta",jpref+"Phi",jpref+"GenJetE", jpref+"GenJetPhi", jpref+"GenJetPt",jpref+"GenJetEta",jpref+"CSVv2",jpref+"DeepCSV",jpref+"PartonFlavour","allExtra"),
        )
    )

jet8CHS=  cms.PSet(
        label = jetak8label,
        prefix = cms.string("jetAK8CHS"),
        maxInstances = cms.untracked.int32(10),
        categories = cms.vstring(),
        scanCuts = cms.vstring(),
        systCats = cms.vstring(sysJet),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetE"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetPt"),
            #cms.InputTag("jetsAK8CHS","jetAK8CHSMass"),                                                                                                 
            cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetEta"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetPhi"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSE"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSPt"),
            #cms.InputTag("jetsAK8CHS","jetAK8CHSMass"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSEta"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSPhi"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSPartonFlavour"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSHadronFlavour"),
            #cms.InputTag("jetsAK8CHS","jetAK8CHSMass"),
            #cms.InputTag("jetsAK8CHS","jetAK8CHSminmass"),
            #cms.InputTag("jetsAK8CHS","jetAK8CHSnSubJets"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSvSubjetIndex0"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSvSubjetIndex1"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSprunedMassCHS"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSsoftDropMassCHS"),
            cms.InputTag("jetsAK8CHS","jetAK8CHStau1CHS"),
            cms.InputTag("jetsAK8CHS","jetAK8CHStau2CHS"),
            cms.InputTag("jetsAK8CHS","jetAK8CHStau3CHS"),
            #cms.InputTag("jetsAK8CHS","jetAK8CHStopMass"),
            cms.InputTag("jetsAK8CHS","jetAK8CHStrimmedMassCHS"),
            #cms.InputTag("jetsAK8CHS","jetAK8CHSwMass"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSjecFactor0"),
            cms.InputTag("jetsAK8CHS","jetAK8CHSjetArea"),
            ),
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        toSave = cms.vstring("jetAK8CHSGenJetE","jetAK8CHSGenJetPt","jetAK8CHSGenJetEta","jetAK8CHSGenJetPhi","jetAK8CHSE","jetAK8CHSPt","jetAK8CHSEta","jetAK8CHSPhi", "jetAK8CHSprunedMassCHS", "jetAK8CHStau1CHS", "jetAK8CHStau2CHS", "jetAK8CHStau3CHS",  "jetAK8CHStrimmedMassCHS", "jetAK8CHSPartonFlavour","jetAK8CHSHadronFlavour", "allExtra"),
        )




subjetCHS= cms.PSet(
        label = subjetak8label,
        prefix = cms.string(sjpref),
        maxInstances = cms.untracked.int32(10),
        categories = cms.vstring(),
        scanCuts = cms.vstring(),
        systCats = cms.vstring(sysJet),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag(sj,sjpref+"E"),
            cms.InputTag(sj,sjpref+"Pt"),
            #cms.InputTag(sj,sjpref+"Mass"),
            cms.InputTag(sj,sjpref+"Eta"),
            cms.InputTag(sj,sjpref+"Phi"),
            cms.InputTag(sj,sjpref+"PartonFlavour"),
            cms.InputTag(sj,sjpref+"CSVv2"),
            cms.InputTag(sj,sjpref+"DeepCSV"),
            ),                         
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        toSave = cms.vstring(sjpref+"E",sjpref+"Pt",sjpref+"Eta",sjpref+"Phi", "allExtra"),
        )
    

jet8Puppi= cms.PSet(
        label = jetak8puplabel,
        prefix = cms.string("jetAK8Puppi"),
        maxInstances = cms.untracked.int32(10),
        categories = cms.vstring(),
        scanCuts = cms.vstring(),
        systCats = cms.vstring(sysJet),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetE"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetPt"),
            #cms.InputTag("jetsAK8Puppi","jetAK8PuppiMass"),                                                                                                 
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetEta"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetPhi"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiE"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiPt"),
            #cms.InputTag("jetsAK8Puppi","jetAK8PuppiMass"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiEta"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiPhi"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiPartonFlavour"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiHadronFlavour"),
            #cms.InputTag("jetsAK8Puppi","jetAK8PuppiMass"),
            #cms.InputTag("jetsAK8Puppi","jetAK8Puppiminmass"),
            #cms.InputTag("jetsAK8Puppi","jetAK8PuppinSubJets"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppivSubjetIndex0"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppivSubjetIndex1"),
#            cms.InputTag("jetsAK8Puppi","jetAK8PuppivSubjetIndex2"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppiprunedMass"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppisoftDropMass"),
            cms.InputTag("jetsAK8Puppi","jetAK8Puppitau1"),
            cms.InputTag("jetsAK8Puppi","jetAK8Puppitau2"),
            cms.InputTag("jetsAK8Puppi","jetAK8Puppitau3"),
            #cms.InputTag("jetsAK8Puppi","jetAK8PuppitopMass"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppitrimmedMass"),
            #cms.InputTag("jetsAK8Puppi","jetAK8PuppiwMass"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppijecFactor0"),
            cms.InputTag("jetsAK8Puppi","jetAK8PuppijetArea"),
            ),
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        toSave = cms.vstring("jetAK8PuppiGenJetE","jetAK8PuppiGenJetPt","jetAK8PuppiGenJetEta","jetAK8PuppiGenJetPhi","jetAK8PuppiE","jetAK8PuppiPt","jetAK8PuppiEta","jetAK8PuppiPhi", "jetAK8PuppiprunedMassCHS", "jetAK8Puppitau1CHS", "jetAK8Puppitau2CHS", "jetAK8Puppitau3CHS",  "jetAK8PuppiPartonFlavour","jetAK8PuppiHadronFlavour","jetAK8PuppitrimmedMassCHS", "allExtra"),
        )


subjetPuppi=cms.PSet(
        label = subjetak8puplabel,
        prefix = cms.string(sjpuppref),
        maxInstances = cms.untracked.int32(10),
        categories = cms.vstring(),
        scanCuts = cms.vstring(),
        systCats = cms.vstring(),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag(sjpup,sjpuppref+"E"),
            cms.InputTag(sjpup,sjpuppref+"Pt"),
            #cms.InputTag(sjpup,sjpuppref+"Mass"),
            cms.InputTag(sjpup,sjpuppref+"Eta"),
            cms.InputTag(sjpup,sjpuppref+"Phi"),
            cms.InputTag(sjpup,sjpuppref+"PartonFlavour"),
            cms.InputTag(sjpup,sjpuppref+"CSVv2"),
            cms.InputTag(sjpup,sjpuppref+"DeepCSV"),
            ),                         
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        toSave = cms.vstring(sjpuppref+"E",sjpuppref+"Pt",sjpuppref+"Eta",sjpuppref+"Phi", "allExtra"),
        )
    


addAK8CHS = True
addAK8PUPPI = True
addAK8CHS = False
addAK8PUPPI = False
if addAK8CHS:
    DMTreesDumper.physicsObjects.append(jet8CHS)
    DMTreesDumper.physicsObjects.append(subjetCHS)
if addAK8PUPPI:
    DMTreesDumper.physicsObjects.append(jet8Puppi)
    DMTreesDumper.physicsObjects.append(subjetPuppi)
