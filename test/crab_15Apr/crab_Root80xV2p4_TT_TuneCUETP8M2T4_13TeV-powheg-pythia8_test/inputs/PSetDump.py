import FWCore.ParameterSet.Config as cms

process = cms.Process("ttDManalysisTrees")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:B2GEDMNtuple_1.root')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.DMTreesDumper = cms.EDAnalyzer("DMAnalysisTreeMaker",
    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"),
    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
    EraLabel = cms.untracked.string('BCD'),
    HBHEIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),
    ak8jetvSubjetIndex0 = cms.InputTag("jetsAK8","jetAK8vSubjetIndex0"),
    ak8jetvSubjetIndex1 = cms.InputTag("jetsAK8","jetAK8vSubjetIndex1"),
    applyRes = cms.untracked.bool(False),
    boostedTopsLabel = cms.string('jetsAK8CHS'),
    boostedTopsSubjetsLabel = cms.string('subjetsAK8CHS'),
    changeJECs = cms.untracked.bool(True),
    channelInfo = cms.PSet(
        HadronPFHT800Triggers = cms.vstring('HLT_PFHT800', 
            'HLT_PFHT800_v0', 
            'HLT_PFHT800_v1', 
            'HLT_PFHT800_v2', 
            'HLT_PFHT800_v3', 
            'HLT_PFHT800_v4', 
            'HLT_PFHT800_v5', 
            'HLT_PFHT800_v6', 
            'HLT_PFHT800_v7', 
            'HLT_PFHT800_v8', 
            'HLT_PFHT800_v9', 
            'HLT_PFHT800_v10', 
            'HLT_PFHT800_v11', 
            'HLT_PFHT800_v12', 
            'HLT_PFHT800_v13', 
            'HLT_PFHT800_v14'),
        HadronPFHT900Triggers = cms.vstring('HLT_PFHT900', 
            'HLT_PFHT900_v0', 
            'HLT_PFHT900_v1', 
            'HLT_PFHT900_v2', 
            'HLT_PFHT900_v3', 
            'HLT_PFHT900_v4', 
            'HLT_PFHT900_v5', 
            'HLT_PFHT900_v6', 
            'HLT_PFHT900_v7', 
            'HLT_PFHT900_v8', 
            'HLT_PFHT900_v9', 
            'HLT_PFHT900_v10', 
            'HLT_PFHT900_v11', 
            'HLT_PFHT900_v12', 
            'HLT_PFHT900_v13', 
            'HLT_PFHT900_v14'),
        HadronPFJet450Triggers = cms.vstring('HLT_PFJet450', 
            'HLT_PFJet450_v0', 
            'HLT_PFJet450_v1', 
            'HLT_PFJet450_v2', 
            'HLT_PFJet450_v3', 
            'HLT_PFJet450_v4', 
            'HLT_PFJet450_v5', 
            'HLT_PFJet450_v6', 
            'HLT_PFJet450_v7', 
            'HLT_PFJet450_v8', 
            'HLT_PFJet450_v9', 
            'HLT_PFJet450_v10', 
            'HLT_PFJet450_v11', 
            'HLT_PFJet450_v12', 
            'HLT_PFJet450_v13', 
            'HLT_PFJet450_v14'),
        PhotonTriggers = cms.vstring(''),
        SingleElTriggers = cms.vstring('HLT_Ele27_eta2p1_WP75_Gsf', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v0', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v1', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v2', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v3', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v4', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v5', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v6', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v7', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v8', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v9', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v10', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v11', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v12', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v13', 
            'HLT_Ele27_eta2p1_WP75_Gsf_v14'),
        SingleMuTriggers = cms.vstring('HLT_Mu50', 
            'HLT_TkMu50', 
            'HLT_Mu50_v0', 
            'HLT_Mu50_v1', 
            'HLT_Mu50_v2', 
            'HLT_Mu50_v3', 
            'HLT_Mu50_v4', 
            'HLT_Mu50_v5', 
            'HLT_Mu50_v6', 
            'HLT_Mu50_v7', 
            'HLT_Mu50_v8', 
            'HLT_Mu50_v9', 
            'HLT_Mu50_v10', 
            'HLT_Mu50_v11', 
            'HLT_Mu50_v12', 
            'HLT_Mu50_v13', 
            'HLT_Mu50_v14', 
            'HLT_TkMu50_v0', 
            'HLT_TkMu50_v1', 
            'HLT_TkMu50_v2', 
            'HLT_TkMu50_v3', 
            'HLT_TkMu50_v4', 
            'HLT_TkMu50_v5', 
            'HLT_TkMu50_v6', 
            'HLT_TkMu50_v7', 
            'HLT_TkMu50_v8', 
            'HLT_TkMu50_v9', 
            'HLT_TkMu50_v10', 
            'HLT_TkMu50_v11', 
            'HLT_TkMu50_v12', 
            'HLT_TkMu50_v13', 
            'HLT_TkMu50_v14'),
        addLHAPDFWeights = cms.untracked.bool(False),
        channel = cms.string('ttDM'),
        crossSection = cms.double(1),
        doTopReweighting = cms.untracked.bool(True),
        doWReweighting = cms.untracked.bool(True),
        getPartonTop = cms.untracked.bool(True),
        getPartonW = cms.untracked.bool(True),
        hadronicTriggers = cms.vstring('HLT_PFMETNoMu120_PFMHTNoMu120_IDTight', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v0', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v1', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v2', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v3', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v4', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v5', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v6', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v7', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v8', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v9', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v10', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v11', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v12', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v13', 
            'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v14', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v0', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v1', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v2', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v3', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v4', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v5', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v6', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v7', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v8', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v9', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v10', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v12', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v13', 
            'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v14'),
        maxWeights = cms.untracked.int32(111),
        metFilters = cms.vstring('Flag_HBHENoiseFilter', 
            'Flag_HBHENoiseIsoFilter', 
            'Flag_EcalDeadCellTriggerPrimitiveFilter', 
            'Flag_goodVertices', 
            'Flag_globalTightHalo2016Filter'),
        originalEvents = cms.double(1),
        useLHE = cms.untracked.bool(True),
        useLHEWeights = cms.untracked.bool(True)
    ),
    cutOnTriggers = cms.untracked.bool(False),
    doPreselection = cms.untracked.bool(False),
    doResolvedTopHad = cms.untracked.bool(False),
    doResolvedTopSemiLep = cms.untracked.bool(False),
    eleLabel = cms.string('electrons'),
    eventLabel = cms.string('eventShapePFVars'),
    eventNumber = cms.InputTag("eventInfo","evtInfoEventNumber"),
    genLabel = cms.string('genPart'),
    genprod = cms.InputTag("generator"),
    getPartonTop = cms.untracked.bool(True),
    getPartonW = cms.untracked.bool(True),
    isData = cms.untracked.bool(False),
    jetKeysAK4CHS = cms.InputTag("jetKeysAK4CHS"),
    jetScanCuts = cms.vdouble(30),
    jetsLabel = cms.string('jetsAK4CHS'),
    lhes = cms.InputTag("externalLHEProducer"),
    lumiBlock = cms.InputTag("eventInfo","evtInfoLumiBlock"),
    metBits = cms.InputTag("METUserData","triggerBitTree"),
    metLabel = cms.string('metFull'),
    metNames = cms.InputTag("METUserData","triggerNameTree"),
    muLabel = cms.string('muons'),
    muonKeys = cms.InputTag("muonKeys"),
    photonLabel = cms.string('photons'),
    physicsObjects = cms.VPSet(cms.PSet(
        categories = cms.vstring(),
        label = cms.string('metFull'),
        maxInstances = cms.untracked.int32(1),
        prefix = cms.string('metFull'),
        saveBaseVariables = cms.untracked.bool(True),
        singleD = cms.VInputTag(),
        singleF = cms.VInputTag(),
        singleI = cms.VInputTag(),
        toSave = cms.vstring(),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(cms.InputTag("metFull","metFulluncorPt"), cms.InputTag("metFull","metFulluncorPhi"), cms.InputTag("metFull","metFullPt"), cms.InputTag("metFull","metFullPhi"), cms.InputTag("metFull","metFullPx"), 
            cms.InputTag("metFull","metFullPy")),
        variablesI = cms.VInputTag()
    ), 
        cms.PSet(
            categories = cms.vstring('Tight', 
                'Loose', 
                'Medium'),
            label = cms.string('muons'),
            maxInstances = cms.untracked.int32(6),
            prefix = cms.string('mu'),
            saveBaseVariables = cms.untracked.bool(False),
            saveNoCat = cms.untracked.bool(True),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('muE', 
                'muPt', 
                'muEta', 
                'muPhi', 
                'muIso04', 
                'muCharge', 
                'muIsMediumMuon', 
                'muIsTightMuon', 
                'muIsLooseMuon', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("muons","muE"), cms.InputTag("muons","muPt"), cms.InputTag("muons","muEta"), cms.InputTag("muons","muPhi"), cms.InputTag("muons","muCharge"), 
                cms.InputTag("muons","muIsLooseMuon"), cms.InputTag("muons","muIsSoftMuon"), cms.InputTag("muons","muIsMediumMuon"), cms.InputTag("muons","muIsTightMuon"), cms.InputTag("muons","muDB"), 
                cms.InputTag("muons","muDBerr"), cms.InputTag("muons","muDz"), cms.InputTag("muons","muDxy"), cms.InputTag("muons","muGenMuonCharge"), cms.InputTag("muons","muGenMuonEta"), 
                cms.InputTag("muons","muGenMuonPt"), cms.InputTag("muons","muGenMuonE"), cms.InputTag("muons","muGenMuonPhi"), cms.InputTag("muons","muGlbTrkNormChi2"), cms.InputTag("muons","muInTrkNormChi2"), 
                cms.InputTag("muons","muIsGlobalMuon"), cms.InputTag("muons","muIsPFMuon"), cms.InputTag("muons","muIsTrackerMuon"), cms.InputTag("muons","muIso04"), cms.InputTag("muons","muNumberMatchedStations"), 
                cms.InputTag("muons","muNumberOfPixelLayers"), cms.InputTag("muons","muNumberOfValidTrackerHits"), cms.InputTag("muons","muNumberTrackerLayers"), cms.InputTag("muons","muNumberValidMuonHits"), cms.InputTag("muons","muNumberValidPixelHits"), 
                cms.InputTag("muons","muSumChargedHadronPt"), cms.InputTag("muons","muSumNeutralHadronPt"), cms.InputTag("muons","muSumPUPt"), cms.InputTag("muons","muSumPhotonPt")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('photons'),
            maxInstances = cms.untracked.int32(6),
            prefix = cms.string('pho'),
            saveBaseVariables = cms.untracked.bool(True),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('phoPt', 
                'phoHoverE', 
                'phoSigmaIEtaIEta', 
                'phoNeutralHadronIso', 
                'phoNeutralHadronIsoEAcorrected', 
                'phoChargedHadronIso', 
                'phoChargedHadronIsoEAcorrected', 
                'phoPhotonIso', 
                'phoPhotonIsoEAcorrected', 
                'phoHasPixelSeed', 
                'phoPassLooseID', 
                'phoPassMediumID', 
                'phoPassTightID', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("photons","phoEta"), cms.InputTag("photons","phoPt"), cms.InputTag("photons","phoHoverE"), cms.InputTag("photons","phoSigmaIEtaIEta"), cms.InputTag("photons","phoNeutralHadronIso"), 
                cms.InputTag("photons","phoChargedHadronIso"), cms.InputTag("photons","phoPhotonIso"), cms.InputTag("photons","phoHasPixelSeed"), cms.InputTag("photons","phoPassLooseID"), cms.InputTag("photons","phoPassMediumID"), 
                cms.InputTag("photons","phoPassTightID")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring('Tight', 
                'Veto'),
            label = cms.string('electrons'),
            maxInstances = cms.untracked.int32(6),
            prefix = cms.string('el'),
            saveBaseVariables = cms.untracked.bool(False),
            saveNoCat = cms.untracked.bool(True),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('elE', 
                'elPt', 
                'elEta', 
                'elPhi', 
                'elIso03', 
                'elisTight', 
                'elisMedium', 
                'elisLoose', 
                'elisVeto', 
                'elscEta', 
                'allExtra', 
                'elvidVeto', 
                'elvidLoose', 
                'elvidMedium', 
                'elvidTight'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("electrons","elE"), cms.InputTag("electrons","elPt"), cms.InputTag("electrons","elEta"), cms.InputTag("electrons","elPhi"), cms.InputTag("electrons","elCharge"), 
                cms.InputTag("electrons","elDz"), cms.InputTag("electrons","elEta"), cms.InputTag("electrons","elSCEta"), cms.InputTag("electrons","elHoE"), cms.InputTag("electrons","elIso03"), 
                cms.InputTag("electrons","eldEtaIn"), cms.InputTag("electrons","eldPhiIn"), cms.InputTag("electrons","elmissHits"), cms.InputTag("electrons","elfull5x5siee"), cms.InputTag("electrons","elooEmooP"), 
                cms.InputTag("electrons","elhasMatchedConVeto"), cms.InputTag("electrons","elvidVeto"), cms.InputTag("electrons","elvidTight"), cms.InputTag("electrons","elvidMedium"), cms.InputTag("electrons","elvidVeto"), 
                cms.InputTag("electrons","elvidLoose"), cms.InputTag("electrons","elvidMedium"), cms.InputTag("electrons","elvidTight")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('eventShapePFVars'),
            maxInstances = cms.untracked.int32(1),
            prefix = cms.string(''),
            saveBaseVariables = cms.untracked.bool(True),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('eventShapePFJetVars'),
            maxInstances = cms.untracked.int32(1),
            prefix = cms.string(''),
            saveBaseVariables = cms.untracked.bool(True),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring('Tight'),
            label = cms.string('jetsAK4CHS'),
            maxInstances = cms.untracked.int32(20),
            prefix = cms.string('jetAK4CHS'),
            saveBaseVariables = cms.untracked.bool(False),
            saveNoCat = cms.untracked.bool(True),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('jetAK4CHSE', 
                'jetAK4CHSPt', 
                'jetAK4CHSEta', 
                'jetAK4CHSPhi', 
                'jetAK4CHSGenJetE', 
                'jetAK4CHSGenJetPhi', 
                'jetAK4CHSGenJetPt', 
                'jetAK4CHSGenJetEta', 
                'jetAK4CHSCSVv2', 
                'jetAK4CHSPartonFlavour', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("jetsAK4CHS","jetAK4CHSE"), cms.InputTag("jetsAK4CHS","jetAK4CHSPt"), cms.InputTag("jetsAK4CHS","jetAK4CHSEta"), cms.InputTag("jetsAK4CHS","jetAK4CHSPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSQGL"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSPartonFlavour"), cms.InputTag("jetsAK4CHS","jetAK4CHSPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSCSVv2"), cms.InputTag("jetsAK4CHS","jetAK4CHSCharge"), cms.InputTag("jetsAK4CHS","jetAK4CHSElectronEnergy"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetCharge"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetE"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetEta"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetPt"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonCharge"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonE"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonEta"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonPt"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSHadronFlavour"), cms.InputTag("jetsAK4CHS","jetAK4CHSchargedEmEnergyFrac"), cms.InputTag("jetsAK4CHS","jetAK4CHSchargedHadronEnergyFrac"), cms.InputTag("jetsAK4CHS","jetAK4CHSneutralEmEnergyFrac"), cms.InputTag("jetsAK4CHS","jetAK4CHSneutralHadronEnergyFrac"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSelectronMultiplicity"), cms.InputTag("jetsAK4CHS","jetAK4CHSmuonMultiplicity"), cms.InputTag("jetsAK4CHS","jetAK4CHSneutralMultiplicity"), cms.InputTag("jetsAK4CHS","jetAK4CHSphotonMultiplicity"), cms.InputTag("jetsAK4CHS","jetAK4CHSchargedMultiplicity"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSneutralEmEnergy"), cms.InputTag("jetsAK4CHS","jetAK4CHSjecFactor0"), cms.InputTag("jetsAK4CHS","jetAK4CHSjetArea")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('jetsAK8CHS'),
            maxInstances = cms.untracked.int32(10),
            prefix = cms.string('jetAK8CHS'),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('jetAK8CHSGenJetE', 
                'jetAK8CHSGenJetPt', 
                'jetAK8CHSGenJetEta', 
                'jetAK8CHSGenJetPhi', 
                'jetAK8CHSE', 
                'jetAK8CHSPt', 
                'jetAK8CHSEta', 
                'jetAK8CHSPhi', 
                'jetAK8CHSprunedMassCHS', 
                'jetAK8CHStau1CHS', 
                'jetAK8CHStau2CHS', 
                'jetAK8CHStau3CHS', 
                'jetAK8CHStrimmedMassCHS', 
                'jetAK8CHSPartonFlavour', 
                'jetAK8CHSHadronFlavour', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetE"), cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetPt"), cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetEta"), cms.InputTag("jetsAK8CHS","jetAK8CHSGenJetPhi"), cms.InputTag("jetsAK8CHS","jetAK8CHSE"), 
                cms.InputTag("jetsAK8CHS","jetAK8CHSPt"), cms.InputTag("jetsAK8CHS","jetAK8CHSEta"), cms.InputTag("jetsAK8CHS","jetAK8CHSPhi"), cms.InputTag("jetsAK8CHS","jetAK8CHSPartonFlavour"), cms.InputTag("jetsAK8CHS","jetAK8CHSHadronFlavour"), 
                cms.InputTag("jetsAK8CHS","jetAK8CHSvSubjetIndex0"), cms.InputTag("jetsAK8CHS","jetAK8CHSvSubjetIndex1"), cms.InputTag("jetsAK8CHS","jetAK8CHSprunedMassCHS"), cms.InputTag("jetsAK8CHS","jetAK8CHSsoftDropMassCHS"), cms.InputTag("jetsAK8CHS","jetAK8CHStau1CHS"), 
                cms.InputTag("jetsAK8CHS","jetAK8CHStau2CHS"), cms.InputTag("jetsAK8CHS","jetAK8CHStau3CHS"), cms.InputTag("jetsAK8CHS","jetAK8CHStrimmedMassCHS"), cms.InputTag("jetsAK8CHS","jetAK8CHSjecFactor0"), cms.InputTag("jetsAK8CHS","jetAK8CHSjetArea")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('subjetsAK8CHS'),
            maxInstances = cms.untracked.int32(10),
            prefix = cms.string('subjetAK8CHS'),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('subjetAK8CHSE', 
                'subjetAK8CHSPt', 
                'subjetAK8CHSEta', 
                'subjetAK8CHSPhi', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("subjetsAK8CHS","subjetAK8CHSE"), cms.InputTag("subjetsAK8CHS","subjetAK8CHSPt"), cms.InputTag("subjetsAK8CHS","subjetAK8CHSEta"), cms.InputTag("subjetsAK8CHS","subjetAK8CHSPhi"), cms.InputTag("subjetsAK8CHS","subjetAK8CHSPartonFlavour"), 
                cms.InputTag("subjetsAK8CHS","subjetAK8CHSCSVv2")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('jetsAK8CHS'),
            maxInstances = cms.untracked.int32(10),
            prefix = cms.string('jetAK8Puppi'),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('jetAK8PuppiGenJetE', 
                'jetAK8PuppiGenJetPt', 
                'jetAK8PuppiGenJetEta', 
                'jetAK8PuppiGenJetPhi', 
                'jetAK8PuppiE', 
                'jetAK8PuppiPt', 
                'jetAK8PuppiEta', 
                'jetAK8PuppiPhi', 
                'jetAK8PuppiprunedMassCHS', 
                'jetAK8Puppitau1CHS', 
                'jetAK8Puppitau2CHS', 
                'jetAK8Puppitau3CHS', 
                'jetAK8PuppiPartonFlavour', 
                'jetAK8PuppiHadronFlavour', 
                'jetAK8PuppitrimmedMassCHS', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetE"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetPt"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetEta"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiGenJetPhi"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiE"), 
                cms.InputTag("jetsAK8Puppi","jetAK8PuppiPt"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiEta"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiPhi"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiPartonFlavour"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiHadronFlavour"), 
                cms.InputTag("jetsAK8Puppi","jetAK8PuppivSubjetIndex0"), cms.InputTag("jetsAK8Puppi","jetAK8PuppivSubjetIndex1"), cms.InputTag("jetsAK8Puppi","jetAK8PuppiprunedMass"), cms.InputTag("jetsAK8Puppi","jetAK8PuppisoftDropMass"), cms.InputTag("jetsAK8Puppi","jetAK8Puppitau1"), 
                cms.InputTag("jetsAK8Puppi","jetAK8Puppitau2"), cms.InputTag("jetsAK8Puppi","jetAK8Puppitau3"), cms.InputTag("jetsAK8Puppi","jetAK8PuppitrimmedMass"), cms.InputTag("jetsAK8Puppi","jetAK8PuppijecFactor0"), cms.InputTag("jetsAK8Puppi","jetAK8PuppijetArea")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('subjetsAK8Puppi'),
            maxInstances = cms.untracked.int32(10),
            prefix = cms.string('subjetAK8Puppi'),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('subjetAK8PuppiE', 
                'subjetAK8PuppiPt', 
                'subjetAK8PuppiEta', 
                'subjetAK8PuppiPhi', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("subjetsAK8Puppi","subjetAK8PuppiE"), cms.InputTag("subjetsAK8Puppi","subjetAK8PuppiPt"), cms.InputTag("subjetsAK8Puppi","subjetAK8PuppiEta"), cms.InputTag("subjetsAK8Puppi","subjetAK8PuppiPhi"), cms.InputTag("subjetsAK8Puppi","subjetAK8PuppiPartonFlavour"), 
                cms.InputTag("subjetsAK8Puppi","subjetAK8PuppiCSVv2")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('genPart'),
            maxInstances = cms.untracked.int32(100),
            prefix = cms.string('genPart'),
            saveBaseVariables = cms.untracked.bool(False),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            toSave = cms.vstring('genpartE', 
                'genpartEta', 
                'genpartID', 
                'genpartMom0ID', 
                'genpartPhi', 
                'genpartPt', 
                'genpartStatus', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(),
            variablesI = cms.VInputTag()
        )),
    runNumber = cms.InputTag("eventInfo","evtInfoRunNumber"),
    systematics = cms.vstring('noSyst', 
        'jer__up', 
        'jer__down', 
        'jes__up', 
        'jes__down'),
    triggerBits = cms.InputTag("TriggerUserData","triggerBitTree"),
    triggerNames = cms.InputTag("TriggerUserData","triggerNameTree"),
    triggerPrescales = cms.InputTag("TriggerUserData","triggerPrescaleTree"),
    useTriggers = cms.untracked.bool(True),
    vertexChi2 = cms.InputTag("vertexInfo","chi"),
    vertexNdof = cms.InputTag("vertexInfo","ndof"),
    vertexRho = cms.InputTag("vertexInfo","rho"),
    vertexZ = cms.InputTag("vertexInfo","z")
)


process.analysisPath = cms.Path(process.DMTreesDumper)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary', 
        'HLTrigReport'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(100)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('treesTTJets.root')
)


process.CastorDbProducer = cms.ESProducer("CastorDbProducer")


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelQualityFromDbRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        ))
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(cms.PSet(
        Label = cms.untracked.string(''),
        NormalizationFactor = cms.untracked.double(1.0),
        Record = cms.string('SiStripApvGainRcd')
    ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiStripDetVOffRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('80X_mcRun2_asymptotic_2016_miniAODv2'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string('frontier://FrontierProd/'),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HERecalibration = cms.bool(False),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalibration = cms.bool(False),
    iLumi = cms.double(-1.0),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths')
)


process.prefer("es_hardcode")

