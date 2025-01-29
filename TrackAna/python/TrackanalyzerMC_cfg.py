import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            #'file:/scratchssd/muraleed/PREanalysis/rawtuples/RelValSingleMuonSIMRECO_140X_NoPU.root'
            '/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/140X_mcRun3_2024_realistic_v3_STD_2024_noPU-v1/2580000/46464e62-58a0-4e19-90cd-9331f366811f.root'
            #'/store/relval/CMSSW_14_0_0/RelValSingleMuPt10/GEN-SIM-RECO/140X_mcRun3_2024_realistic_v3_STD_2024_noPU-v1/2580000/994d0b6f-3a15-458c-aabf-ebad2bf0204e.root'
            #'/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/140X_mcRun3_2024_realistic_v3_RecoOnly_2024_PU-v1/2580000/105af884-df6a-49bc-a543-151d1c57a239.root',  #PU pi
            # '/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/140X_mcRun3_2024_realistic_v3_RecoOnly_2024_PU-v1/2580000/45cf45fe-ff66-434a-af9e-90d46e0204e6.root'   #PU pi
    ),
    secondaryFileNames = cms.untracked.vstring(
        #'file:/scratchssd/muraleed/PREanalysis/rawtuples/RelValSingleMuonSIMDIGI_140X_NoPU.root'
        #'/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM/140X_mcRun3_2024_realistic_v3_STD_2024_PU-v1/2580000/97abb881-e995-45c9-b0b2-dce45231dc15.root' #PU pi
        '/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v3_STD_2024_noPU-v1/2580000/64ad038b-b936-430c-9f46-698722871c01.root',
        '/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v3_STD_2024_noPU-v1/2580000/b6f9bcb0-d4de-4d50-88a9-298d01f019f4.root'
        # '/store/relval/CMSSW_14_0_0/RelValSingleMuPt10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v3_STD_2024_noPU-v1/2580000/7e9ab006-5920-4a46-bb66-bb55bc3ea443.root',
        # '/store/relval/CMSSW_14_0_0/RelValSingleMuPt10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v3_STD_2024_noPU-v1/2580000/d1de063a-ae29-477c-9acb-1e120ff8a8f5.root'

   ),
)

# Other statements
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'140X_mcRun3_2024_realistic_v3' , '')#'140X_dataRun3_Prompt_v3'

process.trackAnalyzer = cms.EDAnalyzer('TrackAnalyzer',
                                       IsMC = cms.untracked.bool(True),
                                       tracks = cms.InputTag('generalTracks','','RECO'),
                                       jets = cms.InputTag('ak4PFJets','','RECO'),
                                       g4SimHitsLowTofTag = cms.InputTag("g4SimHits","TrackerHitsTOBLowTof","SIM"),
                                       g4SimHitsHighTofTag = cms.InputTag("g4SimHits","TrackerHitsTOBHighTof","SIM"),
                                       g4SimHitsTrackTag = cms.InputTag("g4SimHits","","SIM")
                                   )
process.outputTracks = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('ZLIB'),
    compressionLevel = cms.untracked.int32(7),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('file:/scratchssd/muraleed/PREanalysis/flat_ntuples/Flat_RawPrimeIntTrackTest.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',   
    'keep *_*generalTracks*_*_*',
    'keep *_*SiStripClusters2ApproxClusters*_*_*',
    'keep *_*SiStripClusterizerForRawPrime*_*_*'   
    )
)
process.options.numberOfThreads = 8
process.options.numberOfStreams = 0

process.TFileService = cms.Service("TFileService", fileName=cms.string('/scratchssd/muraleed/PREanalysis/flat_ntuples/Flat_RawRelValPionTrackOTSIM.root'))#/scratchssd/muraleed/PREanalysis/flat_ntuples/Flat_RawRelValMuonTrackOT.root

process.p = cms.Path(process.trackAnalyzer)
