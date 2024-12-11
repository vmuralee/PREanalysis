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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_14_0_18/RelValTTbar_14TeV/GEN-SIM-RECO/140X_mcRun3_2022_realistic_v12_STD_noPU_2022_reMC-v1/2590000/493f8853-13f3-483f-87c0-0d8c3d87fb38.root'
        #'/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/140X_mcRun3_2023_realistic_v3_STD_2023_noPU-v1/2580000/6d072f9b-382f-423f-91d5-10d8fd0c8a49.root'
        #'/store/mc/CMSSW_14_1_0_pre1/RelValSingleMuPt10/GEN-SIM-RECO/140X_mcRun3_2024_realistic_v4_STD_2024_noPU-v2/2550000/642e14a0-2cf6-46ec-af58-5af55ecb74a4.root'
        #'/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun3_2023_realistic_v3_RecoOnly_2023_PU-v1/2580000/019e5150-b584-4d8e-803d-58e731884eda.root'
        '/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun3_2023_realistic_v3_RecoOnly_2023_noPU-v1/2580000/2ccff1ec-cc92-45ec-892d-d553a1efd77e.root'
        #'file:../../NTuplizer/step_raw_RAW2DIGI_L1Reco_RECO.root'
        
    ),
    secondaryFileNames = cms.untracked.vstring(
        # '/store/mc/CMSSW_14_1_0_pre1/RelValSingleMuPt10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v4_STD_2024_noPU-v2/2550000/0d3c8d12-2527-44a5-bfbe-4ca5fcbdbe41.root',
        # '/store/mc/CMSSW_14_1_0_pre1/RelValSingleMuPt10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v4_STD_2024_noPU-v2/2550000/9f2b7613-9d58-4fce-86d5-ee775db12e00.root'
        # '/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/140X_mcRun3_2023_realistic_v3_STD_2023_noPU-v1/2580000/45fadc6a-9aae-4068-bf5c-ca297c2dc5cb.root',
        # '/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/140X_mcRun3_2023_realistic_v3_STD_2023_noPU-v1/2580000/7cd0ea6f-79d8-40c2-868c-41bab5fec396.root',
        # '/store/relval/CMSSW_14_0_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/140X_mcRun3_2023_realistic_v3_STD_2023_noPU-v1/2580000/b3026954-ec19-412f-baa2-7c288dfc022f.root'
   ),
)

# Other statements
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'140X_mcRun3_2022_realistic_v12' , '')#'140X_dataRun3_Prompt_v3'

process.trackAnalyzer = cms.EDAnalyzer('TrackAnalyzer',
                                       IsMC = cms.untracked.bool(False),
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

process.TFileService = cms.Service("TFileService", fileName=cms.string('/scratchssd/muraleed/PREanalysis/flat_ntuples/Flat_RawRelValMinBiasTrack.root'))#'Flat_RawPrimeMCMinBiasNoPUTrackOT.root'))

process.p = cms.Path(process.trackAnalyzer)
