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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/home-v/vmuralee/PREanalysis/raw/Muon_outputPhysicsHIPhysicsRawRECOPrime0.root' # File for Tracks from SiStripCluster
        #'file:/eos/home-v/vmuralee/PREanalysis/raw/step2_dump_Rawprime.root'  # File for Tracks from ApproxSiStripCluster
    )
)

# Other statements
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_HLT_v3', '')

process.trackAnalyzer = cms.EDAnalyzer('TrackAnalyzer',
                                       tracks = cms.InputTag('generalTracks'),
                                       #beamSpot = process.offlineBeamSpot#cms.InputTag('offlineBeamSpot')
)
process.outputTracks = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('ZLIB'),
    compressionLevel = cms.untracked.int32(7),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('file:/eos/home-v/vmuralee/PREanalysis/outputFiles/flatTracktest.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',   
    'keep *_*SiStripClusters2ApproxClusters*_*_*',
    'keep *_*SiStripClusterizerForRawPrime*_*_*'   
    )
)

process.TFileService = cms.Service("TFileService", fileName=cms.string('/eos/home-v/vmuralee/PREanalysis/outputFiles/flatTracktest.root'))

process.p = cms.Path(process.trackAnalyzer)
