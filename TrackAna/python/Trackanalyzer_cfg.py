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
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(11))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:../../NTuplizer/step_raw_RAW2DIGI_L1Reco_RECO.root'
        'file:/scratchssd/muraleed/PREanalysis/rereco_ntuples/Muon_RecoRawPrime0Int16bit.root'
        #'file:/scratchssd/muraleed/PREanalysis/rereco_ntuples/Muon_RecoRawPrime0Int8bit.root' # File for Tracks from SiStripCluster
        #'file:/eos/home-v/vmuralee/PREanalysis/raw/step2_dump_Rawprime.root'  # File for Tracks from ApproxSiStripCluster
    )
)

# Other statements
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_Prompt_v3', '')

process.trackAnalyzer = cms.EDAnalyzer('TrackAnalyzer',
                                       tracks = cms.InputTag('generalTracks','','reRECO'),
                                       jets = cms.InputTag('ak4PFJets','','reRECO')
                                       
)
process.outputTracks = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('ZLIB'),
    compressionLevel = cms.untracked.int32(7),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('file:/scratchssd/muraleed/PREanalysis/flat_ntuples/Flat_RawPrimeIntTrackTest.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',   
    'keep *_*SiStripClusters2ApproxClusters*_*_*',
    'keep *_*SiStripClusterizerForRawPrime*_*_*'   
    )
)

process.TFileService = cms.Service("TFileService", fileName=cms.string('/scratchssd/muraleed/PREanalysis/flat_ntuples/Flat_TrackTest.root'))

process.p = cms.Path(process.trackAnalyzer)
