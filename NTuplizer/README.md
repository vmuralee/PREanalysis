# PREanalysis
## instaliing CMSSW

```
cmsrel CMSSW_14_0_11
cd CMSSW_14_0_11/src
cmsenv

```

## Re-run HLT menu for ApproxSiStripCluster
creating config file for Re-run hlt menu.
```
hltGetConfiguration /users/vmuralee/PREmenu/V9 \
--globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
--data --unprescale --max-events 100 --eras Run3 \
--input /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root \
> prehlt.py
```

## Re-run RECO:Tracks
Follow the cmsDriver commands
```
cmsDriver.py  step3 --scenario pp --conditions 140X_dataRun3_Prompt_v3 -s REPACK:DigiToApproxClusterRaw --datatier RECO --eventcontent FEVTDEBUGHLT --era Run2_2018_pp_on_AA -n 10 --procModifiers approxSiStripClusters --customise_commands process.hltSiStripRawToDigi.ProductLabel='rawDataCollector';process.hltScalersRawToDigi.scalersInputTag='rawDataCollector' --process REHLT --filein /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root
```
rereco
```
cmsDriver.py step4 -s RAW2DIGI,L1Reco,RECO --conditions auto:phase1_2022_realistic_hi -n 10 --datatier RECO --eventcontent RECO --era Run3_pp_on_PbPb_approxSiStripClusters --filein file:step3_REPACK.root --hltProcess REHLT --process reRECO
```
change the output commands inside the config file

```
outputCommands = cms.untracked.vstring( 'drop *',
                                            'keep *_*siStripClusters*_*_*',
                                            'keep *_*generalTracks*_*_*',
      'keep *_hltSiStripClusters2ApproxClusters_*_*',
      'keep DetIds_hltSiStripRawToDigi_*_*',
      'keep FEDRawDataCollection_rawPrimeDataRepacker_*_*',
      'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep triggerTriggerEvent_*_*_*' ),
```
