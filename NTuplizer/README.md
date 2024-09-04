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
Follow the cmsDriver commands for Raw' dataset
```
cmsDriver.py step2 --scenario pp --conditions auto:run3_data_prompt -s REPACK:DigiToApproxClusterRaw --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --era Run3 --eventcontent REPACKRAW -n 100 --customise_commands "process.rawPrimeDataRepacker.src='rawDataRepacker'" --repacked --process ReHLT --filein file:/gpfs/ddn/cms/user/muraleed/PREana/rawprime/Muon_outputPhysicsRawPrimeUint16check_t.root
```
The jets get empty in `--era Run3_pp_on_PbPb_approxSiStripClusters` but not in the `--era Run3`
rereco
```
cmsDriver.py step3 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3 --filein file:/gpfs/ddn/cms/user/muraleed/PREana/outputFiles/step2_REPACK.root
```
And change the outputCommands for the step3 configurateion by,
```
'drop *',
'keep *_*siStripClusters*_*_*',
'keep *_*generalTracks*_*_*',
'keep *_hltSiStripClusters2ApproxClusters_*_*',
'keep DetIds_hltSiStripRawToDigi_*_ReHLT',
'keep FEDRawDataCollection_raw*_*_ReHLT',
'keep FEDRawDataCollection_hltSiStripDigiToZSRaw_*_ReHLT',
'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_ReHLT',
'keep edmTriggerResults_*_*_ReHLT',
'keep triggerTriggerEvent_*_*_ReHLT'
```

Finally the RAW dataset also re-run the reco::Tracks
```
cmsDriver.py step5 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3 --nThreads 254 --filein /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root
```
