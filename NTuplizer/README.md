# PREanalysis
## instaliing CMSSW

```
cmsrel CMSSW_14_0_11
cd CMSSW_14_0_11/src
cmsenv

git cms-merge-topic vmuralee:myDev
scram b -j 32

```
The approximated SiStripCluster has the compressed barycenter scaled to 8 bit value maximum. 
## Re-run HLT menu for ApproxSiStripCluster
creating config file for Re-run hlt menu.
```
hltGetConfiguration /users/vmuralee/PREmenu/V9 \
--globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
--data --unprescale --max-events 100 --eras Run3 \
--input /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root \
> prehlt.py
```
The `--globaltag 140X_dataRun3_Prompt_v3` is used for offline cluster and track matching. To store only the important output modules,

```
outputCommands = cms.untracked.vstring('drop *',
      'keep *_*siStripClusters*_*_*',
      'keep *_*generalTracks*_*_*',
      'keep *_hltSiStripClusters2ApproxClusters_*_*',
      'keep DetIds_hltSiStripRawToDigi_*_HLTX',
      'keep FEDRawDataCollection_raw*_*_HLTX',
      'keep FEDRawDataCollection_hltSiStripDigiToZSRaw_*_HLTX',
      'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_HLTX',
      'keep edmTriggerResults_*_*_HLTX',
      'keep triggerTriggerEvent_*_*_HLTX')

```
The Reruning HLT by `cmsRun prehlt.py`.
## Re-run RECO:Tracks 
Follow the cmsDriver commands for Raw' dataset
```
cmsDriver.py step2 --scenario pp --conditions auto:run3_data_prompt -s REPACK:DigiToApproxClusterRaw --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --era Run3_pp_on_PbPb_approxSiStripClusters --eventcontent REPACKRAW -n 100 --customise_commands "process.rawPrimeDataRepacker.src='rawDataRepacker'" --repacked --process ReHLT --filein file:/gpfs/ddn/cms/user/muraleed/PREana/rawprime/Muon_outputPhysicsRawPrimeUint16check_t.root --no_exec 
```
change the era condition in `step2_REPACK.py` ,as 

```
from Configuration.Eras.Era_Run2024_pp_on_PbPb_approxSiStripCluster import Run3_pp_on_PbPb_approxSiStripClusters_2024
process = cms.Process('ReHLT',Run3_pp_on_PbPb_approxSiStripClusters_2024)
```
and reRun HLT.
rereco
```
cmsDriver.py step3 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3_pp_on_PbPb_approxSiStripClusters --filein file:/gpfs/ddn/cms/user/muraleed/PREana/outputFiles/step2_REPACK.root
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
change the era condition in `step3_RAW2DIGI_L1Reco_RECO.py` ,as 

```
from Configuration.Eras.Era_Run2024_pp_on_PbPb_approxSiStripCluster import Run3_pp_on_PbPb_approxSiStripClusters_2024
process = cms.Process('ReHLT',Run3_pp_on_PbPb_approxSiStripClusters_2024)
```
and reRun RECO.
Finally the RAW dataset also re-run the reco::Tracks
```
cmsDriver.py step5 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3 --nThreads 254 --filein /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root
```
