# PREanalysis
## instaliing CMSSW

```
cmsrel CMSSW_14_1_0_pre6
cd CMSSW_14_1_0_pre6/src/
cmsenv
git cms-init
git pull git@github.com:vmuralee/cmssw.git myDev
git cms-addpkg Configuration/Eras DataFormats/SiStripCluster
git clone https://github.com/vmuralee/PREanalysis.git -b instructions
scram b -j 32

```
The approximated SiStripCluster has the compressed barycenter scaled to 8 bit value maximum. 
## Re-run HLT menu for ApproxSiStripCluster
Go to `cd PREanalysis/NTuplizer` directory and creating the config file for Re-run hlt menu.

```
hltGetConfiguration /users/vmuralee/PREmenu/V9 \
--globaltag 140X_dataRun3_Prompt_v3\
--data --unprescale --max-events 100 --eras Run3 \
--input /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root \
> prehlt.py
```
change the output commands as,

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
Or simply Rerun existing `prehlt.py` file.
## Re-run the reco step in RAW' dataset
Follow the cmsDriver commands for Raw' dataset
```
cmsDriver.py step2 --scenario pp --conditions auto:run3_data_prompt -s REPACK:DigiToApproxClusterRaw --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --era Run3_pp_on_PbPb_approxSiStripClusters --eventcontent REPACKRAW -n 100 --customise_commands "process.rawPrimeDataRepacker.src='rawDataRepacker'" --repacked --process ReHLT --filein file:Muon_outputPhysicsHIPhysicsRawPrime0.root --no_exec 
```
command out the following lines in `step2_REPACK.py`,
```
from Configuration.Eras.Era_Run3_pp_on_PbPb_approxSiStripClusters_cff import Run3_pp_on_PbPb_approxSiStripClusters
process = cms.Process('ReHLT',Run3_pp_on_PbPb_approxSiStripClusters)
```
and include the era condition, as 

```
from Configuration.Eras.Era_Run2024_pp_on_PbPb_approxSiStripCluster import Run3_pp_on_PbPb_approxSiStripClusters_2024
process = cms.Process('ReHLT',Run3_pp_on_PbPb_approxSiStripClusters_2024)
```
and proceed with `cmsRun step2_REPACK.py`. The output file is used for reRECO process. 
```
cmsDriver.py step3 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3_pp_on_PbPb_approxSiStripClusters --filein file:step2_REPACK.root --no_exec 
```
And change the outputCommands for the step3 configurateion `step3_RAW2DIGI_L1Reco_RECO.py`  by,
```
'drop *',
'keep *_ak4PFJets_*_*',
'keep *_*pfMet*_*_*',
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
command out the following lines in `step3_RAW2DIGI_L1Reco_RECO.py`,
```
from Configuration.Eras.Era_Run3_pp_on_PbPb_approxSiStripClusters_cff import Run3_pp_on_PbPb_approxSiStripClusters
process = cms.Process('reRECO',Run3_pp_on_PbPb_approxSiStripClusters)
```
and include the era condition, as 

```
from Configuration.Eras.Era_Run2024_pp_on_PbPb_approxSiStripCluster import Run3_pp_on_PbPb_approxSiStripClusters_2024
process = cms.Process('reRECO',Run3_pp_on_PbPb_approxSiStripClusters_2024)
```
After, running `step3_RAW2DIGI_L1Reco_RECO.py` file will produced reconstructed tracks,jets and met. 

## Re-run the reco step in RAW dataset.

Finally the RAW dataset also re-run with the condition `auto:run3_data_prompt`. 
```
cmsDriver.py step5 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3_2024 --nThreads 254 --filein /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root
```
