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
cmsDriver.py step2 --scenario pp --conditions auto:run3_data_prompt -s REPACK:DigiToApproxClusterRaw --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --era Run3_pp_on_PbPb_approxSiStripClusters --eventcontent REPACKRAW -n 100 --customise_commands process.rawPrimeDataRepacker.src='rawDataRepacker' --repacked --process ReHLT --filein file:/gpfs/ddn/cms/user/muraleed/PREana/rawprime/Muon_outputPhysicsRawPrimeFull.root
```
rereco
```
cmsDriver.py step3 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3_pp_on_PbPb_approxSiStripClusters --filein file:/gpfs/ddn/cms/user/muraleed/PREana/outputFiles/step2_REPACK.root
```
Finally the RAW dataset also re-run the reco::Tracks
```
cmsDriver.py step5 --conditions auto:run3_data_prompt -s RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --data --process reRECO --scenario pp -n 100 --repacked --era Run3_pp_on_PbPb --nThreads 254 --filein /store/data/Run2024F/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/382/216/00000/aadd1ab9-4eb8-4fb2-ac62-bdd1bebe882e.root
```
