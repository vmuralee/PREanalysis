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


## FED studies
To study the event size for each detector FED contributions are follow.

```
cmsRun python/Trackanalyzer_cfg.py 
```
