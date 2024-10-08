# PREanalysis
## instaliing CMSSW

```
cmsrel CMSSW_14_1_0_pre6
cd CMSSW_14_1_0_pre6/src/
cmsenv
git cms-init
git clone https://github.com/vmuralee/PREanalysis.git -b FED_study
scram b -j 32

```
## Re-run HLT menu for FED study
Go to `cd PREanalysis/NTuplizer` directory, The config file used for the study of each FED contribution is `prehlt_FED.py`. The file **FEDlist.py** contains the pin info of each detector. To run the config file,
```
cmsRun prehlt_FED.py <detector>
```
The list of detector key words are,
- Pixel   : Pixel contribution
- Strips  : Strip contribution
- ECAL    : ECAL contribution
- HCAL    : HCAL contribution
- Others  : Other minor contribution

By running the cmsRun command will produce the output file, which only has the contribution of given detector FED. The average compressed size of  **FEDRawDataCollection_rawPrimeDataRepacker__HLTX** can find by `edmEventSize -v <outputFile>`. 
  

