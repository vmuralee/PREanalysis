# PREanalysis

The tool used for the Next Generation Trigger (NGT) working pack 3 of Task 3.3 **Reduction of the RAW data size for HLT**.  This tool used for the following studies,
    - The RAW event size in Run-3 pp-collisions 
    - The Strip size compression and analysis
    - The Track clusters studies

The following section will go through each topic, how the packages are 

## instaliing CMSSW setup

```
cmsrel CMSSW_14_1_0_pre6
cd CMSSW_14_1_0_pre6/src
cmsenv
git cms-init
git pull git@github.com:vmuralee/cmssw.git myDev
git clone https://github.com/vmuralee/PREanalysis.git
scram b -j 32

```

## The Track clusters studies

```
cd $CMSSW_BASE/src/TrackAna
```
To perform the MC 
