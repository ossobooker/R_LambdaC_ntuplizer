# nanoAOD producer customized for BParking analysis 

Forked from: https://github.com/CMSBParking/BParkingNANO

The focus is on R($\Lambda_{c}^{(*)}$) analyses.

## More info

https://codimd.web.cern.ch/mKQBQ3ErQ9uKMr4TdS1udA?both

## Setup Env

```shell
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_27
cd CMSSW_10_2_27/src
cmsenv
source /cvmfs/cms.cern.ch/common/crab-setup.sh
git cms-init
pip install --user yapf
```
<!-- ## Add low-pT energy ID and regression

The ID model is `2020Sept15` (depth=15, ntrees=1000).

```shell
git cms-merge-topic -u CMSBParking:from-CMSSW_10_2_15_2020Sept15_v1
git clone --single-branch --branch from-CMSSW_10_2_15_2020Sept15 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

To run on CRAB, the following three lines __must__ be executed:

```shell
git cms-addpkg RecoEgamma/ElectronIdentification
mkdir -p $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data/LowPtElectrons
cp $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data/LowPtElectrons
``` -->

<!-- ## Add support for GBRForest to parse ROOT files

```shell
git cms-merge-topic -u CMSBParking:convertXMLToGBRForestROOT
``` -->

<!-- ## Add the modification needed to use post-fit quantities for electrons  

```shell
git cms-merge-topic -u CMSBParking:GsfTransientTracks # unsafe checkout (no checkdeps), but suggested here
``` -->

## Add the modification needed to use the KinematicParticleVertexFitter  

```shell
git cms-merge-topic -u CMSBParking:fixKinParticleVtxFitter # unsafe checkout (no checkdeps), but suggested here
```

## Add the BParkingNano package and build everything

```shell
git clone git@github.com:CMSBParking/BParkingNANO.git ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```

## To run on a test file

```shell
cd PhysicsTools/BParkingNano/test/
cmsenv 
cmsRun run_nano_cfg.py
```
## Quick Start

```
source setup_env.sh
```
