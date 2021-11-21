#!/bin/bash

for i in /eos/cms/store/user/ftorresd/buffer/BParkingNANO_v01_2021Nov07/BsToPhiJpsi_ToKKMuMu_probefilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BsToPhiJpsi_ext1-v1/211107_173851/0000/BParkNANO_mc_v01_2021Nov07_*.root
do 
    echo "--->>>> $i "
    root -l -b -q $i
done            




