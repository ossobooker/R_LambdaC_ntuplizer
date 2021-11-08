#!/bin/bash

for i in BParkingNANO_v01_2021Nov07/crab_* 
do 
    crab status -d $i
    crab resubmit -d $i
done            
