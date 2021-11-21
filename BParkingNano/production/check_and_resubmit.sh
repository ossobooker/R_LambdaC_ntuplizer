#!/bin/bash

for i in BParkingNANO_v02_2021Nov10/crab_* 
do 
    crab status -d $i
    crab resubmit -d $i
    # crab kill -d $i
done            


rm -rf /eos/cms/store/user/ftorresd/buffer/BParkingNANO_v01_2021Nov07/ParkingBPH1
rm -rf /eos/cms/store/user/ftorresd/buffer/BParkingNANO_v01_2021Nov07/ParkingBPH2
rm -rf /eos/cms/store/user/ftorresd/buffer/BParkingNANO_v01_2021Nov07/ParkingBPH3
rm -rf /eos/cms/store/user/ftorresd/buffer/BParkingNANO_v01_2021Nov07/ParkingBPH5
