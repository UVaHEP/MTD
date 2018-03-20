# DRS-UVa

This is a starting point to do reconstructions of TestBeam data for DRS readout of LYSO+SiPM with Ratios
It is based on

https://github.com/CaltechPrecisionTiming/FNALTestbeam_052017/tree/master/Analysis/BTL


Clean install and running
```
git init
git clone https://github.com/UVaHEP/MTD.git
cd DRS-UVa
make clean
make
./BTL_Analysis --inputRootFile=../data/Run1816.root --mergeChannels=<true/false>


--inputRootFile --> Specify your own data file
--mergeChannels --> boolean to tell code whether to merge channels. Leaving this option off defaults to False
```
Edit
```
BTL_Analysis.cc 
```
to choose a DRS channel to investigate or calibrate
