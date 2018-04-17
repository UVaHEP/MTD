# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 6, 2018
#Purpose: Script to combine executables for 1) converting testbeam .dat files to .root, and 2) analyzing root file and producing timing-relevant information by channel

from barClass import barClass
from ROOT import TFile, TTree

f0 = TFile('/home/btannenw/Desktop/MTD/testbeam_FNAL_03-2018/all5exposure/DataCMSVMETiming_5barExposure.root', 'READ') # all 5
#f0 = TFile('/home/btannenw/Desktop/MTD/testbeam_FNAL_03-2018/lowBias/bottombars_66V.root', 'READ') # low bias (66 V) bars 1 and 2
t0 = f0.pulse

barClass(t0, 'all5exposure', '04-16-18_plots')
