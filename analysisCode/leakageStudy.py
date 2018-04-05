# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 5, 2018
#Purpose: Script to combine executables for 1) converting testbeam .dat files to .root, and 2) analyzing root file and producing timing-relevant information by channel

import os,sys, argparse
from ROOT import gROOT, TH1, TFile, TTree, TChain, TCanvas

c1 = TCanvas('c1','c1',200,10,750,940)

f0 = TFile('../../testbeam_FNAL_03-2018/all5exposure/DataCMSVMETiming_Run916.root', 'READ')
t0 = f0.pulse

print t0.GetEntries()

for event in t0:
    print event.amp[0]
