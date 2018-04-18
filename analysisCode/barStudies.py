# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 6, 2018
#Purpose: Script to combine executables for 1) converting testbeam .dat files to .root, and 2) analyzing root file and producing timing-relevant information by channel

from barClass import barClass
from ROOT import TFile, TTree
import sys, argparse


# *** 0. setup parser for command line
parser = argparse.ArgumentParser()
parser.add_argument("--test", help="flag for running over only 10k events", action='store_true')
parser.add_argument("--vetoOpt", help="veto decision logic option: none/singleAdj/doubleAdj/allAdj/all")
args = parser.parse_args()

if(args.vetoOpt is None):
    print "#### Seting --vetoOpt singlAdj by default ####"
    args.vetoOpt = 'singleAdj'
else:
    if( not(args.vetoOpt == "none" or args.vetoOpt == "singleAdj" or args.vetoOpt == "doubleAdj" or args.vetoOpt == "allAdj" or args.vetoOpt == "all") ):
        print "#### Please use none/singleAdj/doubleAdj/allAdj/all when setting --vetoOpt <option>. Supplied value ({0}) does not match ####\nEXITING".format(args.vetoOpt)
        quit()
    else:
        print '-- Setting vetoOpt = {0}'.format(args.vetoOpt)

if(args.test is not None):
    print "#### TEST MODE - running over 10k events only ####"



topDir = '04-18-18_plots_{0}'.format(args.vetoOpt)


# Ben Local
f0 = TFile('/home/btannenw/Desktop/MTD/testbeam_FNAL_03-2018/all5exposure/DataCMSVMETiming_5barExposure.root', 'READ') # all 5
f1 = TFile('/home/btannenw/Desktop/MTD/testbeam_FNAL_03-2018/lowBias/bottombars_66V.root', 'READ') # low bias (66 V) bars 1 and 2
#f2 = TFile('/home/btannenw/Desktop/MTD/testbeam_FNAL_03-2018/lowBias/topbars_66V.root', 'READ') # low bias (66 V) bars 3, 4, and 5

# LPC
#f0 = TFile('/eos/uscms/store/user/mjoyce/BTL/FNAL_TB_Mar2018/combined/Run902_to_953.root', 'READ') # all 5
#f1 = TFile('/eos/uscms/store/user/mjoyce/BTL/FNAL_TB_Mar2018/combined/bottombars_66V.root', 'READ') # low bias (66 V) bars 1 and 2
#f2 = TFile('/eos/uscms/store/user/mjoyce/BTL/FNAL_TB_Mar2018/combined/topbars_66V.root', 'READ') # low bias (66 V) bars 3, 4, and 5

t0 = f0.pulse
t1 = f1.pulse
#t2 = f2.pulse


if(args.test is not None):
    barClass(t1, 'bottomBars_66V', topDir, args.vetoOpt, True)
else:
    barClass(t0, 'all5exposure', topDir, args.vetoOpt)
    barClass(t1, 'bottomBars_66V', topDir, vetoOpt)
    barClass(t2, 'topBars_66V', topDir, vetoOpt)
