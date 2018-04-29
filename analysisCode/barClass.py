# !/usr/bin/python

#Author: Ben Tannenwald
#Date: April 16, 2018
#Purpose: Class for handling testbeam bar data

import os,sys, argparse
from ROOT import gROOT, TH1D, TFile, TTree, TChain, TCanvas, TH2D, TLegend, gStyle, TLatex, TProfile, TF1, TGraph, TMath, TPad, TLine, TObjArray

class barClass:
    def __init__(self, tree, runType, topDir, vetoOpt, test=False):
        # *** 0. Top-level options and objects
        self.tree = tree
        self.topDir = topDir
        self.signalThreshold = 0
        self.vetoThreshold = 0
        self.xBoundaries = []
        self.yBoundaries = []
        self.yIntegralOffset = 0
        self.setVarsByRunType(runType)
        self.vetoOpt = vetoOpt
        self.isTest = test
        self.fitPercentThreshold = 0.06
        #self.fitVoltageThreshold = 150 # in mV
        #self.fitVoltageForTiming = 200 # in mV
        self.fitVoltageThreshold = 30 # in mV
        self.fitVoltageForTiming = 50 # in mV #100 gave 150 ps
        self.fitMCPVoltageThreshold = 20 # in mV
        self.fitMCPVoltageForTiming = 40 # in mV
        self.fitSignalThreshold = 400 # in mV
        self.fitTimeWindow = 7
        self.fitMCPTimeWindow = 4
        self.fitFunction = "landau"
        self.histArray = TObjArray() 
        self.xSlice = 5 # in mm

        # *** 1. Define all histograms
        h_b1 = TH2D("h_b1", "h_b1", 40, -5, 35, 35, 0, 35)
        h_b2 = TH2D("h_b2", "h_b2", 40, -5, 35, 35, 0, 35)
        h_b3 = TH2D("h_b3", "h_b3", 40, -5, 35, 35, 0, 35)
        h_b4 = TH2D("h_b4", "h_b4", 40, -5, 35, 35, 0, 35)
        h_b5 = TH2D("h_b5", "h_b5", 40, -5, 35, 35, 0, 35)
        h_b1_t = TH2D("h_b1_t", "h_b1_t", 40, -5, 35, 35, 0, 35)
        h_b2_t = TH2D("h_b2_t", "h_b2_t", 40, -5, 35, 35, 0, 35)
        h_b3_t = TH2D("h_b3_t", "h_b3_t", 40, -5, 35, 35, 0, 35)
        h_b4_t = TH2D("h_b4_t", "h_b4_t", 40, -5, 35, 35, 0, 35)
        h_b5_t = TH2D("h_b5_t", "h_b5_t", 40, -5, 35, 35, 0, 35)
        
        h_mcp0_ch1 = TH1D("h_mcp0_ch1", "h_mcp0_ch1", 25, 0, 500);
        h_mcp0_ch3 = TH1D("h_mcp0_ch3", "h_mcp0_ch3", 25, 0, 500);
        h_mcp0_ch5 = TH1D("h_mcp0_ch5", "h_mcp0_ch5", 25, 0, 500);
        h_mcp1_ch10 = TH1D("h_mcp1_ch10", "h_mcp1_ch10", 25, 0, 500);
        h_mcp1_ch12 = TH1D("h_mcp1_ch12", "h_mcp1_ch12", 25, 0, 500);
        
        h_ch1_vs_ch2 = TH2D("h_ch1_vs_ch2", "h_ch1_vs_ch2", 22, 0, 1100, 22, 0, 1100)
        h_ch3_vs_ch4 = TH2D("h_ch3_vs_ch4", "h_ch3_vs_ch4", 22, 0, 1100, 22, 0, 1100)
        h_ch5_vs_ch6 = TH2D("h_ch5_vs_ch6", "h_ch5_vs_ch6", 22, 0, 1100, 22, 0, 1100)
        h_ch10_vs_ch11 = TH2D("h_ch10_vs_ch11", "h_ch10_vs_ch11", 22, 0, 1100, 22, 0, 1100)
        h_ch12_vs_ch13 = TH2D("h_ch12_vs_ch13", "h_ch12_vs_ch13", 22, 0, 1100, 22, 0, 1100)
        
        h_ch1_ch2_x_vs_ratio = TProfile("h_ch1_ch2_x_vs_ratio", "h_ch1_ch2_x_vs_ratio", 40, -5, 35, 0, 2)
        h_ch3_ch4_x_vs_ratio = TProfile("h_ch3_ch4_x_vs_ratio", "h_ch3_ch4_x_vs_ratio", 40, -5, 35, 0, 2)
        h_ch5_ch6_x_vs_ratio = TProfile("h_ch5_ch6_x_vs_ratio", "h_ch5_ch6_x_vs_ratio", 40, -5, 35, 0, 2)
        h_ch10_ch11_x_vs_ratio = TProfile("h_ch10_ch11_x_vs_ratio", "h_ch10_ch11_x_vs_ratio", 40, -5, 35, 0, 2)
        h_ch12_ch13_x_vs_ratio = TProfile("h_ch12_ch13_x_vs_ratio", "h_ch12_ch13_x_vs_ratio", 40, -5, 35, 0, 2)

        h_ch1_x_vs_amp = TProfile("h_ch1_x_vs_amp", "h_ch1_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch2_x_vs_amp = TProfile("h_ch2_x_vs_amp", "h_ch2_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch3_x_vs_amp = TProfile("h_ch3_x_vs_amp", "h_ch3_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch4_x_vs_amp = TProfile("h_ch4_x_vs_amp", "h_ch4_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch5_x_vs_amp = TProfile("h_ch5_x_vs_amp", "h_ch5_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch6_x_vs_amp = TProfile("h_ch6_x_vs_amp", "h_ch6_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch10_x_vs_amp = TProfile("h_ch10_x_vs_amp", "h_ch10_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch11_x_vs_amp = TProfile("h_ch11_x_vs_amp", "h_ch11_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch12_x_vs_amp = TProfile("h_ch12_x_vs_amp", "h_ch12_x_vs_amp", 40, -5, 35, 0, 1000)
        h_ch13_x_vs_amp = TProfile("h_ch13_x_vs_amp", "h_ch13_x_vs_amp", 40, -5, 35, 0, 1000)

        h_ch1_x_vs_time = TProfile("h_ch1_x_vs_time", "h_ch1_x_vs_time", 40, -5, 35, 0, 100)
        h_ch2_x_vs_time = TProfile("h_ch2_x_vs_time", "h_ch2_x_vs_time", 40, -5, 35, 0, 100)
        h_ch3_x_vs_time = TProfile("h_ch3_x_vs_time", "h_ch3_x_vs_time", 40, -5, 35, 0, 100)
        h_ch4_x_vs_time = TProfile("h_ch4_x_vs_time", "h_ch4_x_vs_time", 40, -5, 35, 0, 100)
        h_ch5_x_vs_time = TProfile("h_ch5_x_vs_time", "h_ch5_x_vs_time", 40, -5, 35, 0, 100)
        h_ch6_x_vs_time = TProfile("h_ch6_x_vs_time", "h_ch6_x_vs_time", 40, -5, 35, 0, 100)
        h_ch10_x_vs_time = TProfile("h_ch10_x_vs_time", "h_ch10_x_vs_time", 40, -5, 35, 0, 100)
        h_ch11_x_vs_time = TProfile("h_ch11_x_vs_time", "h_ch11_x_vs_time", 40, -5, 35, 0, 100)
        h_ch12_x_vs_time = TProfile("h_ch12_x_vs_time", "h_ch12_x_vs_time", 40, -5, 35, 0, 100)
        h_ch13_x_vs_time = TProfile("h_ch13_x_vs_time", "h_ch13_x_vs_time", 40, -5, 35, 0, 100)

        h_allChannel_timing = TH1D("h_allChannel_timing", "h_allChannel_timing", 60, 0, 60)
        h_allChannel_timingLogic = TH1D("h_allChannel_timingLogic", "h_allChannel_timingLogic", 4, 0, 4)
        h_allChannel_timingRes = TH1D("h_allChannel_timingRes", "h_allChannel_timingRes", 300, -1500, 1500)
        #h_allChannel_mcpRef_timingRes = TH1D("h_allChannel_mcpRef_timingRes", "h_allChannel_mcpRef_timingRes", 300, -1500, 1500)
        h_allChannel_mcpRef_timingRes = TH1D("h_allChannel_mcpRef_timingRes", "h_allChannel_mcpRef_timingRes", 350, -3500, 0)
        h_allChannel_x_vs_timingRes = TProfile("h_allChannel_x_vs_timingRes", "h_allChannel_x_vs_timingRes", 40, -5, 35, -1500, 1500)
        #h_allChannel_x_vs_mcpRef_timingRes = TProfile("h_allChannel_x_vs_mcpRef_timingRes", "h_allChannel_x_vs_mcpRef_timingRes", 40, -5, 35, -1500, 1500)
        h_allChannel_x_vs_mcpRef_timingRes = TProfile("h_allChannel_x_vs_mcpRef_timingRes", "h_allChannel_x_vs_mcpRef_timingRes", 40, -5, 35, -2600, -1800)
        #h_allChannel_fitSlope_vs_mcpRef_timingRes = TProfile("h_allChannel_fitSlope_vs_mcpRef_timingRes", "h_allChannel_fitSlope_vs_mcpRef_timingRes", 26, 50, 180, -2600, -1800)
        h_allChannel_fitSlope_vs_mcpRef_timingRes = TProfile("h_allChannel_fitSlope_vs_mcpRef_timingRes", "h_allChannel_fitSlope_vs_mcpRef_timingRes", 150, -3000, -1500, 60, 180)
        h_allChannel_mcpRef_timingRes_ampWalkCorrected = TH1D("h_allChannel_mcpRef_timingRes_ampWalkCorrected", "h_allChannel_mcpRef_timingRes_ampWalkCorrected", 350, -3500, 0)
        h_allChannel_mcpRef_ampWalkCorrection = TH1D("h_allChannel_mcpRef_ampWalkCorrection", "h_allChannel_mcpRef_ampWalkCorrection", 200, -100, 100)

        
        # *** 2. Add all histograms to array
        self.histArray.AddLast(h_b1)
        self.histArray.AddLast(h_b2)
        self.histArray.AddLast(h_b3)
        self.histArray.AddLast(h_b4)
        self.histArray.AddLast(h_b5)
        self.histArray.AddLast(h_b1_t)
        self.histArray.AddLast(h_b2_t)
        self.histArray.AddLast(h_b3_t)
        self.histArray.AddLast(h_b4_t)
        self.histArray.AddLast(h_b5_t)
        self.histArray.AddLast(h_mcp0_ch1)
        self.histArray.AddLast(h_mcp0_ch3)
        self.histArray.AddLast(h_mcp0_ch5)
        self.histArray.AddLast(h_mcp1_ch10)
        self.histArray.AddLast(h_mcp1_ch12)

        self.histArray.AddLast(h_ch1_vs_ch2)
        self.histArray.AddLast(h_ch3_vs_ch4)
        self.histArray.AddLast(h_ch5_vs_ch6)
        self.histArray.AddLast(h_ch10_vs_ch11)
        self.histArray.AddLast(h_ch12_vs_ch13)
        self.histArray.AddLast(h_ch1_ch2_x_vs_ratio)
        self.histArray.AddLast(h_ch3_ch4_x_vs_ratio)
        self.histArray.AddLast(h_ch5_ch6_x_vs_ratio)
        self.histArray.AddLast(h_ch10_ch11_x_vs_ratio)
        self.histArray.AddLast(h_ch12_ch13_x_vs_ratio)
        
        self.histArray.AddLast(h_ch1_x_vs_amp)
        self.histArray.AddLast(h_ch2_x_vs_amp)
        self.histArray.AddLast(h_ch3_x_vs_amp)
        self.histArray.AddLast(h_ch4_x_vs_amp)
        self.histArray.AddLast(h_ch5_x_vs_amp)
        self.histArray.AddLast(h_ch6_x_vs_amp)
        self.histArray.AddLast(h_ch10_x_vs_amp)
        self.histArray.AddLast(h_ch11_x_vs_amp)
        self.histArray.AddLast(h_ch12_x_vs_amp)
        self.histArray.AddLast(h_ch13_x_vs_amp)

        self.histArray.AddLast(h_ch1_x_vs_time)
        self.histArray.AddLast(h_ch2_x_vs_time)
        self.histArray.AddLast(h_ch3_x_vs_time)
        self.histArray.AddLast(h_ch4_x_vs_time)
        self.histArray.AddLast(h_ch5_x_vs_time)
        self.histArray.AddLast(h_ch6_x_vs_time)
        self.histArray.AddLast(h_ch10_x_vs_time)
        self.histArray.AddLast(h_ch11_x_vs_time)
        self.histArray.AddLast(h_ch12_x_vs_time)
        self.histArray.AddLast(h_ch13_x_vs_time)

        self.histArray.AddLast(h_allChannel_timing)
        self.histArray.AddLast(h_allChannel_timingRes)
        self.histArray.AddLast(h_allChannel_timingLogic)
        self.histArray.AddLast(h_allChannel_mcpRef_timingRes)
        self.histArray.AddLast(h_allChannel_x_vs_timingRes)
        self.histArray.AddLast(h_allChannel_x_vs_mcpRef_timingRes)
        self.histArray.AddLast(h_allChannel_fitSlope_vs_mcpRef_timingRes)
        self.histArray.AddLast(h_allChannel_mcpRef_timingRes_ampWalkCorrected)
        self.histArray.AddLast(h_allChannel_mcpRef_ampWalkCorrection)


        # *** 3. Canvases and style
        self.histArray = self.addSlicedHistograms(self.histArray, 1)
        self.histArray = self.addSlicedHistograms(self.histArray, 2)
        self.histArray = self.addSlicedHistograms(self.histArray, 3)
        self.histArray = self.addSlicedHistograms(self.histArray, 4)
        self.histArray = self.addSlicedHistograms(self.histArray, 5)

        # *** 4. Canvases and style
        self.c1 = TCanvas("c1", "c1", 800, 800)
        self.c2 = TCanvas("c2", "c2", 800, 800)
        self.c3 = TCanvas("c3", "c3", 800, 800)
        self.c4 = TCanvas("c4", "c4", 800, 800)

        gStyle.SetOptStat(0000)

        # *** 5. Make some directories if not already existent
        if not os.path.isdir(self.topDir):
            os.system( 'mkdir {0}'.format(self.topDir) )
        self.topDir = '{0}/{1}'.format(topDir, runType)
        if not os.path.isdir( '{0}/{1}'.format(topDir, runType) ):
            os.system( 'mkdir {0}'.format(self.topDir) )

        # *** 6. Run analysis
        self.loopEvents()

    # =============================

    def addSlicedHistograms(self, arr, barNum):
        """ function to add histograms sliced by x"""
        
        # calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        x = -5
        while x < 35:
            xLow = str(x).replace('-','n')
            xHigh = str(x+self.xSlice).replace('-','n')
            xBase = '{0}_to_{1}'.format(xLow, xHigh)
            #print xBase
            h_r = TH1D('h_ch{0}_timingRes_{1}'.format(rightSiPMchannel, xBase), 'h_ch{0}_timingRes_{1}'.format(rightSiPMchannel, xBase), 300, -1500, 1500)
            h_l = TH1D('h_ch{0}_timingRes_{1}'.format(leftSiPMchannel, xBase), 'h_ch{0}_timingRes_{1}'.format(leftSiPMchannel, xBase), 300, -1500, 1500)
            h_r_mcp = TH1D('h_ch{0}_mcpRef_timingRes_{1}'.format(rightSiPMchannel, xBase), 'h_ch{0}_mcpRef_timingRes_{1}'.format(rightSiPMchannel, xBase), 350, -3500, 0)
            h_l_mcp = TH1D('h_ch{0}_mcpRef_timingRes_{1}'.format(leftSiPMchannel, xBase), 'h_ch{0}_mcpRef_timingRes_{1}'.format(leftSiPMchannel, xBase), 350, -3500, 0)
            h_b = TH1D('h_b{0}_timingRes_{1}'.format(barNum, xBase), 'h_b{0}_timingRes_{1}'.format(barNum, xBase), 300, -1500, 1500)
            h_b_mcp = TH1D('h_b{0}_mcpRef_timingRes_{1}'.format(barNum, xBase), 'h_b{0}_mcpRef_timingRes_{1}'.format(barNum, xBase), 350, -3500, 0)
            arr.AddLast(h_r)
            arr.AddLast(h_l)
            arr.AddLast(h_r_mcp)
            arr.AddLast(h_l_mcp)
            arr.AddLast(h_b)
            arr.AddLast(h_b_mcp)
            x += self.xSlice

        return arr
    # =============================

    def setVarsByRunType(self, runType):
        """ set various constants as function of runType"""

        if runType == "all5exposure":
            self.signalThreshold = 800
            self.vetoThreshold   = 100
            self.xBoundaries = [-2, 6, 15, 24, 33]
            self.yBoundaries = [7.5, 12.5, 16.5, 20.5, 24.5]
            self.yIntegralOffset = 2.5
        elif runType == "bottomBars_66V":
            self.signalThreshold = 30
            self.vetoThreshold   = 30
            self.xBoundaries = [17, 19, 21, 23, 25]
            self.yBoundaries = [4.5, 9.5, 13.5, 16.5, 20.5] #5, 9, 13
            self.yIntegralOffset = 2.5
        elif runType == "topBars_66V":
            self.signalThreshold = 30
            self.vetoThreshold   = 30
            self.xBoundaries = [17, 19, 21, 23, 25]
            self.yBoundaries = [25, 25, 2.5, 7.5, 11.5] # 6-9, 1-4 (bins) 10-13
            self.yIntegralOffset = 1.5

    # =============================
    
    def returnVetoDecision(self, event, barNum, vetoOption):
        """ function to return veto decision given option vetoOption = 'None'/'singleAdj'/'doubleAdj'/'allAdj'/'all'"""

        if vetoOption == 'none':
            return False

        if vetoOption == 'singleAdj':
            if barNum == 1 and event.amp[3] < self.vetoThreshold:
                return False
            if barNum == 2 and event.amp[5] < self.vetoThreshold:
                return False
            if barNum == 3 and event.amp[10] < self.vetoThreshold:
                return False
            if barNum == 4 and event.amp[5] < self.vetoThreshold:
                return False
            if barNum == 5 and event.amp[10] < self.vetoThreshold:
                return False
        if vetoOption == 'doubleAdj':
            if barNum == 1 and event.amp[3] < self.vetoThreshold and event.amp[4] < self.vetoThreshold:
                return False
            if barNum == 2 and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold:
                return False
            if barNum == 3 and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold:
                return False
            if barNum == 4 and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold:
                return False
            if barNum == 5 and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold:
                return False

        if vetoOption == 'allAdj':
            if barNum == 1 and event.amp[3] < self.vetoThreshold and event.amp[4] < self.vetoThreshold:
                return False
            if barNum == 2 and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold and event.amp[1] < self.vetoThreshold and event.amp[2] < self.vetoThreshold:
                return False
            if barNum == 3 and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold and event.amp[3] < self.vetoThreshold and event.amp[4] < self.vetoThreshold:
                return False
            if barNum == 4 and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold and event.amp[12] < self.vetoThreshold and event.amp[13] < self.vetoThreshold:
                return False
            if barNum == 5 and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold:
                return False

        if vetoOption == 'all':
            if barNum == 1 and event.amp[3] < self.vetoThreshold and event.amp[4] < self.vetoThreshold and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold and event.amp[12] < self.vetoThreshold and event.amp[13] < self.vetoThreshold:
                return False
            if barNum == 2 and event.amp[1] < self.vetoThreshold and event.amp[2] < self.vetoThreshold and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold and event.amp[12] < self.vetoThreshold and event.amp[13] < self.vetoThreshold:
                return False
            if barNum == 3 and event.amp[1] < self.vetoThreshold and event.amp[2] < self.vetoThreshold and event.amp[3] < self.vetoThreshold and event.amp[4] < self.vetoThreshold and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold and event.amp[12] < self.vetoThreshold and event.amp[13] < self.vetoThreshold:
                return False
            if barNum == 4 and event.amp[1] < self.vetoThreshold and event.amp[2] < self.vetoThreshold and event.amp[3] < self.vetoThreshold and event.amp[4] < self.vetoThreshold and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold and event.amp[12] < self.vetoThreshold and event.amp[13] < self.vetoThreshold:
                return False
            if barNum == 5 and event.amp[1] < self.vetoThreshold and event.amp[2] < self.vetoThreshold and event.amp[3] < self.vetoThreshold and event.amp[4] < self.vetoThreshold and event.amp[5] < self.vetoThreshold and event.amp[6] < self.vetoThreshold and event.amp[10] < self.vetoThreshold and event.amp[11] < self.vetoThreshold:
                return False


        # if it gets here, means nothing returned, i.e. no logic was settled. veto event
        return True

    
    
    # =============================
    
    def returnChannelNumbers(self, barNum):
        """ single function to return mcp channel, time channel, and right/left SiPM channels for a given bar"""
        leftSiPMchannel  = -1
        rightSiPMchannel = -1
        mcpChannel       = -1
        timeChannel      = -1

        if barNum == 1 or barNum == 2 or barNum == 3:
            rightSiPMchannel = 2*barNum - 1
            leftSiPMchannel  = 2*barNum
            mcpChannel  = 0
            timeChannel = 0
                            
        elif barNum == 4 or barNum == 5:
            rightSiPMchannel = 2*barNum + 2
            leftSiPMchannel  = 2*barNum + 3
            mcpChannel = 9
            timeChannel = 1

        return rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel

    # =============================

    def fillChannelPlots(self, event, barNum, arr):
        """ function to fill bar-specific plots"""
        
        # calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        # logic to see if event should be vetoed based on signals in other bars
        doVetoEvent = self.returnVetoDecision(event, barNum, self.vetoOpt) # vetoOpt = none, singleAdj, doubleAdj, allAdj, all

        if (event.amp[rightSiPMchannel] > self.signalThreshold and event.amp[leftSiPMchannel] > self.signalThreshold and not doVetoEvent and abs(event.xSlope) < 0.0004 and abs(event.ySlope) < 0.0004):
            # leakage histograms and profiles
            arr.FindObject('h_b{0}'.format(barNum)).Fill(event.x_dut[2], event.y_dut[2])
            arr.FindObject('h_mcp{0}_ch{1}'.format(timeChannel, rightSiPMchannel)).Fill( event.amp[mcpChannel] )
            arr.FindObject('h_ch{0}_vs_ch{1}'.format(rightSiPMchannel, leftSiPMchannel)).Fill( event.amp[rightSiPMchannel], event.amp[leftSiPMchannel] )
            arr.FindObject('h_ch{0}_ch{1}_x_vs_ratio'.format(rightSiPMchannel, leftSiPMchannel)).Fill(event.x_dut[2], event.amp[rightSiPMchannel] / event.amp[leftSiPMchannel] )
            arr.FindObject('h_ch{0}_x_vs_amp'.format(rightSiPMchannel)).Fill(event.x_dut[2], event.amp[rightSiPMchannel] )
            arr.FindObject('h_ch{0}_x_vs_amp'.format(leftSiPMchannel)).Fill(event.x_dut[2], event.amp[leftSiPMchannel] )

            # test area for hit integral defintion
            if (abs(event.y_dut[2] - self.yBoundaries[barNum - 1]) <= self.yIntegralOffset and event.x_dut[2]>=self.xBoundaries[0] and event.x_dut[2]<=self.xBoundaries[ len(self.xBoundaries)-1 ]):
                arr.FindObject('h_b{0}_t'.format(barNum)).Fill(event.x_dut[2], event.y_dut[2])
                    
            # timing stuff
            #mipTime_R, fitSlope_R, mipTime_R_ampWalkCorrected = self.getTimingForChannel(event.time, event.channel, timeChannel, rightSiPMchannel, event.i_evt)
            #mipTime_L, fitSlope_L, mipTime_L_ampWalkCorrected = self.getTimingForChannel(event.time, event.channel, timeChannel, leftSiPMchannel, event.i_evt)
            fitStartTime_R, fitStartVoltage_R, fitSlope_R = self.getTimingForChannel(event.time, event.channel, timeChannel, rightSiPMchannel, event.i_evt)
            fitStartTime_L, fitStartVoltage_L, fitSlope_L = self.getTimingForChannel(event.time, event.channel, timeChannel, leftSiPMchannel, event.i_evt)
            mipTime_MCP = event.t_peak[mcpChannel]

            if fitSlope_R == 0:
                mipTime_R = 0
            else:
                mipTime_R = fitStartTime_R + (self.fitVoltageForTiming - fitStartVoltage_R)/fitSlope_R

            if fitSlope_L == 0:
                mipTime_L = 0
            else:
                mipTime_L = fitStartTime_L + (self.fitVoltageForTiming - fitStartVoltage_L)/fitSlope_L

            if mipTime_R != 0 and mipTime_L != 0:
                arr.FindObject('h_ch{0}_x_vs_time'.format(rightSiPMchannel)).Fill(event.x_dut[2], mipTime_R)
                arr.FindObject('h_ch{0}_x_vs_time'.format(leftSiPMchannel)).Fill(event.x_dut[2], mipTime_L)

                arr.FindObject('h_allChannel_timingLogic').Fill("Both",1)

                deltaT = 1000*(mipTime_L - mipTime_R) # multiple by 1000 to transfer from ns to ps
                arr.FindObject('h_allChannel_x_vs_timingRes').Fill(event.x_dut[2], deltaT)
                #if event.x_dut[2] > 5 and event.x_dut[2] < 25 and event.amp[rightSiPMchannel] > self.fitSignalThreshold: # keep it central
                arr.FindObject('h_allChannel_timingRes').Fill( deltaT ) 
                arr.FindObject('h_allChannel_timing').Fill(mipTime_L)
                arr.FindObject('h_allChannel_timing').Fill(mipTime_R)
                # *** fill x-slice plots
                x = -5
                while x < 35:
                    xLow = str(x).replace('-','n')
                    xHigh = str(x+self.xSlice).replace('-','n')
                    xBase = '{0}_to_{1}'.format(xLow, xHigh)
                    if event.x_dut[2] >= x and event.x_dut[2] < x+self.xSlice:
                        arr.FindObject('h_ch{0}_timingRes_{1}'.format(rightSiPMchannel, xBase)).Fill( deltaT )
                        arr.FindObject('h_ch{0}_timingRes_{1}'.format(leftSiPMchannel, xBase)).Fill( deltaT )
                        #print 'h_b{0}_timingRes_{1}'.format(barNum, xBase)
                        arr.FindObject('h_b{0}_timingRes_{1}'.format(barNum, xBase)).Fill( deltaT )
                    x += self.xSlice
                
                # *** Do timing resolution with MCP info ***
                if mipTime_MCP != 0 and event.amp[mcpChannel] > 80 and event.amp[mcpChannel] < 160:
                    deltaT_mcp = 1000*(((mipTime_R + mipTime_L)/2) - mipTime_MCP) # multiple by 1000 to transfer from ns to ps
                    arr.FindObject('h_allChannel_x_vs_mcpRef_timingRes').Fill(event.x_dut[2], deltaT_mcp)
                    #if event.x_dut[2] > 5 and event.x_dut[2] < 25 and event.amp[rightSiPMchannel] > self.fitSignalThreshold: # keep it central
                    arr.FindObject('h_allChannel_mcpRef_timingRes').Fill( deltaT_mcp )
                    arr.FindObject('h_allChannel_timing').Fill(mipTime_MCP)
                    #h_all_fitSlope_vs_mcpRef_timeRes.Fill(fitSlope_R, deltaT_mcp)
                    #h_all_fitSlope_vs_mcpRef_timeRes.Fill(fitSlope_L, deltaT_mcp)
                    arr.FindObject('h_allChannel_fitSlope_vs_mcpRef_timingRes').Fill(deltaT_mcp, fitSlope_R)
                    arr.FindObject('h_allChannel_fitSlope_vs_mcpRef_timingRes').Fill(deltaT_mcp, fitSlope_L)

                    # *** fill x-slice plots w/ mcp data
                    x = -5
                    while x < 35:
                        xLow = str(x).replace('-','n')
                        xHigh = str(x+self.xSlice).replace('-','n')
                        xBase = '{0}_to_{1}'.format(xLow, xHigh)
                        if event.x_dut[2] >= x and event.x_dut[2] < x+self.xSlice:
                            arr.FindObject('h_ch{0}_mcpRef_timingRes_{1}'.format(rightSiPMchannel, xBase)).Fill( deltaT_mcp )
                            arr.FindObject('h_ch{0}_mcpRef_timingRes_{1}'.format(leftSiPMchannel, xBase)).Fill( deltaT_mcp )
                            arr.FindObject('h_b{0}_mcpRef_timingRes_{1}'.format(barNum, xBase)).Fill( deltaT_mcp )
                        x += self.xSlice
                        
            
                    # ****   amp-walk corrected plots   ****
                    f_ampWalkCorrected = TF1()
                    if self.vetoOpt == 'singleAdj':
                        f_ampWalkCorrected = TF1("slope_ampWalkCorrected", "(-0.004)*x*x*x + (0.357)*x*x + (-13.681)*x + (-2077.244)") # from 50k run using singleAdj veto
                    if self.vetoOpt == 'doubleAdj':
                        f_ampWalkCorrected = TF1("slope_ampWalkCorrected", "(-0.031)*x + (23.392)") # from 50k run using doubleAdj veto

                    mipTime_R_ampWalkCorrected = fitStartTime_R + (self.fitVoltageForTiming - fitStartVoltage_R)/f_ampWalkCorrected.Eval(deltaT_mcp)
                    mipTime_L_ampWalkCorrected = fitStartTime_L + (self.fitVoltageForTiming - fitStartVoltage_L)/f_ampWalkCorrected.Eval(deltaT_mcp)
                    deltaT_mcp_ampWalkCorrected = 1000*(((mipTime_R_ampWalkCorrected + mipTime_L_ampWalkCorrected)/2) - mipTime_MCP) # multiple by 1000 to transfer from ns to ps
                    #print "fitSlope_R: {0}, walkSlope: {1}".format(fitSlope_R, f_ampWalkCorrected.Eval(deltaT_mcp))
                    #print "fitSlope_L: {0}, walkSlope: {1}".format(fitSlope_L, f_ampWalkCorrected.Eval(deltaT_mcp))
                    arr.FindObject('h_allChannel_mcpRef_ampWalkCorrection').Fill( 1000*(mipTime_R - mipTime_R_ampWalkCorrected) )
                    arr.FindObject('h_allChannel_mcpRef_ampWalkCorrection').Fill( 1000*(mipTime_L - mipTime_L_ampWalkCorrected) )
                    arr.FindObject('h_allChannel_mcpRef_timingRes_ampWalkCorrected').Fill( deltaT_mcp_ampWalkCorrected )

                if mipTime_R == mipTime_L and mipTime_MCP != 0 and event.amp[mcpChannel] > 80 and event.amp[mcpChannel] < 160:
                    print 'mipTime_R = {0}, mipTime_L = {1}, event: {2}'.format(mipTime_R, mipTime_L, event.i_evt)

            if mipTime_R != 0 and mipTime_L == 0:
                arr.FindObject('h_allChannel_timingLogic').Fill("R only",1)
            if mipTime_R == 0 and mipTime_L != 0:
                arr.FindObject('h_allChannel_timingLogic').Fill("L only",1)
            if mipTime_R == 0 and mipTime_L == 0:
                arr.FindObject('h_allChannel_timingLogic').Fill("None",1)

        return arr

    # =============================

    def getTimingForChannel(self, time, channel, drs_time, drs_channel, i_evt):
        """ function to calculate and return information about waveform for fitting"""
        voltageFromFunction = 0
        timeStep = 0
        evalFit = 0
        fitSlope = 0
        fitRes = 0
        slopeRes = 0
        i = 0
        l_time = []
        l_channel = []
        # *** 1. Make arrays of trace information
        while i < 1024:
            l_time.append(time[drs_time*1024 + i])
            l_channel.append(-1*channel[drs_channel*1024 + i])
            i+=1            

        # *** 2. Then store a TGraph --> why not in same step? because it doesn't work for unknown reasons
        i=0
        g = TGraph()       
        while i < 1024:
            g.SetPoint(i, l_time[i], l_channel[i])        
            i=i+1

        fn1 = TF1("fn1", self.fitFunction) # first degree polynomial --> this is just a choice atm
        startFit = self.getWaveformInfo_TOFPET(l_time, l_channel)
        timeWindow = self.fitTimeWindow
        if drs_channel == 0 or drs_channel == 9:
            startFit = self.getWaveformInfo_MCP(l_time, l_channel)
            timeWindow = self.fitMCPTimeWindow
        
        #if i_evt == 2703:
        #    print '1) {0}'.format(timeStep)

        if startFit > 1 and (startFit + timeWindow) < 1024: # protection against weird waveforms
            fn1.SetRange(l_time[startFit], l_time[startFit + timeWindow])
            g.Fit("fn1", "QR")

            fnR = g.GetFunction("fn1")
            # numerical solution for timing info
            fitStop = self.fitVoltageForTiming
            timeStep = l_time[startFit] # timestamp to start at (use startFit cuz... duh)

            if fnR.Eval(timeStep) > self.fitVoltageThreshold: # sometimes we need to push back start when first point of fit is waay beyond fitThresholdVoltage due to steep risetime
                #if i_evt == 2703:
                #    print "REWIND, evt {0}, timeStep: ({1:0.3f}, {2:0.3f})".format(i_evt, timeStep, fnR.Eval(timeStep))
                while fnR.Eval(timeStep) > self.fitVoltageThreshold: 
                    timeStep = timeStep - 0.001 
                    if timeStep < 0:
                        break
                #if i_evt == 2703:
                #    print "REWIND, evt {0}, evalFit: ({1:0.3f}, {2:0.3f})".format(i_evt, timeStep, fnR.Eval(timeStep))
            evalFit = timeStep

            if drs_channel == 0 or drs_channel == 9:
                fitStop = self.fitMCPVoltageForTiming
                
            #if i_evt == 2703:
            #    print '2) {0}'.format(timeStep)

            fitRes, fitSlope = self.numericalSolve(fnR, evalFit, 0.001, fitStop)
            if fitSlope != 0:
                timeDiff = (50.000 - fnR.Eval(evalFit))/fitSlope
                slopeRes = evalFit + timeDiff

            if i_evt%500 == 0:
                c5 = TCanvas("c5", "c5", 800, 800)
                p1 = TPad("p1", "p1", 0.0, 0.0, 1.0, 1.0)
                p1.Draw()
                p1.cd()
                g.Draw()
                
                # Make inset for leading edge
                ltxIn = TLatex()
                ltxIn.SetTextAlign(9)
                ltxIn.SetTextFont(62)
                ltxIn.SetTextSize(0.021)
                ltxIn.SetNDC()
                ltxIn.DrawLatex(0.66, 0.81, "Leading Edge Inset")
                p2 =  TPad("p1", "p1", 0.6, 0.5, 0.9, 0.8)
                p2.Draw()
                p2.cd()
                # produce inset graph directly from data
                inset = TGraph()       
                n_inset = 0
                i = 0
                while i < 1024:
                    if l_time[i] > timeStep - 2 and l_time[i] < timeStep + 0.3:
                        inset.SetPoint(n_inset, l_time[i], l_channel[i])        
                        n_inset += 1
                    i=i+1

                inset.Draw()
                #fnR.Draw("same")
                fnR.DrawF1(timeStep - 2.3, timeStep + 0.2, "same")

                l1= TLine(inset.GetXaxis().GetXmin(),0,inset.GetXaxis().GetXmax(),0);
                l1.SetLineStyle(2);
                l1.SetLineWidth(3);
                l1.SetLineColor(600+2);
                l1.Draw("same");

                c5.Update()
                c5.Print( "{0}/waveformPlusPol1Fit_Ch{1}_Evt{2}.png".format(self.topDir, drs_channel, i_evt) )

                print 'Fit Result = {0}, Fit Slope = {1}'.format(fitRes, fitSlope)
                print 'Lin Result = {0}'.format(slopeRes)
                print 'Lin Amp = {0}, Fit Amp = {1}'.format(fnR.Eval(slopeRes), fnR.Eval(fitRes))
                print 'Chi2 / NDF = {0:0.1f} / {1:0.1f} = {2:0.1f}'.format(fnR.GetChisquare(), fnR.GetNDF(), fnR.GetChisquare()/fnR.GetNDF() ) 

            #if i_evt == 2703:
            #    print '2) {0}'.format(fitRes)

            if (fitRes - evalFit) < 0.002:
                print "only one step, evt {0}, startFit: ({1:0.3f}, {2:0.3f}), fitted: ({3:0.3f}, {4:0.3f})".format(i_evt, evalFit, fnR.Eval(evalFit), fitRes, fnR.Eval(fitRes))

        #if i_evt == 2703:
        #    print '3) {0}'.format(fitRes)

        if fitRes <= evalFit or fitRes >= evalFit + 35: # something wonky --> not good
            return 0, 0, 0
        else: # good timing result!
            #if i_evt == 2703:
            #    print '4) {0}'.format(fitRes)
            #return fitRes, fitSlope
            #return slopeRes, fitSlope, result_ampWalkCorrected
            return evalFit, fnR.Eval(evalFit), fitSlope

    # =============================

    def getWaveformInfo(self, time, channel):
        """ function to calculate and return information about waveform for fitting"""

        # there are probably better ways to do all of this in python...
        maxADC = maxADC_stamp = threshold_stamp = -1
        i = 0

        # first find max
        for reading in channel:
            if reading > maxADC:
                maxADC = reading
                maxADC_stamp = i
            i += 1

        print maxADC, 'at', maxADC_stamp
        
        i=0
        # second find first reading above percent threshold i.e. place to start fit
        for reading in channel:
            if reading/maxADC > self.fitPercentThreshold and threshold_stamp == -1:
                threshold_stamp = i
            i += 1

        #print 'thresh met at', threshold_stamp, 'with', channel[threshold_stamp]

        return threshold_stamp

    # =============================

    def getWaveformInfo_TOFPET(self, time, channel):
        """ function to calculate and return information about waveform for fitting"""

        # there are probably better ways to do all of this in python...
        threshold_stamp = -1
        i = 0

        # find first reading above voltage threshold i.e. place to start fit
        for reading in channel:
            if abs(reading) > self.fitVoltageThreshold and threshold_stamp == -1:
                threshold_stamp = i
            i = i + 1

        #print 'thresh met at', threshold_stamp, 'with', channel[threshold_stamp]

        return threshold_stamp

    # =============================

    def getWaveformInfo_MCP(self, time, channel):
        """ function to calculate and return information about waveform for fitting"""

        # there are probably better ways to do all of this in python...
        threshold_stamp = -1
        i = 0

        # find first reading above voltage threshold i.e. place to start fit
        for reading in channel:
            if abs(reading) > self.fitMCPVoltageThreshold and threshold_stamp == -1:
                threshold_stamp = i
            i = i + 1

        #print 'thresh met at', threshold_stamp, 'with', channel[threshold_stamp]

        return threshold_stamp

    # =============================
    
    def loopEvents(self):
        """ function looping over all events in file"""

        print self.tree.GetEntries()
        nTotal=0

        for event in self.tree:        
            if nTotal > 1000 and self.isTest:
                break

            nTotal += 1
            if (nTotal % 10000 == 0):
                print nTotal, "processed"

            
            # fill bar information
            self.histArray = self.fillChannelPlots(event, 1, self.histArray) # bar 1
            self.histArray = self.fillChannelPlots(event, 2, self.histArray) # bar 2
            self.histArray = self.fillChannelPlots(event, 3, self.histArray) # bar 3
            self.histArray = self.fillChannelPlots(event, 4, self.histArray) # bar 4
            self.histArray = self.fillChannelPlots(event, 5, self.histArray) # bar 5

       # end filling loop     

        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b1'), 1)
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b2'), 2)
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b3'), 3)
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b4'), 4)
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b5'), 5)
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b1_t'), 1, "test")
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b2_t'), 2, "test")
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b3_t'), 3, "test")
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b4_t'), 4, "test")
        self.draw2Dbar(self.c1, self.histArray.FindObject('h_b5_t'), 5, "test")
        
        self.c2.cd()
        self.c2.SetLeftMargin(0.15);
        self.c2.SetRightMargin(0.05);
        self.c2.SetBottomMargin(0.10);
        self.c2.SetTopMargin(0.05);
        # kBlack == 1, kRed == 632, kBlue == 600, kGreen == 416, kMagenta == 616
        self.histArray.FindObject('h_mcp0_ch1').SetLineColor(1) # kBlack
        self.histArray.FindObject('h_mcp0_ch3').SetLineColor(600) # kBlue
        self.histArray.FindObject('h_mcp0_ch5').SetLineColor(632) # kRed
        self.histArray.FindObject('h_mcp1_ch10').SetLineColor(416+2) # kGreen+2
        self.histArray.FindObject('h_mcp1_ch12').SetLineColor(616-3) # kMagenta-3
        
        self.histArray.FindObject('h_mcp0_ch1').SetLineWidth(3) 
        self.histArray.FindObject('h_mcp0_ch3').SetLineWidth(3) 
        self.histArray.FindObject('h_mcp0_ch5').SetLineWidth(3) 
        self.histArray.FindObject('h_mcp1_ch10').SetLineWidth(3) 
        self.histArray.FindObject('h_mcp1_ch12').SetLineWidth(3) 
        
        self.histArray.FindObject('h_mcp1_ch12').SetYTitle("Noramlized Entries / 20 mV")
        self.histArray.FindObject('h_mcp1_ch12').SetXTitle("MCP Amplitude [mV]")
        self.histArray.FindObject('h_mcp1_ch12').SetTitle("")
        self.histArray.FindObject('h_mcp1_ch12').DrawNormalized()
        self.histArray.FindObject('h_mcp1_ch12').GetYaxis().SetRangeUser(0,1.6)
        
        self.histArray.FindObject('h_mcp0_ch3').DrawNormalized("same")
        self.histArray.FindObject('h_mcp0_ch5').DrawNormalized("same")
        self.histArray.FindObject('h_mcp1_ch10').DrawNormalized("same")
        self.histArray.FindObject('h_mcp0_ch1').DrawNormalized("same")
        
        leg = TLegend(0.5, 0.4, .85, .7);
        leg.AddEntry(self.histArray.FindObject('h_mcp0_ch1'), "MCP Amplitude: Bar 1 Signal", "l");
        leg.AddEntry(self.histArray.FindObject('h_mcp0_ch3'), "MCP Amplitude: Bar 2 Signal", "l");
        leg.AddEntry(self.histArray.FindObject('h_mcp0_ch5'), "MCP Amplitude: Bar 3 Signal", "l");
        leg.AddEntry(self.histArray.FindObject('h_mcp1_ch10'), "MCP Amplitude: Bar 4 Signal", "l");
        leg.AddEntry(self.histArray.FindObject('h_mcp1_ch12'), "MCP Amplitude: Bar 5 Signal", "l");
        leg.Draw("same");
        
        self.c2.Print("{0}/mcp_amplitudes.png".format(self.topDir) )

        
        self.drawLvsRinBar(self.c3, self.histArray.FindObject('h_ch1_vs_ch2'), 1)
        self.drawLvsRinBar(self.c3, self.histArray.FindObject('h_ch3_vs_ch4'), 2)
        self.drawLvsRinBar(self.c3, self.histArray.FindObject('h_ch5_vs_ch6'), 3)
        self.drawLvsRinBar(self.c3, self.histArray.FindObject('h_ch10_vs_ch11'), 4)
        self.drawLvsRinBar(self.c3, self.histArray.FindObject('h_ch12_vs_ch13'), 5)
        
        
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch1_ch2_x_vs_ratio'), 1, 'Right/Left')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch3_ch4_x_vs_ratio'), 2, 'Right/Left')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch5_ch6_x_vs_ratio'), 3, 'Right/Left')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch10_ch11_x_vs_ratio'), 4, 'Right/Left')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch12_ch13_x_vs_ratio'), 5, 'Right/Left')

        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch1_x_vs_amp'), 1, 'Bar 1 Amplitude (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch2_x_vs_amp'), 1, 'Bar 1 Amplitude (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch3_x_vs_amp'), 2, 'Bar 2 Amplitude (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch4_x_vs_amp'), 2, 'Bar 2 Amplitude (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch5_x_vs_amp'), 3, 'Bar 3 Amplitude (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch6_x_vs_amp'), 3, 'Bar 3 Amplitude (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch10_x_vs_amp'), 4, 'Bar 4 Amplitude (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch11_x_vs_amp'), 4, 'Bar 4 Amplitude (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch12_x_vs_amp'), 5, 'Bar 5 Amplitude (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch13_x_vs_amp'), 5, 'Bar 5 Amplitude (Left SiPM)')

        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch1_x_vs_time'), 1, 'Bar 1 Time (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch2_x_vs_time'), 1, 'Bar 1 Time (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch3_x_vs_time'), 2, 'Bar 2 Time (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch4_x_vs_time'), 2, 'Bar 2 Time (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch5_x_vs_time'), 3, 'Bar 3 Time (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch6_x_vs_time'), 3, 'Bar 3 Time (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch10_x_vs_time'), 4, 'Bar 4 Time (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch11_x_vs_time'), 4, 'Bar 4 Time (Left SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch12_x_vs_time'), 5, 'Bar 5 Time (Right SiPM)')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_ch13_x_vs_time'), 5, 'Bar 5 Time (Left SiPM)')

        self.c4.cd()
        self.histArray.FindObject('h_allChannel_timing').Draw()
        self.c4.Print( "{0}/h_allChannel_timing.png".format(self.topDir) )

        self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_timingRes'), False)
        self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_mcpRef_timingRes'), True)
        self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_mcpRef_timingRes_ampWalkCorrected'), True, True)

        self.c4.cd()
        self.histArray.FindObject('h_allChannel_mcpRef_ampWalkCorrection').SetXTitle("Amp Walk Correction [ps]")
        self.histArray.FindObject('h_allChannel_mcpRef_ampWalkCorrection').SetYTitle("Entries / 1 ps")
        self.histArray.FindObject('h_allChannel_mcpRef_ampWalkCorrection').Draw()
        self.c4.Print( "{0}/h_allChannel_mcpRef_ampWalkCorrection.png".format(self.topDir) )
        
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_allChannel_x_vs_timingRes'), 0, 't_{right SiPM} - t_{left SiPM}')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_allChannel_x_vs_mcpRef_timingRes'), 0, '(t_{right SiPM} + t_{left SiPM})/2 - t_{MCP}')
        self.drawSingleProfile(self.c4, self.histArray.FindObject('h_allChannel_fitSlope_vs_mcpRef_timingRes'), 0, '(t_{right SiPM} + t_{left SiPM})/2 - t_{MCP}', opt='slope')

        self.c4.cd()
        self.histArray.FindObject('h_allChannel_timingLogic').Draw("TEXT")
        self.c4.Print( "{0}/h_allChannel_timingLogic.png".format(self.topDir) )

        self.drawTimingResXSlices(self.c4, self.histArray, 1, usingMCP=False)
        self.drawTimingResXSlices(self.c4, self.histArray, 1, usingMCP=True)
        self.drawTimingResXSlices(self.c4, self.histArray, 2, usingMCP=False)
        self.drawTimingResXSlices(self.c4, self.histArray, 2, usingMCP=True)
        self.drawTimingResXSlices(self.c4, self.histArray, 3, usingMCP=False)
        self.drawTimingResXSlices(self.c4, self.histArray, 3, usingMCP=True)
        self.drawTimingResXSlices(self.c4, self.histArray, 4, usingMCP=False)
        self.drawTimingResXSlices(self.c4, self.histArray, 4, usingMCP=True)
        self.drawTimingResXSlices(self.c4, self.histArray, 5, usingMCP=False)
        self.drawTimingResXSlices(self.c4, self.histArray, 5, usingMCP=True)

        # === function graveyard. keep for reference"
        #self.drawXquadrants(self.c4, self.histArray.FindObject('h_ch1_ch2_ratio_x1, self.histArray.FindObject('h_ch1_ch2_ratio_x2, self.histArray.FindObject('h_ch1_ch2_ratio_x3, self.histArray.FindObject('h_ch1_ch2_ratio_x4, 1, 'Right/Left')
        #self.drawXquadrants(self.c4, self.histArray.FindObject('h_ch1_x1, self.histArray.FindObject('h_ch1_x2, self.histArray.FindObject('h_ch1_x3, self.histArray.FindObject('h_ch1_x4, 1, 'Bar 1 [Right SiPM]')
        
    # =============================

    def drawResolutionPlot(self, c0, h0, isMCPref, isCorrected=False):
        """function to fit resolution plot with gaussian and print resutls"""
        xTitle = 't_{right SiPM} - t_{left SiPM} [ps]'
        plotName = 'h_allChannel_timingRes'
        xMax = 50
        xMin = -450
        if isMCPref:
            xTitle = '(t_{right SiPM} + t_{left SiPM})/2 - t_{MCP} [ps]'
            plotName = 'h_allChannel_mcpRef_timingRes'
            if isCorrected:
                xTitle = 'Amp-Walk Corrected (t_{right SiPM} + t_{left SiPM})/2 - t_{MCP} [ps]'
                plotName = 'h_allChannel_mcpRef_timingRes_ampWalkCorrected'

            xMax = -2100
            xMin = -2500

        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.10);
        h0.SetTitle("")
        h0.SetXTitle(xTitle)
        h0.SetYTitle("Entries / 10 ps")
        h0.Draw()
        f_res = TF1("f_res", "gaus") # gaussian
        f_res.SetRange(xMin, xMax)
        h0.Fit("f_res", "R") # should be "R" to impose range
        print "a"
        ltx1 = TLatex()
        ltx1.SetTextAlign(9)
        ltx1.SetTextFont(62)
        ltx1.SetTextSize(0.035)
        ltx1.SetNDC()
        ltx1.DrawLatex(0.20, 0.85, "UNBELIEVABLY PRELIMINARY")
        ltx1 = TLatex()
        ltx1.SetTextAlign(9)
        ltx1.SetTextFont(62)
        ltx1.SetTextSize(0.035)
        ltx1.SetNDC()
        ltx1.DrawLatex(0.70, 0.55, '#sigma_{{t}} = {0:0.1f} ps'.format(f_res.GetParameter(2)))
        c0.Print( "{0}/{1}.png".format(self.topDir, plotName) )
        

    # =============================

    def draw2Dbar(self, c0, h0, barNum, test=""):
        """ function to recieve canvas (c0), histogram (h0), and name for 2D bar plot"""
        c0.cd()
        h0.SetTitle( "Bar {0}".format(barNum) )
        h0.SetXTitle("X [mm]")
        h0.SetYTitle("Y [mm]")
        h0.Draw("colz")
        
        #h0.Integral(x1, x2, y1, y2)
        #7-10, 12-14, 16
        #ymin = 4*barNum + self.yIntegralOffset
        #ymax = 4*barNum + self.yIntegralOffset + 4
        ymin = int(self.yBoundaries[barNum -1] - self.yIntegralOffset)
        ymax = int(self.yBoundaries[barNum -1] + self.yIntegralOffset)
        window = h0.Integral(4, 38, ymin, ymax)
        
        if h0.GetEntries() > 0:
            print "In bar {0}: {1}\t Total: {2}\t % in Bar: {3}".format(barNum, window, h0.GetEntries(), window/h0.GetEntries())
            ltx1 = TLatex()
            ltx1.SetTextAlign(9)
            ltx1.SetTextFont(62)
            ltx1.SetTextSize(0.035)
            ltx1.SetNDC()
            ltx1.DrawLatex(0.15, 0.80, ("%Hits in Bar: {0:0.1f}%".format(100*window/h0.GetEntries())) )
            ltx1.DrawLatex(0.15, 0.84, ("N_{{Hits}} in Bar: {0}".format(window)) )
        else:
            print "no entries"

        c0.Print( "{0}/bar{1}_python{2}.png".format(self.topDir, barNum, test) )

    # =============================

    def drawLvsRinBar(self, c0, h0, barNum):
        """ function to recieve canvas (c0), histogram (h0), and name for 2D bar plot"""
        c0.cd()
        c0.SetLeftMargin(0.15);
        h0.SetTitle( "Bar {0}: R vs L Amplitudes".format(barNum) )
        h0.SetXTitle("Right SiPM [mV]")
        h0.SetYTitle("Left SiPM [mV]")
        h0.Draw("colz")
        
        c0.Print( "{0}/bar{1}_rightVleft.png".format(self.topDir, barNum) )

    # =============================

    def drawXquadrants(self, c0, h_x1, h_x2, h_x3, h_x4, barNum, varName):
        """ function to recieve canvas (c0), histograms (h_xN [N=1,2,3,4], and bar number"""
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.10);
        # kBlack == 1, kRed == 632, kBlue == 600, kGreen == 416, kMagenta == 616
        h_x1.SetLineColor(1) # kBlack
        h_x2.SetLineColor(600) # kBlue
        h_x3.SetLineColor(632) # kRed
        h_x4.SetLineColor(416+2) # kGreen+2
        
        h_x1.SetLineWidth(3)
        h_x2.SetLineWidth(3)
        h_x3.SetLineWidth(3)
        h_x4.SetLineWidth(3)
        
        if varName == 'Right/Left':
            xtitle = "Right SiPM Amplitude / Left SiPM Amplitude"
        else:
            xtitle = '{0} Amplitude'.format(varName)
        
        h_x4.SetYTitle("Noramlized Entries / Bin")
        h_x4.SetXTitle(xtitle)
        h_x4.SetXTitle(xtitle)
        h_x4.SetTitle("Bar {0}".format(barNum))
        h_x4.DrawNormalized()

        if h_x4.GetBinContent(h_x4.GetMaximumBin()) > h_x1.GetBinContent( h_x1.GetMaximumBin()):
            maxVal = h_x4.GetBinContent( h_x4.GetMaximumBin()) 
        else:
            maxVal = h_x1.GetBinContent( h_x1.GetMaximumBin())

        h_x4.GetYaxis().SetRangeUser(0,1.6*maxVal)
        h_x4.DrawNormalized()
        
        h_x3.DrawNormalized("same")
        h_x2.DrawNormalized("same")
        h_x1.DrawNormalized("same")
        
        leg = TLegend(0.6, 0.5, .85, .7);
        leg.AddEntry(h_x1, "{0} Amplitude: Q1".format(varName), "l");
        leg.AddEntry(h_x2, "{0} Amplitude: Q2".format(varName), "l");
        leg.AddEntry(h_x3, "{0}: Q3".format(varName), "l");
        leg.AddEntry(h_x4, "{0} Amplitude: Q4".format(varName), "l");
        leg.Draw("same");

        if varName == 'Right/Left':
            filename = '{0}/bar{1}_RL_ratio_x.png'.format(self.topDir, barNum)
        else:
            filename = '{0}/{1}_x.png'.format(self.topDir, varName.replace(' ', '_').replace('[','').replace(']',''))

        c0.Print(filename)

    # =============================

    def drawSingleProfile(self, c0, h_p, barNum, varName, opt='x'):
        """ function to recieve canvas (c0), TProfile (h_p), and bar number"""
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.10);
        
        if varName == 'Right/Left':
            ytitle = "Right SiPM Amplitude / Left SiPM Amplitude"
            filename = '{0}/bar{1}_RL_ratio_x.png'.format(self.topDir, barNum)
        elif 'Amplitude' in varName:
            ytitle = "{0} [mV]".format(varName)
        else:
            ytitle = '{0} [ps]'.format(varName)
            if opt == 'slope':
                ytitle = "Fit Slope"

        if opt =='x':
            h_p.SetXTitle("X [mm]")
        elif opt == 'slope':
            h_p.SetXTitle = '{0} [ps]'.format(varName)
        else:
            h_p.SetXTitle("barf")

        h_p.SetYTitle(ytitle)
        h_p.SetTitle("Bar {0}".format(barNum))
        if barNum == 0:
            h_p.SetTitle("All Bars".format(barNum))      


        if "t_{" in varName:

            filename = '{0}/allChannel_timing_vs_{1}.png'.format(self.topDir, opt)
            if 'MCP' in varName and opt == 'x':
                h_p.SetMinimum(-1800)
                h_p.SetMaximum(-2600)
                h_p.Draw()

            elif 'MCP' in varName and opt == 'slope':
                h_p.SetMinimum(60)
                h_p.SetMaximum(180)
                xMin = -2800 
                xMax = -1800 
                h_p.Draw()
                
                h_p1 = h_p.Clone("h_p1")
                h_p2 = h_p.Clone("h_p2")
                h_p3 = h_p.Clone("h_p3")
                h_p4 = h_p.Clone("h_p4")
                h_p5 = h_p.Clone("h_p5")
                ampWalk_pol1 = TF1('ampWalk_pol1', 'pol1') # first degree polynomial --> this is just a choice atm
                ampWalk_pol2 = TF1('ampWalk_pol2', 'pol2')
                ampWalk_pol3 = TF1('ampWalk_pol3', 'pol3')
                ampWalk_pol4 = TF1('ampWalk_pol4', 'pol4')
                ampWalk_pol5 = TF1('ampWalk_pol5', 'pol5')

                c_pol1 = 632 # kRed
                c_pol2 = 600 # kBlue
                c_pol3 = 416+2 # kGreen+2
                c_pol4 = 616-3 # kMagenta-3
                c_pol5 = 1 # kBlack

                h_p1.Fit('ampWalk_pol1', 'R', '', xMin, xMax)
                h_p1.GetFunction('ampWalk_pol1').SetLineColor(c_pol1)
                h_p2.Fit('ampWalk_pol2', 'R', '', xMin, xMax)
                h_p2.GetFunction('ampWalk_pol2').SetLineColor(c_pol2)
                h_p3.Fit('ampWalk_pol3', 'R', '', xMin, xMax)
                h_p3.GetFunction('ampWalk_pol3').SetLineColor(c_pol3)
                h_p4.Fit('ampWalk_pol4', 'R', '', xMin, xMax)
                h_p4.GetFunction('ampWalk_pol4').SetLineColor(c_pol4) 
                h_p5.Fit('ampWalk_pol5', 'R', '', xMin, xMax)
                h_p5.GetFunction('ampWalk_pol5').SetLineColor(c_pol5) 

                c0.cd()
                h_p1.Draw()
                h_p2.Draw("same")
                h_p3.Draw("same")
                h_p4.Draw("same")
                h_p5.Draw("same")
                
                ltx1 = TLatex()
                ltx1.SetTextAlign(9)
                ltx1.SetTextFont(62)
                ltx1.SetTextSize(0.025)
                ltx1.SetNDC()

                ltx1.SetTextColor(c_pol1)
                ltx1.DrawLatex(.6, .73, 'pol1: #chi^{{2}}/NDF = {0:0.2f}'.format(ampWalk_pol1.GetChisquare() / ampWalk_pol1.GetNDF()) )
                if self.vetoOpt == 'doubleAdj':
                    ltx1.DrawLatex(.40, .8, 'y = {0:0.3f}x + {1:0.3f}'.format(ampWalk_pol1.GetParameter(1), ampWalk_pol1.GetParameter(0)) )
                ltx1.SetTextColor(c_pol2)
                ltx1.DrawLatex(.6, .70, 'pol2: #chi^{{2}}/NDF = {0:0.2f}'.format(ampWalk_pol2.GetChisquare() / ampWalk_pol2.GetNDF()) )
                if self.vetoOpt == 'singleAdj':
                    ltx1.DrawLatex(.40, .8, 'y = {0:0.6f}x^{{2}} + {1:0.3f}x + {2:0.1f}'.format(ampWalk_pol2.GetParameter(2), ampWalk_pol2.GetParameter(1), ampWalk_pol2.GetParameter(0)) )
                ltx1.SetTextColor(c_pol3)
                ltx1.DrawLatex(.6, .67, 'pol3: #chi^{{2}}/NDF = {0:0.2f}'.format(ampWalk_pol3.GetChisquare() / ampWalk_pol3.GetNDF()) )
                ltx1.SetTextColor(c_pol4)
                ltx1.DrawLatex(.6, .64, 'pol4: #chi^{{2}}/NDF = {0:0.2f}'.format(ampWalk_pol4.GetChisquare() / ampWalk_pol4.GetNDF()) )
                ltx1.SetTextColor(c_pol5)
                ltx1.DrawLatex(.6, .61, 'pol5: #chi^{{2}}/NDF = {0:0.2f}'.format(ampWalk_pol5.GetChisquare() / ampWalk_pol5.GetNDF()) )
                
                #print 'y = {0:0.2f}*x + {1:0.2f}'.format(ampWalk.GetParameter(1), ampWalk.GetParameter(0))
                #h_p.Print("all")
                filename = '{0}/allChannel_mcpRef_timing_vs_{1}.png'.format(self.topDir, opt)
        else:
            filename = '{0}/{1}_{2}.png'.format(self.topDir, varName.replace(' ', '_').replace('(','').replace(')',''), opt)
            h_p.Draw()


        c0.Print(filename)

    # =============================

    def drawTimingResXSlices(self, c0, arr, barNum, usingMCP=False):
        """ function to make pretty canvas of x-slices of timing res """
        
        # calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        hname = 'h_b{0}_timingRes'.format(barNum)
        if usingMCP:
            hname = 'h_b{0}_mcpRef_timingRes'.format(barNum)

        leg = TLegend(0.68, 0.68, .93, .88);
        legendCol = [1, 600, 632, 416+2, 616-3, 800+7, 432-3, 900-3] # kBlack, kBlue, kRed, kGreen+2, kMagenta-3, kOrange+7, kCyan-3, kPink-3
        nSlices = 0
        x = -5
        c0.cd()

        while x < 35:
            xLow = str(x).replace('-','n')
            xHigh = str(x+self.xSlice).replace('-','n')
            xBase = '{0}_to_{1}'.format(xLow, xHigh)
            #h_b = TH1D('{0}_{1}'.format(hname, xBase), '{0}_{1}'.format(hname, xBase), 300, -1500, 1500)
            h_b = arr.FindObject('{0}_{1}'.format(hname, xBase))
            h_b.SetLineColor(legendCol[nSlices])
            leg.AddEntry(h_b, xBase.replace('n','-').replace('_',' '), "l")
            if nSlices == 0:
                h_b.GetYaxis().SetRangeUser(0, 1.4*h_b.GetBinContent(h_b.GetMaximumBin()))
                h_b.Draw()
            else:
                h_b.Draw("same")

            nSlices += 1
            x += self.xSlice
            
        leg.Draw("same")
        c0.Print('{0}/{1}_byXslice.png'.format(self.topDir, hname))

    # =============================

    def numericalSolve(self, func, x_eval, x_step, y_end):
        """ function to recieve TF1, find + return numerical solution by incrementing x_eval by x_step to point where function equals y_end"""
        y_calc = 0
        x_start = x_eval

        while abs(y_calc) < abs(y_end):
            store = y_calc
            y_calc = func.Eval(x_eval)
            #if i_evt == 2703:
            #    print 'old: {0}, new: {1}'.format(store, voltageFromFunction)
            x_eval += x_step
            if x_eval > x_start + 35: # scan over window of 35 ps
                break
                
                
        f1 = func.Eval(x_start)
        f2 = func.Eval(x_eval)
        fitSlope = (f2 - f1)/(x_eval - x_start)

        return x_eval, fitSlope

    # =============================
