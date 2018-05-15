# !/usr/bin/python

#Author: Ben Tannenwald
#Date: April 16, 2018
#Purpose: Class for handling testbeam bar data

import os,sys, argparse, ROOT
from ROOT import gROOT, TH1D, TFile, TTree, TChain, TCanvas, TH2D, TLegend, gStyle, TLatex, TProfile, TF1, TGraph, TMath, TPad, TLine, TObjArray

class barClass:
    def __init__(self, tree, runType, topDir, vetoOpt, doTiming=True, test=False, batch=False):
        # *** 0. Top-level options and objects
        self.tree = tree
        self.topDir = topDir
        self.doTiming = doTiming
        self.signalThreshold = 0
        self.vetoThreshold = 0
        self.xBoundaries = []
        self.yBoundaries = []
        self.yIntegralOffset = 0
        self.setVarsByRunType(runType)
        self.runType = runType
        self.vetoOpt = vetoOpt
        self.isTest = test
        self.runBatch = batch
        self.fitVoltageThreshold = 30 # in mV
        self.fitVoltageForTiming = 50 # in mV #100 gave 150 ps
        self.fitMCPVoltageThreshold = 20 # in mV
        self.fitMCPVoltageForTiming = 40 # in mV
        self.fitSignalThreshold = 400 # in mV
        self.fitTimeWindow = 4
        self.fitMCPTimeWindow = 4
        self.fitFunction = "gaus"
        self.histArray = TObjArray() 
        self.xSlice = 5 # in mm
        self.slopeSlice = 15 # arb units
        self.peakFitFunction = "landau"
        self.peakFitRiseThreshold = 400 # in mV
        self.peakFitFallThreshold = 850 # in mV
        self.peakFitVoltageVeto = 875 # in mV
        self.fitPercentThreshold = 0.06
        self.fitPercentForTiming = 0.10
                
        # *** 1. Define all histograms
        # ** A. Leakage histograms and general position info of hits
        self.histArray = self.addLeakageHistograms(self.histArray)
        # ** B. Add profiles
        self.histArray = self.addProfiles(self.histArray)
        # ** C. Add x-sliced histograms
        self.histArray = self.addSlicedHistograms(self.histArray)
        # ** D. Add timing histograms --> need more automatic way to incorporate...
        h_allChannel_timing = TH1D("h_allChannel_timing", "h_allChannel_timing", 60, 0, 60)
        h_allChannel_fracFit_timing = TH1D("h_allChannel_fracFit_timing", "h_allChannel_fracFit_timing", 60, 0, 60)
        h_allChannel_timingLogic = TH1D("h_allChannel_timingLogic", "h_allChannel_timingLogic", 4, 0, 4)
        h_allChannel_timingRes = TH1D("h_allChannel_timingRes", "h_allChannel_timingRes", 300, -1500, 1500)
        h_allChannel_fracFit_timingRes = TH1D("h_allChannel_fracFit_timingRes", "h_allChannel_fracFit_timingRes", 300, -1500, 1500)
        #h_allChannel_mcpRef_timingRes = TH1D("h_allChannel_mcpRef_timingRes", "h_allChannel_mcpRef_timingRes", 300, -1500, 1500)
        h_allChannel_mcpRef_timingRes = TH1D("h_allChannel_mcpRef_timingRes", "h_allChannel_mcpRef_timingRes", 350, -3500, 0)
        h_allChannel_mcpRef_fracFit_timingRes = TH1D("h_allChannel_mcpRef_fracFit_timingRes", "h_allChannel_mcpRef_fracFit_timingRes", 350, -3500, 0)
        h_allChannel_x_vs_timingRes = TProfile("h_allChannel_x_vs_timingRes", "h_allChannel_x_vs_timingRes", 40, -5, 35, -1500, 1500)
        #h_allChannel_x_vs_mcpRef_timingRes = TProfile("h_allChannel_x_vs_mcpRef_timingRes", "h_allChannel_x_vs_mcpRef_timingRes", 40, -5, 35, -1500, 1500)
        h_allChannel_x_vs_mcpRef_timingRes = TProfile("h_allChannel_x_vs_mcpRef_timingRes", "h_allChannel_x_vs_mcpRef_timingRes", 40, -5, 35, -2600, -1800)
        h_allChannel_fitSlope_vs_mcpRef_timingRes = TProfile("h_allChannel_fitSlope_vs_mcpRef_timingRes", "h_allChannel_fitSlope_vs_mcpRef_timingRes", 26, 50, 180, -2600, -1800)
        #h_allChannel_fitSlope_vs_mcpRef_timingRes = TProfile("h_allChannel_fitSlope_vs_mcpRef_timingRes", "h_allChannel_fitSlope_vs_mcpRef_timingRes", 150, -3000, -1500, 60, 180)
        h_allChannel_mcpRef_timingRes_ampWalkCorrected = TH1D("h_allChannel_mcpRef_timingRes_ampWalkCorrected", "h_allChannel_mcpRef_timingRes_ampWalkCorrected", 350, -3500, 0)
        h_allChannel_mcpRef_ampWalkCorrection = TH1D("h_allChannel_mcpRef_ampWalkCorrection", "h_allChannel_mcpRef_ampWalkCorrection", 200, -100, 100)
        h_allChannel_ampFit_percentError = TH1D("h_allChannel_ampFit_percentError", "h_allChannel_ampFit_percentError", 100, -50, 50)

        # *** 2. Add non-automated histograms histograms to array
        self.histArray.AddLast(h_allChannel_timing)
        self.histArray.AddLast(h_allChannel_timingRes)
        self.histArray.AddLast(h_allChannel_timingLogic)
        self.histArray.AddLast(h_allChannel_mcpRef_timingRes)
        self.histArray.AddLast(h_allChannel_x_vs_timingRes)
        self.histArray.AddLast(h_allChannel_x_vs_mcpRef_timingRes)
        self.histArray.AddLast(h_allChannel_fitSlope_vs_mcpRef_timingRes)
        self.histArray.AddLast(h_allChannel_mcpRef_timingRes_ampWalkCorrected)
        self.histArray.AddLast(h_allChannel_mcpRef_ampWalkCorrection)

        self.histArray.AddLast(h_allChannel_ampFit_percentError)
        self.histArray.AddLast(h_allChannel_fracFit_timing)
        self.histArray.AddLast(h_allChannel_fracFit_timingRes)
        self.histArray.AddLast(h_allChannel_mcpRef_fracFit_timingRes)

        # *** 3. Canvases and style
        if self.runBatch:
            ROOT.gROOT.SetBatch(True)
            
        self.c1 = TCanvas("c1", "c1", 800, 800)
        self.c2 = TCanvas("c2", "c2", 800, 800)
        self.c3 = TCanvas("c3", "c3", 800, 800)
        self.c4 = TCanvas("c4", "c4", 800, 800)

        gStyle.SetOptStat(0000)

        # *** 4. Make some directories if not already existent
        if not os.path.isdir(self.topDir):
            os.system( 'mkdir {0}'.format(self.topDir) )
        self.topDir = '{0}/{1}'.format(topDir, runType)
        if not os.path.isdir( '{0}/{1}'.format(topDir, runType) ):
            os.system( 'mkdir {0}'.format(self.topDir) )
  
        # *** 5. Run analysis
        self.loopEvents()

    # =============================

    def addSlicedHistograms(self, arr):
        """ function to add histograms sliced by x"""
        
        # calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        #rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        x = -5
        while x < 35:
            xLow = str(x).replace('-','n')
            xHigh = str(x+self.xSlice).replace('-','n')
            xBase = '{0}_to_{1}'.format(xLow, xHigh)
            #print xBase
        
            h_b1 = TH1D('h_b1_timingRes_{0}'.format(xBase), 'h_b1_timingRes_{0}'.format(xBase), 150, -1500, 1500)
            h_b2 = TH1D('h_b2_timingRes_{0}'.format(xBase), 'h_b2_timingRes_{0}'.format(xBase), 150, -1500, 1500)
            h_b3 = TH1D('h_b3_timingRes_{0}'.format(xBase), 'h_b3_timingRes_{0}'.format(xBase), 150, -1500, 1500)
            h_b4 = TH1D('h_b4_timingRes_{0}'.format(xBase), 'h_b4_timingRes_{0}'.format(xBase), 150, -1500, 1500)
            h_b5 = TH1D('h_b5_timingRes_{0}'.format(xBase), 'h_b5_timingRes_{0}'.format(xBase), 150, -1500, 1500)
            h_b1_mcp = TH1D('h_b1_mcpRef_timingRes_{0}'.format(xBase), 'h_b1_mcpRef_timingRes_{0}'.format(xBase), 175, -3500, 0)
            h_b2_mcp = TH1D('h_b2_mcpRef_timingRes_{0}'.format(xBase), 'h_b2_mcpRef_timingRes_{0}'.format(xBase), 175, -3500, 0)
            h_b3_mcp = TH1D('h_b3_mcpRef_timingRes_{0}'.format(xBase), 'h_b3_mcpRef_timingRes_{0}'.format(xBase), 175, -3500, 0)
            h_b4_mcp = TH1D('h_b4_mcpRef_timingRes_{0}'.format(xBase), 'h_b4_mcpRef_timingRes_{0}'.format(xBase), 175, -3500, 0)
            h_b5_mcp = TH1D('h_b5_mcpRef_timingRes_{0}'.format(xBase), 'h_b5_mcpRef_timingRes_{0}'.format(xBase), 175, -3500, 0)
        
            h_all = TH1D('h_allBar_timingRes_{0}'.format(xBase), 'h_allBar_timingRes_{0}'.format(xBase), 150, -1500, 1500)
            h_all_mcp = TH1D('h_allBar_mcpRef_timingRes_{0}'.format(xBase), 'h_allBar_mcpRef_timingRes_{0}'.format(xBase), 175, -3500, 0)

            arr.AddLast(h_b1)
            arr.AddLast(h_b2)
            arr.AddLast(h_b3)
            arr.AddLast(h_b4)
            arr.AddLast(h_b5)
            arr.AddLast(h_b1_mcp)
            arr.AddLast(h_b2_mcp)
            arr.AddLast(h_b3_mcp)
            arr.AddLast(h_b4_mcp)
            arr.AddLast(h_b5_mcp)

            arr.AddLast(h_all)
            arr.AddLast(h_all_mcp)

            x += self.xSlice


        m = 60
        while m < 150:
            mLow = str(m)
            mHigh = str(m+self.slopeSlice)
            mBase = '{0}_to_{1}'.format(mLow, mHigh)
            #print xBase
         
            h_b1_m = TH1D('h_b1_timingRes_byFitSlope_{0}'.format(mBase), 'h_b1_timingRes_byFitSlope_{0}'.format(mBase), 150, -1500, 1500)
            h_b2_m = TH1D('h_b2_timingRes_byFitSlope_{0}'.format(mBase), 'h_b2_timingRes_byFitSlope_{0}'.format(mBase), 150, -1500, 1500)
            h_b3_m = TH1D('h_b3_timingRes_byFitSlope_{0}'.format(mBase), 'h_b3_timingRes_byFitSlope_{0}'.format(mBase), 150, -1500, 1500)
            h_b4_m = TH1D('h_b4_timingRes_byFitSlope_{0}'.format(mBase), 'h_b4_timingRes_byFitSlope_{0}'.format(mBase), 150, -1500, 1500)
            h_b5_m = TH1D('h_b5_timingRes_byFitSlope_{0}'.format(mBase), 'h_b5_timingRes_byFitSlope_{0}'.format(mBase), 150, -1500, 1500)
            h_b1_mcp_m = TH1D('h_b1_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 'h_b1_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 175, -3500, 0)
            h_b2_mcp_m = TH1D('h_b2_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 'h_b2_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 175, -3500, 0)
            h_b3_mcp_m = TH1D('h_b3_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 'h_b3_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 175, -3500, 0)
            h_b4_mcp_m = TH1D('h_b4_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 'h_b4_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 175, -3500, 0)
            h_b5_mcp_m = TH1D('h_b5_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 'h_b5_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 175, -3500, 0)

            h_all_m = TH1D('h_allBar_timingRes_byFitSlope_{0}'.format(mBase), 'h_allBar_timingRes_byFitSlope_{0}'.format(mBase), 150, -1500, 1500)
            h_all_mcp_m = TH1D('h_allBar_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 'h_allBar_mcpRef_timingRes_byFitSlope_{0}'.format(mBase), 175, -3500, 0)

            arr.AddLast(h_b1_m)
            arr.AddLast(h_b2_m)
            arr.AddLast(h_b3_m)
            arr.AddLast(h_b4_m)
            arr.AddLast(h_b5_m)
            arr.AddLast(h_b1_mcp_m)
            arr.AddLast(h_b2_mcp_m)
            arr.AddLast(h_b3_mcp_m)
            arr.AddLast(h_b4_mcp_m)
            arr.AddLast(h_b5_mcp_m)

            arr.AddLast(h_all_m)
            arr.AddLast(h_all_mcp_m)

            m += self.slopeSlice

        return arr

    # =============================

    def addLeakageHistograms(self, arr):
        """ function to add histograms sliced by x"""
        
        # calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        #rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        barNum = 0
        while barNum < 5:
            barNum +=1
        
            # ** A. Calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
            rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

            # ** B. General plots
            h_b = TH2D('h_b{0}'.format(barNum), 'h_b{0}'.format(barNum), 40, -5, 35, 35, 0, 35)
            h_b_t = TH2D('h_b{0}_t'.format(barNum), 'h_b{0}_t'.format(barNum), 40, -5, 35, 35, 0, 35)
            h_mcp_chR = TH1D('h_mcp{0}_ch{1}'.format(timeChannel, rightSiPMchannel), 'h_mcp{0}_ch{1}'.format(timeChannel, rightSiPMchannel), 25, 0, 500)
            h_chR_vs_chL = TH2D('h_ch{0}_vs_ch{1}'.format(rightSiPMchannel, leftSiPMchannel), 'h_ch{0}_vs_ch{1}'.format(rightSiPMchannel, leftSiPMchannel), 22, 0, 1100, 22, 0, 1100)
            h_chR_chL_x_vs_ratio = TProfile('h_ch{0}_ch{1}_x_vs_ratio'.format(rightSiPMchannel, leftSiPMchannel), 'h_ch{0}_ch{1}_x_vs_ratio'.format(rightSiPMchannel, leftSiPMchannel), 40, -5, 35, 0, 2)

            # ** C. Focused leakage studies
            h_b_left = TH1D('h_trackIn_b{0}_leftSignalInOtherBars'.format(barNum), 'h_trackIn_b{0}_leftSignalInOtherBars'.format(barNum), 25, 0, 100)
            h_b_right = TH1D('h_trackIn_b{0}_rightSignalInOtherBars'.format(barNum), 'h_trackIn_b{0}_rightSignalInOtherBars'.format(barNum), 25, 0, 100)
            h_b_sum = TH1D('h_trackIn_b{0}_sumSignalInOtherBars'.format(barNum), 'h_trackIn_b{0}_sumSignalInOtherBars'.format(barNum), 50, 0, 200)
            h_b_diff = TH1D('h_trackIn_b{0}_diffSignalInOtherBars'.format(barNum), 'h_trackIn_b{0}_diffSignalInOtherBars'.format(barNum), 50, -100, 100)

            h_b_left_out = TH1D('h_trackOut_b{0}_leftSignalInOtherBars'.format(barNum), 'h_trackOut_b{0}_leftSignalInOtherBars'.format(barNum), 25, 0, 100)
            h_b_right_out = TH1D('h_trackOut_b{0}_rightSignalInOtherBars'.format(barNum), 'h_trackOut_b{0}_rightSignalInOtherBars'.format(barNum), 25, 0, 100)
            h_b_sum_out = TH1D('h_trackOut_b{0}_sumSignalInOtherBars'.format(barNum), 'h_trackOut_b{0}_sumSignalInOtherBars'.format(barNum), 50, 0, 200)
            h_b_diff_out = TH1D('h_trackOut_b{0}_diffSignalInOtherBars'.format(barNum), 'h_trackOut_b{0}_diffSignalInOtherBars'.format(barNum), 50, -100, 100)

            arr.AddLast(h_b)
            arr.AddLast(h_b_t)
            arr.AddLast(h_mcp_chR)
            arr.AddLast(h_chR_vs_chL)
            arr.AddLast(h_chR_chL_x_vs_ratio)
            arr.AddLast(h_b_left)
            arr.AddLast(h_b_right)
            arr.AddLast(h_b_sum)
            arr.AddLast(h_b_diff)
            arr.AddLast(h_b_left_out)
            arr.AddLast(h_b_right_out)
            arr.AddLast(h_b_sum_out)
            arr.AddLast(h_b_diff_out)

        return arr

    # =============================

    def addProfiles(self, arr):
        """ function to add histograms sliced by x"""
        
        barNum = 0
        while barNum < 5:
            barNum +=1
            rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

            h_chR_x_vs_amp  = TProfile('h_ch{0}_x_vs_amp'.format(rightSiPMchannel), 'h_ch{0}_x_vs_amp'.format(rightSiPMchannel), 40, -5, 35, 0, 1000)
            h_chR_x_vs_time = TProfile('h_ch{0}_x_vs_time'.format(rightSiPMchannel), 'h_ch{0}_x_vs_time'.format(rightSiPMchannel), 40, -5, 35, 35, 45)
            h_chR_minusMCP_x_vs_time = TProfile('h_ch{0}_minusMCP_x_vs_time'.format(rightSiPMchannel), 'h_ch{0}_minusMCP_x_vs_time'.format(rightSiPMchannel), 40, -5, 35, 35, 45)
            h_chL_x_vs_amp  = TProfile('h_ch{0}_x_vs_amp'.format(leftSiPMchannel), 'h_ch{0}_x_vs_amp'.format(leftSiPMchannel), 40, -5, 35, 0, 1000)
            h_chL_x_vs_time = TProfile('h_ch{0}_x_vs_time'.format(leftSiPMchannel), 'h_ch{0}_x_vs_time'.format(leftSiPMchannel), 40, -5, 35, 35, 45)
            h_chL_minusMCP_x_vs_time = TProfile('h_ch{0}_minusMCP_x_vs_time'.format(leftSiPMchannel), 'h_ch{0}_minusMCP_x_vs_time'.format(leftSiPMchannel), 40, -5, 35, 35, 45)

            h_chR_diffRL_vs_time = TProfile('h_ch{0}_diffRL_vs_time'.format(rightSiPMchannel), 'h_ch{0}_diffRL_vs_time'.format(rightSiPMchannel), 60, -1, 1, 35, 45)
            h_chL_diffRL_vs_time = TProfile('h_ch{0}_diffRL_vs_time'.format(leftSiPMchannel), 'h_ch{0}_diffRL_vs_time'.format(leftSiPMchannel), 60, -1, 1, 35, 45)
            h_chR_chL_over2_diffRL_vs_time = TProfile('h_ch{0}_ch{1}_over2_diffRL_vs_time'.format(rightSiPMchannel, leftSiPMchannel), 'h_ch{0}_ch{1}_over2_diffRL_vs_time'.format(rightSiPMchannel, leftSiPMchannel), 60, -1, 1, 35, 45)
        
            arr.AddLast(h_chR_x_vs_amp)
            arr.AddLast(h_chR_x_vs_time)
            arr.AddLast(h_chR_minusMCP_x_vs_time)
            arr.AddLast(h_chL_x_vs_amp)
            arr.AddLast(h_chL_x_vs_time)
            arr.AddLast(h_chL_minusMCP_x_vs_time)

            arr.AddLast(h_chR_diffRL_vs_time)
            arr.AddLast(h_chL_diffRL_vs_time)
            arr.AddLast(h_chR_chL_over2_diffRL_vs_time)

        return arr

    # =============================

    def setVarsByRunType(self, runType):
        """ set various constants as function of runType"""

        if runType == "all5exposure":
            self.signalThreshold = 800 # should be 800 (using 100 for CERN testbeam comparison)
            self.vetoThreshold   = 100
            self.xBoundaries = [-2, 6, 15, 24, 33]
            self.yBoundaries = [7.5, 12.5, 16.5, 20.5, 24.5]
            self.yIntegralOffset = 2.5
        elif runType == "bottomBars_66V":
            self.signalThreshold = 60
            self.vetoThreshold   = 30
            self.xBoundaries = [17, 19, 21, 23, 25]
            self.yBoundaries = [4.5, 9.5, 13.5, 16.5, 20.5] #5, 9, 13
            self.yIntegralOffset = 2.5

        elif runType == "topBars_66V":
            self.signalThreshold = 60
            self.vetoThreshold   = 30
            self.xBoundaries = [17, 19, 21, 23, 25]
            self.yBoundaries = [25, 25, 2.5, 7.5, 11.5] # 6-9, 1-4 (bins) 10-13
            self.yIntegralOffset = 1.5

    # =============================
    
    def returnVetoDecision(self, event, barNum, vetoOption):
        """ function to return veto decision given option vetoOption = 'None'/'singleAdj'/'doubleAdj'/'allAdj'/'all'"""

        if vetoOption == 'none':
            return False

        # quick veto if >1 track
        if event.ntracks > 1:
            return True
        
        # now go through complete logic
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

    def inBarZone(self, track_x, track_y, barNum):
        """ function to return boolean whether hit inside fill bar-specific plots"""
        
        if (abs(track_y - self.yBoundaries[barNum - 1]) <= self.yIntegralOffset) and (track_x >= self.xBoundaries[0]) and (track_x <= self.xBoundaries[ len(self.xBoundaries)-1 ]):
            return True
        else:
            return False

    # =============================

    def inWhichBar(self, track_x, track_y):
        """ function to return boolean whether hit inside fill bar-specific plots"""
        
        if self.inBarZone(track_x, track_y, 1):
            return 1
        elif self.inBarZone(track_x, track_y, 2):
            return 2
        elif self.inBarZone(track_x, track_y, 3):
            return 3
        elif self.inBarZone(track_x, track_y, 4):
            return 4
        elif self.inBarZone(track_x, track_y, 5):
            return 5

        # if gets here, means hit was not in any bar zone.... weird but okay
        return 0

    # =============================

    def fillChannelPlots(self, event, barNum, arr):
        """ function to fill bar-specific plots"""
        
        # calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        # logic to see if event should be vetoed based on signals in other bars
        doVetoEvent = self.returnVetoDecision(event, barNum, self.vetoOpt) # vetoOpt = none, singleAdj, doubleAdj, allAdj, all

        if (event.amp[rightSiPMchannel] > self.signalThreshold and event.amp[leftSiPMchannel] > self.signalThreshold and not doVetoEvent and abs(event.xSlope) < 0.0004 and abs(event.ySlope) < 0.0004) and event.ntracks == 1:
            # leakage histograms and profiles
            arr.FindObject('h_b{0}'.format(barNum)).Fill(event.x_dut[2], event.y_dut[2])
            arr.FindObject('h_mcp{0}_ch{1}'.format(timeChannel, rightSiPMchannel)).Fill( event.amp[mcpChannel] )
            arr.FindObject('h_ch{0}_vs_ch{1}'.format(rightSiPMchannel, leftSiPMchannel)).Fill( event.amp[rightSiPMchannel], event.amp[leftSiPMchannel] )
            arr.FindObject('h_ch{0}_ch{1}_x_vs_ratio'.format(rightSiPMchannel, leftSiPMchannel)).Fill(event.x_dut[2], event.amp[rightSiPMchannel] / event.amp[leftSiPMchannel] )
            arr.FindObject('h_ch{0}_x_vs_amp'.format(rightSiPMchannel)).Fill(event.x_dut[2], event.amp[rightSiPMchannel] )
            arr.FindObject('h_ch{0}_x_vs_amp'.format(leftSiPMchannel)).Fill(event.x_dut[2], event.amp[leftSiPMchannel] )

            trackInBar = self.inBarZone(event.x_dut[2], event.y_dut[2], barNum)
            # test area for hit integral defintion
            if trackInBar:
                arr.FindObject('h_b{0}_t'.format(barNum)).Fill(event.x_dut[2], event.y_dut[2])

            # fill leakage breakdown histograms
            if trackInBar:
                arr = self.fillLeakageHistograms(arr, event, barNum, 'trackIn')
            else:
                arr = self.fillLeakageHistograms(arr, event, barNum, 'trackOut')
                # fill trace histogram with signal bar waveforms and waveforms from bar with track
                trackBar = self.inWhichBar(event.x_dut[2], event.y_dut[2])
                if trackBar != 0:
                    rightSiPMchannel_track, leftSiPMchannel_track, mcpChannel_track, timeChannel_track = self.returnChannelNumbers(trackBar)
                    if event.i_evt % 10 == 0:
                        self.drawFourChannelTrace(event.time, event.channel, rightSiPMchannel, leftSiPMchannel, timeChannel, rightSiPMchannel_track, leftSiPMchannel_track, timeChannel_track, event.i_evt) #MIXME --> uncomment for leakage studies
                        # if signal and track bar separated by one bar, look at intermediate bar too
                        if abs(trackBar - barNum) == 2:
                            rightSiPMchannel_mid, leftSiPMchannel_mid, mcpChannel_mid, timeChannel_mid = self.returnChannelNumbers( min(trackBar,barNum)+1 )
                            self.drawSixChannelTrace(event.time, event.channel, rightSiPMchannel, leftSiPMchannel, timeChannel, rightSiPMchannel_track, leftSiPMchannel_track, timeChannel_track, rightSiPMchannel_mid, leftSiPMchannel_mid, timeChannel_mid, event.i_evt) #MIXME --> uncomment for leakage studies
                            
            # *** 2. Timing stuff
            # ** A. Break if no timing analysis requested
            if not self.doTiming:
                return arr

            # ** B. Calculate times
            #mipTime_R, fitSlope_R, mipTime_R_ampWalkCorrected = self.getTimingForChannel(event.time, event.channel, timeChannel, rightSiPMchannel, event.i_evt)
            #mipTime_L, fitSlope_L, mipTime_L_ampWalkCorrected = self.getTimingForChannel(event.time, event.channel, timeChannel, leftSiPMchannel, event.i_evt)
            fitStartTime_R, fitStartVoltage_R, fitSlope_R, ampFitPercentErr_R, mipTime_fracFit_R = self.getTimingForChannel(event.time, event.channel, timeChannel, rightSiPMchannel, event.i_evt)
            fitStartTime_L, fitStartVoltage_L, fitSlope_L, ampFitPercentErr_L, mipTime_fracFit_L = self.getTimingForChannel(event.time, event.channel, timeChannel, leftSiPMchannel, event.i_evt)
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

                arr.FindObject('h_ch{0}_diffRL_vs_time'.format(rightSiPMchannel)).Fill((mipTime_R - mipTime_L), mipTime_R )
                arr.FindObject('h_ch{0}_diffRL_vs_time'.format(leftSiPMchannel)).Fill((mipTime_R - mipTime_L), mipTime_L )
                arr.FindObject('h_ch{0}_ch{1}_over2_diffRL_vs_time'.format(rightSiPMchannel, leftSiPMchannel)).Fill((mipTime_R - mipTime_L), (mipTime_L + mipTime_R)/2 )
                
                arr.FindObject('h_allChannel_timingLogic').Fill("Both",1)
                arr.FindObject('h_allChannel_ampFit_percentError').Fill(100*ampFitPercentErr_L)
                arr.FindObject('h_allChannel_ampFit_percentError').Fill(100*ampFitPercentErr_R)

                deltaT         = 1000*(mipTime_L - mipTime_R) # multiple by 1000 to transfer from ns to ps
                deltaT_fracFit = 1000*(mipTime_fracFit_L - mipTime_fracFit_R) # multiple by 1000 to transfer from ns to ps
                arr.FindObject('h_allChannel_x_vs_timingRes').Fill(event.x_dut[2], deltaT)
                #if event.x_dut[2] > 5 and event.x_dut[2] < 25 and event.amp[rightSiPMchannel] > self.fitSignalThreshold: # keep it central
                arr.FindObject('h_allChannel_timingRes').Fill( deltaT ) 
                arr.FindObject('h_allChannel_timing').Fill(mipTime_L)
                arr.FindObject('h_allChannel_timing').Fill(mipTime_R)

                arr.FindObject('h_allChannel_fracFit_timingRes').Fill( deltaT_fracFit ) 
                arr.FindObject('h_allChannel_fracFit_timing').Fill(mipTime_fracFit_L)
                arr.FindObject('h_allChannel_fracFit_timing').Fill(mipTime_fracFit_R)
                # *** fill x-slice plots
                x = -5
                while x < 35:
                    xLow = str(x).replace('-','n')
                    xHigh = str(x+self.xSlice).replace('-','n')
                    xBase = '{0}_to_{1}'.format(xLow, xHigh)
                    if event.x_dut[2] >= x and event.x_dut[2] < x+self.xSlice:
                        #arr.FindObject('h_ch{0}_timingRes_{1}'.format(rightSiPMchannel, xBase)).Fill( deltaT )
                        arr.FindObject('h_b{0}_timingRes_{1}'.format(barNum, xBase)).Fill( deltaT )
                        arr.FindObject('h_allBar_timingRes_{0}'.format(xBase)).Fill( deltaT )
                    x += self.xSlice

                # *** fill slope-slice plots
                m = 60
                while m < 150:
                    mLow = str(m)
                    mHigh = str(m+self.slopeSlice)
                    mBase = '{0}_to_{1}'.format(mLow, mHigh)

                    if fitSlope_R >= m and fitSlope_R < m+self.slopeSlice:
                        arr.FindObject('h_b{0}_timingRes_byFitSlope_{1}'.format(barNum, mBase)).Fill( deltaT )
                        arr.FindObject('h_allBar_timingRes_byFitSlope_{0}'.format(mBase)).Fill( deltaT )
                    if fitSlope_L >= m and fitSlope_L < m+self.slopeSlice:
                        arr.FindObject('h_b{0}_timingRes_byFitSlope_{1}'.format(barNum, mBase)).Fill( deltaT )
                        arr.FindObject('h_allBar_timingRes_byFitSlope_{0}'.format(mBase)).Fill( deltaT )
                    # !!!! slopes are usually not in same slice! investigate later
                    #if fitSlope_R >= m and fitSlope_R < m+self.slopeSlice and fitSlope_L >= m and fitSlope_L < m+self.slopeSlice:
                    #    arr.FindObject('h_allBar_timingRes_byFitSlope_{0}'.format(mBase)).Fill( deltaT )

                    m += self.slopeSlice
                
                # *** Do timing resolution with MCP info ***
                if mipTime_MCP != 0 and event.amp[mcpChannel] > 80 and event.amp[mcpChannel] < 160:
                    deltaT_mcp = 1000*(((mipTime_R + mipTime_L)/2) - mipTime_MCP) # multiple by 1000 to transfer from ns to ps
                    deltaT_mcp_fracFit = 1000*(((mipTime_fracFit_R + mipTime_fracFit_L)/2) - mipTime_MCP) # multiple by 1000 to transfer from ns to ps
                    arr.FindObject('h_allChannel_x_vs_mcpRef_timingRes').Fill(event.x_dut[2], deltaT_mcp)
                    #if event.x_dut[2] > 5 and event.x_dut[2] < 25 and event.amp[rightSiPMchannel] > self.fitSignalThreshold: # keep it central
                    arr.FindObject('h_allChannel_mcpRef_timingRes').Fill( deltaT_mcp )
                    arr.FindObject('h_allChannel_mcpRef_fracFit_timingRes').Fill( deltaT_mcp_fracFit )
                    arr.FindObject('h_allChannel_timing').Fill(mipTime_MCP)
                    #arr.FindObject('h_allChannel_fitSlope_vs_mcpRef_timingRes').Fill(fitSlope_R, deltaT_mcp)
                    #arr.FindObject('h_allChannel_fitSlope_vs_mcpRef_timingRes').Fill(fitSlope_L, deltaT_mcp)
                    arr.FindObject('h_allChannel_fitSlope_vs_mcpRef_timingRes').Fill( (fitSlope_L+fitSlope_R)/2, deltaT_mcp)
                    
                    arr.FindObject('h_ch{0}_minusMCP_x_vs_time'.format(rightSiPMchannel)).Fill(event.x_dut[2], mipTime_R - mipTime_MCP)
                    arr.FindObject('h_ch{0}_minusMCP_x_vs_time'.format(leftSiPMchannel)).Fill(event.x_dut[2], mipTime_L - mipTime_MCP)

                    # *** fill x-slice plots w/ mcp data
                    x = -5
                    while x < 35:
                        xLow = str(x).replace('-','n')
                        xHigh = str(x+self.xSlice).replace('-','n')
                        xBase = '{0}_to_{1}'.format(xLow, xHigh)
                        if event.x_dut[2] >= x and event.x_dut[2] < x+self.xSlice:
                            arr.FindObject('h_b{0}_mcpRef_timingRes_{1}'.format(barNum, xBase)).Fill( deltaT_mcp )
                            arr.FindObject('h_allBar_mcpRef_timingRes_{0}'.format(xBase)).Fill( deltaT_mcp )
                        x += self.xSlice
                    
                    # *** fill slope-slice plots
                    m = 60
                    while m < 150:
                        mLow = str(m)
                        mHigh = str(m+self.slopeSlice)
                        mBase = '{0}_to_{1}'.format(mLow, mHigh)

                        if fitSlope_R >= m and fitSlope_R < m+self.slopeSlice:
                            arr.FindObject('h_b{0}_mcpRef_timingRes_byFitSlope_{1}'.format(barNum, mBase)).Fill( deltaT_mcp )
                            arr.FindObject('h_allBar_mcpRef_timingRes_byFitSlope_{0}'.format(mBase)).Fill( deltaT_mcp )
                        if fitSlope_L >= m and fitSlope_L < m+self.slopeSlice:
                            arr.FindObject('h_b{0}_mcpRef_timingRes_byFitSlope_{1}'.format(barNum, mBase)).Fill( deltaT_mcp )
                            arr.FindObject('h_allBar_mcpRef_timingRes_byFitSlope_{0}'.format(mBase)).Fill( deltaT_mcp )
                        # !!!! slopes are usually not in same slice! investigate later
                        #if fitSlope_R >= m and fitSlope_R < m+self.slopeSlice and fitSlope_L >= m and fitSlope_L < m+self.slopeSlice:
                        #    arr.FindObject('h_allBar_mcpRef_timingRes_byFitSlope_{0}'.format(mBase)).Fill( deltaT_mcp )

                        m += self.slopeSlice
                        
                    # ****   amp-walk corrected plots   ****
                    f_ampWalkCorrection = TF1()
                    if self.vetoOpt == 'singleAdj':
                        f_ampWalkCorrection = TF1("slope_ampWalkCorrected", "(-0.004)*x*x*x + (0.357)*x*x + (-13.681)*x + (-2077.244)") # from 50k run using singleAdj veto
                    if self.vetoOpt == 'doubleAdj' or self.vetoOpt == 'allAdj' or self.vetoOpt == 'all': #FIXME --> should only be doubleAdj but need fix atm
                        #f_ampWalkCorrection = TF1("slope_ampWalkCorrected", "(-2.3285)*x + (-2053.81)") # from 10k run using doubleAdj veto (using R and L independently)
                        #f_ampWalkCorrection = TF1("slope_ampWalkCorrected", "(-5.5129)*x + (-1731.97)") # from 10k run using doubleAdj veto (using R+L/2 )
                        f_ampWalkCorrection = TF1("slope_ampWalkCorrected", "(-4.871)*x + (-1801.585)") # from full run using doubleAdj veto (using R+L/2 )

                    # old approach
                    """mipTime_R_ampWalkCorrected = fitStartTime_R + (self.fitVoltageForTiming - fitStartVoltage_R)/f_ampWalkCorrected.Eval(deltaT_mcp)
                    mipTime_L_ampWalkCorrected = fitStartTime_L + (self.fitVoltageForTiming - fitStartVoltage_L)/f_ampWalkCorrected.Eval(deltaT_mcp)
                    deltaT_mcp_ampWalkCorrected = 1000*(((mipTime_R_ampWalkCorrected + mipTime_L_ampWalkCorrected)/2) - mipTime_MCP) # multiple by 1000 to transfer from ns to ps
                    #print "fitSlope_R: {0}, walkSlope: {1}".format(fitSlope_R, f_ampWalkCorrected.Eval(deltaT_mcp))
                    #print "fitSlope_L: {0}, walkSlope: {1}".format(fitSlope_L, f_ampWalkCorrected.Eval(deltaT_mcp))
                    arr.FindObject('h_allChannel_mcpRef_ampWalkCorrection').Fill( 1000*(mipTime_R - mipTime_R_ampWalkCorrected) )
                    arr.FindObject('h_allChannel_mcpRef_ampWalkCorrection').Fill( 1000*(mipTime_L - mipTime_L_ampWalkCorrected) )
                    """
                    
                    # new approach
                    deltaT_mcp_ampWalkCorrection = f_ampWalkCorrection(110) - f_ampWalkCorrection( (fitSlope_R + fitSlope_L)/2 ) # in ns
                    arr.FindObject('h_allChannel_mcpRef_ampWalkCorrection').Fill( deltaT_mcp_ampWalkCorrection )
                    deltaT_mcp_ampWalkCorrected = deltaT_mcp + deltaT_mcp_ampWalkCorrection # multiple by 1000 to transfer from ns to ps
                    arr.FindObject('h_allChannel_mcpRef_timingRes_ampWalkCorrected').Fill( deltaT_mcp_ampWalkCorrected )
                    #print "deltaT_mcp: {0}, m~: {1}, correction: {2}, deltaT_mcp_corrected: {3}".format(deltaT_mcp, (fitSlope_R + fitSlope_L)/2, deltaT_mcp_ampWalkCorrection, deltaT_mcp_ampWalkCorrected)

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
        peakFit_percentError = 0
        fracTime = 0
        i = 0
        l_time = []
        l_channel = []

        # *** 1. Make arrays of trace information
        while i < 1024:
            l_time.append(time[drs_time*1024 + i])
            l_channel.append(-1*channel[drs_channel*1024 + i])
            i+=1            

        
        # *** 2. Then store a TGraph --> why not in same step? because it doesn't work for unknown reasons
        g, g_peakFit = self.returnWaveformGraph(time, channel, drs_time, drs_channel, vetoClipping=True)
        
        # *** 3. Set up fitting function
        fn1 = TF1("fn1", self.fitFunction) # get function from user config
        startFit = self.getWaveformInfo_TOFPET(l_time, l_channel)
        timeWindow = self.fitTimeWindow
        if drs_channel == 0 or drs_channel == 9:
            startFit = self.getWaveformInfo_MCP(l_time, l_channel)
            timeWindow = self.fitMCPTimeWindow
        

        # *** 4. Perform fit and extract timing if 'good' waveform
        if startFit > 1 and (startFit + timeWindow) < 1024: # protection against weird waveforms
            fn1.SetRange(l_time[startFit], l_time[startFit + timeWindow])
            g.Fit("fn1", "QR")

            # ** A. Get extrapolated peak amplitude and perform fit
            fnPeak = TF1("fnPeak", self.peakFitFunction)
            startPeakFit, endPeakFit = self.getWaveformInfo_peakFit(l_time, l_channel, isLowBias=True)
            fnPeak.SetRange(l_time[startPeakFit], l_time[endPeakFit])
            g_peakFit.Fit("fnPeak", "QR")
            peakAmp = fnPeak.Eval(fnPeak.GetParameter(1))
            peakFit_percentError = ( max(l_channel)-peakAmp ) / max(l_channel) 
            fnR = g.GetFunction("fn1")

            # ** B. Numerical solution for timing info
            fitStop = self.fitVoltageForTiming
            timeStep = l_time[startFit] # timestamp to start at (use startFit cuz... duh)

            if fnR.Eval(timeStep) > self.fitVoltageThreshold: # sometimes we need to push back start when first point of fit is waay beyond fitThresholdVoltage due to steep risetime
                while fnR.Eval(timeStep) > self.fitVoltageThreshold: 
                    timeStep = timeStep - 0.001 
                    if timeStep < 0:
                        break
                        
            evalFit = timeStep

            if drs_channel == 0 or drs_channel == 9:
                fitStop = self.fitMCPVoltageForTiming
                
            fitRes, fitSlope = self.numericalSolve(fnR, evalFit, 0.001, fitStop)
            fitVoltage = fnR.Eval(evalFit)
            if fitSlope != 0:
                timeDiff = (50.000 - fitVoltage)/fitSlope
                slopeRes = evalFit + timeDiff

            # ** C. Dump waveform + fit info for visual inspections    
            if i_evt%500 == 0:
                self.dumpFitInfo(g, fnR, l_time, l_channel, timeStep, drs_channel, i_evt, g_peakFit, fnPeak, fitRes, fitSlope, slopeRes)
                
            if (fitRes - evalFit) < 0.002:
                print "only one step, evt {0}, startFit: ({1:0.3f}, {2:0.3f}), fitted: ({3:0.3f}, {4:0.3f})".format(i_evt, evalFit, fitVoltage, fitRes, fnR.Eval(fitRes))

            # ** D. START CONST FRAC FIT
            if peakAmp != 0: # protection against no peak found
                startFit_constFrac = self.getWaveformInfo_constFrac(l_time, l_channel, peakAmp)
                fnFrac = TF1("fnFrac", self.fitFunction) # get function from user config
                fnFrac.SetRange(l_time[startFit_constFrac], l_time[startFit_constFrac + timeWindow])
                g.Fit("fnFrac", "QR")
                fnFrac = g.GetFunction("fnFrac")
                # * i. Numerical solution for timing info
                fitStop = self.fitPercentForTiming
                timeStep = l_time[startFit_constFrac] # timestamp to start at (use startFit cuz... duh)
                
                # * ii. Get appropriate starting point
                if fnFrac.Eval(timeStep)/peakAmp > self.fitPercentThreshold: # sometimes we need to push back start when first point of fit is waay beyond fitPercentThreshold due to steep risetime
                    while fnFrac.Eval(timeStep)/peakAmp > self.fitPercentThreshold: 
                        timeStep = timeStep - 0.001 
                        if timeStep < 0:
                            break
                evalFracFit = timeStep

                # * iii. Get numerical solution
                fracStop = self.fitPercentForTiming
                fracFitRes, fracFitSlope = self.numericalFracSolve(fnFrac, evalFracFit, 0.001, peakAmp, self.fitPercentForTiming)

        if fitRes <= evalFit or fitRes >= evalFit + 35: # something wonky --> not good
            return 0, 0, 0, 0, 0
        else: # good timing result!
            return evalFit, fitVoltage, fitSlope, peakFit_percentError, fracFitRes

    # =============================

    def dumpFitInfo(self, graph, fit, time, channel, timeStep, drs_channel, i_evt, graphPeak, fitPeak, fitRes, fitSlope, slopeRes):
        """ function to print canvas with lots of fit stuff """

        c5 = TCanvas("c5", "c5", 800, 800)
        p1 = TPad("p1", "p1", 0.0, 0.0, 1.0, 1.0)
        p1.Draw()
        p1.cd()
        graph.Draw()
        
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
            if time[i] > timeStep - 2 and time[i] < timeStep + 0.3:
                inset.SetPoint(n_inset, time[i], channel[i])        
                n_inset += 1
            i=i+1

        inset.Draw()
        #fit.Draw("same")
        fit.DrawF1(timeStep - 2.3, timeStep + 0.2, "same")
        
        l1= TLine(inset.GetXaxis().GetXmin(),0,inset.GetXaxis().GetXmax(),0);
        l1.SetLineStyle(2);
        l1.SetLineWidth(3);
        l1.SetLineColor(600+2);
        l1.Draw("same");
        
        c5.Update()
        c5.Print( "{0}/waveformPlusPol1Fit_Ch{1}_Evt{2}.png".format(self.topDir, drs_channel, i_evt) )
        
        # draw graph without saturated peak
        c5.cd()
        graphPeak.Draw()
        graphPeak.GetYaxis().SetRangeUser(1.4*min(channel), 1.2*fitPeak.Eval(fitPeak.GetParameter(1)))
        fitPeak.Draw("same")
        ltxIn = TLatex()
        ltxIn.SetTextAlign(9)
        ltxIn.SetTextFont(62)
        ltxIn.SetTextSize(0.021)
        ltxIn.SetNDC()
        ltxIn.DrawLatex(0.66, 0.81, "t(max): {0:0.3f} [Fit]".format(fitPeak.GetParameter(1)))
        ltxIn.DrawLatex(0.66, 0.785, "A(max): {0:0.1f} [Fit]".format(fitPeak.Eval(fitPeak.GetParameter(1))))
        ltxIn.DrawLatex(0.66, 0.76, "A(max): {0:0.1f} [Real]".format(max(channel)))
        
        c5.Print( "{0}/waveformWithoutSaturation_Ch{1}_Evt{2}.png".format(self.topDir, drs_channel, i_evt) )
        
        print 'Fit Result = {0}, Fit Slope = {1}'.format(fitRes, fitSlope)
        print 'Lin Result = {0}'.format(slopeRes)
        print 'Lin Amp = {0}, Fit Amp = {1}'.format(fit.Eval(slopeRes), fit.Eval(fitRes))
        print 'Chi2 / NDF = {0:0.1f} / {1:0.1f} = {2:0.1f}'.format(fit.GetChisquare(), fit.GetNDF(), fit.GetChisquare()/fit.GetNDF() ) 
        
        
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

    def getWaveformInfo_constFrac(self, time, channel, maxAmp):
        """ function to calculate and return information about waveform for fitting"""

        # there are probably better ways to do all of this in python...
        threshold_stamp = -1
        i = 0

        # find first reading above voltage threshold i.e. place to start fit
        for reading in channel:
            if abs(reading)/maxAmp > self.fitPercentThreshold and threshold_stamp == -1:
                threshold_stamp = i
            i = i + 1

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

    def getWaveformInfo_peakFit(self, time, channel, isLowBias=False):
        """ function to calculate and return information about waveform for fitting"""

        # there are probably better ways to do all of this in python...
        rise_timestamp = -1
        fall_timestamp = -1
        previousPoint = 0
        i = 0
        
        if isLowBias:
            maxVal = max(channel)
            riseThreshold = 0.4*maxVal
            fallThreshold = 0.85*maxVal
        else:
            riseThreshold = self.peakFitRiseThreshold
            fallThreshold = self.peakFitFallThreshold

        # find first reading above voltage threshold i.e. place to start fit
        for reading in channel:
            if abs(reading) > riseThreshold and rise_timestamp == -1:
                rise_timestamp = i
                fall_timestamp = 0 # start to look for falling
            #if fall_timestamp == 0 and abs(reading) < fallThreshold and abs(reading) < abs(previousPoint):
            if fall_timestamp == 0 and abs(reading) < fallThreshold and time[i] > time[channel.index(max(channel))]:
                fall_timestamp = i
        
            previousPoint = reading

            i = i + 1

        #print 'rise at', rise_timestamp, 'with', channel[rise_timestamp]
        #print 'fall at', fall_timestamp, 'with', channel[fall_timestamp]

        return rise_timestamp, fall_timestamp

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
            if nTotal > 50000 and self.isTest:
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
        self.c2.SetLeftMargin(0.15)
        self.c2.SetRightMargin(0.05)
        self.c2.SetBottomMargin(0.10)
        self.c2.SetTopMargin(0.05)
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

        self.drawTripleProfile(self.c4, 1)
        self.drawTripleProfile(self.c4, 2)
        self.drawTripleProfile(self.c4, 3)
        self.drawTripleProfile(self.c4, 4)
        self.drawTripleProfile(self.c4, 5)

        self.drawBarSplits(self.c2, 'h_trackIn_b_rightSignalInOtherBars')
        self.drawBarSplits(self.c2, 'h_trackIn_b_leftSignalInOtherBars')
        self.drawBarSplits(self.c2, 'h_trackIn_b_sumSignalInOtherBars')
        self.drawBarSplits(self.c2, 'h_trackIn_b_diffSignalInOtherBars')
        self.drawBarSplits(self.c2, 'h_trackOut_b_rightSignalInOtherBars')
        self.drawBarSplits(self.c2, 'h_trackOut_b_leftSignalInOtherBars')
        self.drawBarSplits(self.c2, 'h_trackOut_b_sumSignalInOtherBars')
        self.drawBarSplits(self.c2, 'h_trackOut_b_diffSignalInOtherBars')

        if self.doTiming:
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
            
            self.c4.cd()
            self.histArray.FindObject('h_allChannel_fracFit_timing').Draw()
            self.c4.Print( "{0}/h_allChannel_fracFit_timing.png".format(self.topDir) )
            
            self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_timingRes'), False)
            self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_mcpRef_timingRes'), True)
            self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_mcpRef_timingRes_ampWalkCorrected'), True, isCorrected=True)

            self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_fracFit_timingRes'), False, False, isFracFit=True)
            self.drawResolutionPlot(self.c4, self.histArray.FindObject('h_allChannel_mcpRef_fracFit_timingRes'), True, False, isFracFit=True)

            self.c4.cd()
            self.histArray.FindObject('h_allChannel_ampFit_percentError').SetXTitle("(Fit Amp - Real Amp) / Real Amp")
            self.histArray.FindObject('h_allChannel_ampFit_percentError').SetYTitle("Entries / 1%")
            self.histArray.FindObject('h_allChannel_ampFit_percentError').SetTitle("")
            f_err = TF1("f_err", "gaus", -20, 10) # gaussian
            self.histArray.FindObject('h_allChannel_ampFit_percentError').Fit("f_err", "QR") # should be "R" to impose range
            self.histArray.FindObject('h_allChannel_ampFit_percentError').Draw()
            ltxE = TLatex()
            ltxE.SetTextAlign(9)
            ltxE.SetTextFont(62)
            ltxE.SetTextSize(0.021)
            ltxE.SetNDC()
            ltxE.DrawLatex(0.75, 0.81, "Mean: {0:0.2}".format(f_err.GetParameter(1)))
            ltxE.DrawLatex(0.75, 0.785, "Sigma: {0:0.2}".format(f_err.GetParameter(2)))
            self.c4.Print( "{0}/h_allChannel_ampFit_percentError.png".format(self.topDir) )
            
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

            self.drawTimingResSlices(self.c4, self.histArray, 1, slicedBy='X', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 1, slicedBy='X', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 2, slicedBy='X', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 2, slicedBy='X', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 3, slicedBy='X', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 3, slicedBy='X', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 4, slicedBy='X', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 4, slicedBy='X', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 5, slicedBy='X', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 5, slicedBy='X', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 0, slicedBy='X', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 0, slicedBy='X', usingMCP=True)
            
            self.drawTimingResSlices(self.c4, self.histArray, 1, slicedBy='Slope', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 1, slicedBy='Slope', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 2, slicedBy='Slope', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 2, slicedBy='Slope', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 3, slicedBy='Slope', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 3, slicedBy='Slope', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 4, slicedBy='Slope', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 4, slicedBy='Slope', usingMCP=True)
            self.drawTimingResSlices(self.c4, self.histArray, 5, slicedBy='Slope', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 5, slicedBy='Slope', usingMCP=True)

            self.drawTimingResSlices(self.c4, self.histArray, 0, slicedBy='Slope', usingMCP=False)
            self.drawTimingResSlices(self.c4, self.histArray, 0, slicedBy='Slope', usingMCP=True)

        # === function graveyard. keep for reference"
        #self.drawXquadrants(self.c4, self.histArray.FindObject('h_ch1_ch2_ratio_x1, self.histArray.FindObject('h_ch1_ch2_ratio_x2, self.histArray.FindObject('h_ch1_ch2_ratio_x3, self.histArray.FindObject('h_ch1_ch2_ratio_x4, 1, 'Right/Left')
        #self.drawXquadrants(self.c4, self.histArray.FindObject('h_ch1_x1, self.histArray.FindObject('h_ch1_x2, self.histArray.FindObject('h_ch1_x3, self.histArray.FindObject('h_ch1_x4, 1, 'Bar 1 [Right SiPM]')
        
    # =============================

    def drawResolutionPlot(self, c0, h0, isMCPref, isCorrected=False, isFracFit=False):
        """function to fit resolution plot with gaussian and print resutls"""
        xTitle = 't_{right SiPM} - t_{left SiPM} [ps]'
        plotName = 'h_allChannel_timingRes'
        xMax = 50
        xMin = -450
        if isMCPref:
            xTitle = '(t_{right SiPM} + t_{left SiPM})/2 - t_{MCP} [ps]'
            plotName = 'h_allChannel_mcpRef_timingRes'
            xMax = -2100
            xMin = -2500
            if isCorrected:
                xTitle = 'Amp-Walk Corrected (t_{right SiPM} + t_{left SiPM})/2 - t_{MCP} [ps]'
                plotName = 'h_allChannel_mcpRef_timingRes_ampWalkCorrected'
                xMax = -2100
                xMin = -2500

        if isFracFit:
            plotName = plotName.replace('timingRes', 'fracFit_timingRes')
            xTitle = 'Frac Fit ' + xTitle
            if isMCPref:
                xMax = -1500
                xMin = -2000

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
                #h_p.SetMinimum(60)
                #h_p.SetMaximum(180)
                #xMin = -2800 
                #xMax = -1800 
                h_p.SetMinimum(-2500)
                h_p.SetMaximum(-2000)
                xMin = 60
                xMax = 180 
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
                h_p1.SetXTitle('Fit Slope')
                h_p1.SetYTitle('{0} [ps]'.format(varName))
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
                
                if self.runType == "all5exposure":
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

    def drawTripleProfile(self, c0, barNum, opt='x'):
        """ function to recieve canvas (c0), TProfile (h_p), and bar number"""
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.10);
        
        rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        pR = self.histArray.FindObject('h_ch{0}_diffRL_vs_time'.format(rightSiPMchannel))
        pL = self.histArray.FindObject('h_ch{0}_diffRL_vs_time'.format(leftSiPMchannel))
        pAvg = self.histArray.FindObject('h_ch{0}_ch{1}_over2_diffRL_vs_time'.format(rightSiPMchannel, leftSiPMchannel))

        
        ytitle = '{0} [ns]'.format("#Delta t")
        filename = '{0}/bar{1}_tripleProfile_timing_vs_x.png'.format(self.topDir, barNum)
        pR.SetXTitle("t_{R} - t_{L} [ns]")
        pR.SetYTitle(ytitle)
        pR.SetTitle("Bar {0}".format(barNum))

        pAvg.SetLineColor(1) # kBlack
        pR.SetLineColor(600) # kBlue
        pL.SetLineColor(632) # kRed

        pAvg.SetLineWidth(3)
        pR.SetLineWidth(3)
        pL.SetLineWidth(3)

        pR.SetMinimum(37)
        pR.SetMaximum(41)
        pR.Draw()
        pL.Draw('same')
        pAvg.Draw('same')

        leg = TLegend(0.38, 0.68, .53, .88);
        leg.AddEntry(pR, "t_{R}", "l");
        leg.AddEntry(pL, "t_{L}", "l");
        leg.AddEntry(pAvg, "t_{R}+t_{L}/2", "l");
        leg.Draw("same")

        c0.Print(filename)

    # =============================

    def drawTimingResSlices(self, c0, arr, barNum, slicedBy, usingMCP=False):
        """ function to make pretty canvas of x-slices of timing res """
        
        # ** 0. Calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        rightSiPMchannel, leftSiPMchannel, mcpChannel, timeChannel = self.returnChannelNumbers(barNum)

        # ** 1. Histogram naming stuff
        hname = 'h_b{0}_timingRes'.format(barNum)
        if usingMCP:
            hname = 'h_b{0}_mcpRef_timingRes'.format(barNum)
        if slicedBy == 'Slope':
            hname += '_byFitSlope'
        if barNum == 0:
            hname = hname.replace('b0', 'allBar')
            
        leg = TLegend(0.68, 0.68, .93, .88);
        legendCol = [1, 600, 632, 416+2, 616-3, 800+7, 432-3, 900-3] # kBlack, kBlue, kRed, kGreen+2, kMagenta-3, kOrange+7, kCyan-3, kPink-3
        nSlices = 0
        maxVal = 0

        # ** 2. Set things depending on slicedBy
        s = 0
        sMin = 0
        sMax = 0
        sStep = 0

        if slicedBy=='X':
            s = -5
            sMax = 35
            sStep = self.xSlice
        if slicedBy=='Slope':
            s = 60
            sMax = 150
            sStep = self.slopeSlice

        # ** 3. Start loop
        c0.cd()
        while s < sMax:
            sLow = str(s).replace('-','n')
            sHigh = str(s+sStep).replace('-','n')
            sBase = '{0}_to_{1}'.format(sLow, sHigh)

            #h_b = TH1D('{0}_{1}'.format(hname, xBase), '{0}_{1}'.format(hname, xBase), 300, -1500, 1500)
            h_b = arr.FindObject('{0}_{1}'.format(hname, sBase))
            h_b.SetLineColor(legendCol[nSlices])
            leg.AddEntry(h_b, sBase.replace('n','-').replace('_',' '), "l")
            if nSlices == 0:
                maxVal = h_b.GetBinContent(h_b.GetMaximumBin())
                h_b.SetXTitle('t_{right SiPM} - t_{left SiPM} [ps]')
                if usingMCP:
                    h_b.SetXTitle('(t_{right SiPM} + t_{left SiPM})/2 - t_{MCP} [ps]')
                h_b.GetYaxis().SetRangeUser(0, 1.4*maxVal)
                h_b.SetYTitle("Normalized Entries / 20 ps")
                h_b.SetTitle("")
                h_b.DrawNormalized()
            else:
                h_b.DrawNormalized("same")

            nSlices += 1
            s += sStep
            
        leg.Draw("same")
        hname = hname.strip('_byFitSlope')
        c0.Print('{0}/{1}_by{2}Slice.png'.format(self.topDir, hname, slicedBy))

    # =============================

    def numericalSolve(self, func, x_eval, x_step, y_end):
        """ function to recieve TF1, find + return numerical solution by incrementing x_eval by x_step to point where function equals y_end"""
        y_calc = 0
        x_start = x_eval

        while abs(y_calc) < abs(y_end):
            store = y_calc
            y_calc = func.Eval(x_eval)
            x_eval += x_step
            if x_eval > x_start + 35: # scan over window of 35 ps
                break
                
                
        f1 = func.Eval(x_start)
        f2 = func.Eval(x_eval)
        fitSlope = (f2 - f1)/(x_eval - x_start)

        return x_eval, fitSlope

    # =============================

    def numericalFracSolve(self, func, x_eval, x_step, fracDenom, y_end):
        """ function to recieve TF1, find + return numerical solution by incrementing x_eval by x_step to point where function equals y_end"""
        y_calc = 0
        x_start = x_eval

        while abs(y_calc) < abs(y_end):
            store = y_calc
            y_calc = func.Eval(x_eval)/fracDenom
            x_eval += x_step
            if x_eval > x_start + 35: # scan over window of 35 ps
                break
                
                
        f1 = func.Eval(x_start)
        f2 = func.Eval(x_eval)
        fitSlope = (f2 - f1)/(x_eval - x_start)

        return x_eval, fitSlope
    
    # =============================

    def drawTwoChannelTrace(self, time, channel, ch1, ch2, chTime, i_evt, altColor=False):
        """ function to draw 2-channel traces --> probably just for leakage studies"""

        c0 = TCanvas("c0", "c0", 800, 800)
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.05);
        
        g1 = self.returnWaveformGraph(time, channel, chTime, ch1)
        g2 = self.returnWaveformGraph(time, channel, chTime, ch2)
        g1.SetLineColor(600) #kBlue
        g2.SetLineColor(632) #kRed

        c0.cd()
        g1.Draw()
        g2.Draw("same")
        c0.Print( "{0}/trackOutWaveforms/waveformTrackOut_ch{1}_ch{2}_Evt{3}.png".format(self.topDir, ch1, ch2, i_evt) )

    # =============================

    def drawFourChannelTrace(self, time, channel, ch1, ch2, chTime, chA, chB, chTime2, i_evt):
        """ function to draw 2-channel traces --> probably just for leakage studies"""

        c0 = TCanvas("c0", "c0", 800, 800)
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.05);
        
        g1 = self.returnWaveformGraph(time, channel, chTime, ch1)
        g2 = self.returnWaveformGraph(time, channel, chTime, ch2)
        gA = self.returnWaveformGraph(time, channel, chTime2, chA)
        gB = self.returnWaveformGraph(time, channel, chTime2, chB)
        g1.SetLineColor(1) #kBlack
        g2.SetLineColor(416+2) #kGreen+2
        gA.SetLineColor(600) #kBlue
        gB.SetLineColor(632) #kRed

        c0.cd()
        g1.Draw()
        g2.Draw("same")
        gA.Draw("same")
        gB.Draw("same")
        c0.Print( "{0}/trackOutWaveforms/waveformTrackOut_4channel_sig_ch{1}_ch{2}_track_ch{3}_ch{4}_Evt{5}.png".format(self.topDir, ch1, ch2, chA, chB, i_evt) )

    # =============================

    def drawSixChannelTrace(self, time, channel, ch1, ch2, chTime, chA, chB, chTime2, chI, chII, chTime3, i_evt):
        """ function to draw 2-channel traces --> probably just for leakage studies"""

        c0 = TCanvas("c0", "c0", 800, 800)
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.05);
        
        g1  = self.returnWaveformGraph(time, channel, chTime, ch1)
        g2  = self.returnWaveformGraph(time, channel, chTime, ch2)
        gA  = self.returnWaveformGraph(time, channel, chTime2, chA)
        gB  = self.returnWaveformGraph(time, channel, chTime2, chB)
        gI  = self.returnWaveformGraph(time, channel, chTime3, chI)
        gII = self.returnWaveformGraph(time, channel, chTime3, chII)
        g1.SetLineColor(1) #kBlack
        g2.SetLineColor(416+2) #kGreen+2
        gA.SetLineColor(600) #kBlue
        gB.SetLineColor(632) #kRed
        gI.SetLineColor(616-3) #kMagenta-3
        gII.SetLineColor(432-3) #kCyan-3
        #legendCol = [1, 600, 632, 416+2, 616-3, 800+7, 432-3, 900-3] # kBlack, kBlue, kRed, kGreen+2, kMagenta-3, kOrange+7, kCyan-3, kPink-3

        c0.cd()
        g1.Draw()
        g2.Draw("same")
        gA.Draw("same")
        gB.Draw("same")
        gI.Draw("same")
        gII.Draw("same")
        c0.Print( "{0}/trackOutWaveforms/waveformTrackOut_6channel_sig_ch{1}_ch{2}_track_ch{3}_ch{4}_middle_ch{5}_ch{6}_Evt{7}.png".format(self.topDir, ch1, ch2, chA, chB, chI, chII, i_evt) )

    # =============================

    def fillLeakageHistograms(self, arr, event, barNum, trackIn):
        """ function to fill histograms breaking down information about leakage"""
        
        if not os.path.isdir( '{0}/trackOutWaveforms'.format(self.topDir) ):
            os.system( 'mkdir {0}/trackOutWaveforms'.format(self.topDir) )
  
        if barNum != 1:
            arr.FindObject('h_{0}_b{1}_rightSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[1])
            arr.FindObject('h_{0}_b{1}_leftSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[2])
            arr.FindObject('h_{0}_b{1}_sumSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[1] + event.amp[2])
            arr.FindObject('h_{0}_b{1}_diffSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[1] - event.amp[2])
            #if event.i_evt % 500 == 0:
            #self.drawTwoChannelTrace(event.time, event.channel, 1, 2, 0, event.i_evt)
        if barNum != 2:
            arr.FindObject('h_{0}_b{1}_rightSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[3])
            arr.FindObject('h_{0}_b{1}_leftSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[4])
            arr.FindObject('h_{0}_b{1}_sumSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[3] + event.amp[4])
            arr.FindObject('h_{0}_b{1}_diffSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[3] - event.amp[4])
            #if event.i_evt % 500 == 0:
            #self.drawTwoChannelTrace(event.time, event.channel, 3, 4, 0, event.i_evt)
        if barNum != 3:
            arr.FindObject('h_{0}_b{1}_rightSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[5])
            arr.FindObject('h_{0}_b{1}_leftSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[6])
            arr.FindObject('h_{0}_b{1}_sumSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[5] + event.amp[6])
            arr.FindObject('h_{0}_b{1}_diffSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[5] - event.amp[6])
            #if event.i_evt % 500 == 0:
            #self.drawTwoChannelTrace(event.time, event.channel, 5, 6, 0, event.i_evt)
        if barNum != 4:
            arr.FindObject('h_{0}_b{1}_rightSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[10])
            arr.FindObject('h_{0}_b{1}_leftSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[11])
            arr.FindObject('h_{0}_b{1}_sumSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[10] + event.amp[11])
            arr.FindObject('h_{0}_b{1}_diffSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[10] - event.amp[11])
            #if event.i_evt % 500 == 0:
            #self.drawTwoChannelTrace(event.time, event.channel, 10, 11, 1, event.i_evt)
        if barNum != 5:
            arr.FindObject('h_{0}_b{1}_rightSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[12])
            arr.FindObject('h_{0}_b{1}_leftSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[13])
            arr.FindObject('h_{0}_b{1}_sumSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[12] + event.amp[13])
            arr.FindObject('h_{0}_b{1}_diffSignalInOtherBars'.format(trackIn, barNum)).Fill(event.amp[12] - event.amp[13])
            #if event.i_evt % 500 == 0:
            #self.drawTwoChannelTrace(event.time, event.channel, 12, 13, 1, event.i_evt)

        return arr
        
    # =============================

    def drawBarSplits(self, c0, hname):
        """ function to produce multi-color histogram split by bar """
        
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.05);

        h1 = self.histArray.FindObject(hname.replace('_b_', '_b1_'))
        h2 = self.histArray.FindObject(hname.replace('_b_', '_b2_'))        
        h3 = self.histArray.FindObject(hname.replace('_b_', '_b3_'))
        h4 = self.histArray.FindObject(hname.replace('_b_', '_b4_'))        
        h5 = self.histArray.FindObject(hname.replace('_b_', '_b5_'))
        
        # kBlack == 1, kRed == 632, kBlue == 600, kGreen == 416, kMagenta == 616
        h1.SetLineColor(1) # kBlack
        h2.SetLineColor(600) # kBlue
        h3.SetLineColor(632) # kRed
        h4.SetLineColor(416+2) # kGreen+2
        h5.SetLineColor(616-3) # kMagenta-3
        
        h1.SetLineWidth(3) 
        h2.SetLineWidth(3) 
        h3.SetLineWidth(3) 
        h4.SetLineWidth(3) 
        h5.SetLineWidth(3) 
        
        h1.SetYTitle("Noramlized Entries / 4 mV")
        h1.SetXTitle("Leakage Amplitude [mV]")
        h1.SetTitle("")
        h1.DrawNormalized()
        h1.GetYaxis().SetRangeUser(0,2.6)
        
        h2.DrawNormalized("same")
        h3.DrawNormalized("same")
        h4.DrawNormalized("same")
        h5.DrawNormalized("same")
        
        leg = TLegend(0.58, 0.63, .93, .93);
        leg.AddEntry(h1, "Bar 1 Signal", "l");
        leg.AddEntry(h2, "Bar 2 Signal", "l");
        leg.AddEntry(h3, "Bar 3 Signal", "l");
        leg.AddEntry(h4, "Bar 4 Signal", "l");
        leg.AddEntry(h5, "Bar 5 Signal", "l");
        leg.Draw("same");
        
        c0.Print("{0}/{1}.png".format(self.topDir, hname.replace('_b_', '_byBar_') ))
    
    # =============================

    def returnWaveformGraph(self, time, channel, drs_time, drs_channel, vetoClipping=False):
        """ function to produce and return graph of waveform """

        l_time = []
        l_channel = []
        i=0

        # *** 1. Make arrays of trace information
        while i < 1024:
            l_time.append(time[drs_time*1024 + i])
            l_channel.append(-1*channel[drs_channel*1024 + i])
            i+=1            
        
        # *** 2. Then store a TGraph --> why not in same step? because it doesn't work for unknown reasons
        i=0
        g = TGraph()       
        i_peakFit=0
        g_peakFit = TGraph()       
        voltageVeto = self.peakFitVoltageVeto
        if '66V' in self.runType:
            voltageVeto = 0.875*max(l_channel)

        while i < 1024:
            g.SetPoint(i, l_time[i], l_channel[i])        
            if l_channel[i] < voltageVeto:
                g_peakFit.SetPoint(i_peakFit, l_time[i], l_channel[i])        
                i_peakFit += 1
            i=i+1

        
        if not vetoClipping:
            return g
        else:
            return g, g_peakFit
        
    # =============================
    
