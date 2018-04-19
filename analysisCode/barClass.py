# !/usr/bin/python

#Author: Ben Tannenwald
#Date: April 16, 2018
#Purpose: Class for handling testbeam bar data

import os,sys, argparse
from ROOT import gROOT, TH1D, TFile, TTree, TChain, TCanvas, TH2D, TLegend, gStyle, TLatex, TProfile, TF1, TGraph

class barClass:
    def __init__(self, tree, runType, topDir, vetoOpt, test=False):
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
        self.fitVoltageThreshold = 200 # in mV
        self.fitVoltageForTiming = 250 # in mV
        self.fitTimeWindow = 4
        self.fitFunction = "pol1"

        self.h_b1 = TH2D("h_b1", "h_b1", 40, -5, 35, 35, 0, 35)
        self.h_b2 = TH2D("h_b2", "h_b2", 40, -5, 35, 35, 0, 35)
        self.h_b3 = TH2D("h_b3", "h_b3", 40, -5, 35, 35, 0, 35)
        self.h_b4 = TH2D("h_b4", "h_b4", 40, -5, 35, 35, 0, 35)
        self.h_b5 = TH2D("h_b5", "h_b5", 40, -5, 35, 35, 0, 35)
        self.h_b1_t = TH2D("h_b1_t", "h_b1_t", 40, -5, 35, 35, 0, 35)
        self.h_b2_t = TH2D("h_b2_t", "h_b2_t", 40, -5, 35, 35, 0, 35)
        self.h_b3_t = TH2D("h_b3_t", "h_b3_t", 40, -5, 35, 35, 0, 35)
        self.h_b4_t = TH2D("h_b4_t", "h_b4_t", 40, -5, 35, 35, 0, 35)
        self.h_b5_t = TH2D("h_b5_t", "h_b5_t", 40, -5, 35, 35, 0, 35)
        
        self.h_mcp0_ch1 = TH1D("h_mcp0_bar1", "h_mcp0_ch1", 25, 0, 500);
        self.h_mcp0_ch3 = TH1D("h_mcp0_bar3", "h_mcp0_ch3", 25, 0, 500);
        self.h_mcp0_ch5 = TH1D("h_mcp0_bar5", "h_mcp0_ch5", 25, 0, 500);
        self.h_mcp1_ch10 = TH1D("h_mcp1_ch10", "h_mcp1_ch10", 25, 0, 500);
        self.h_mcp1_ch12 = TH1D("h_mcp1_ch12", "h_mcp1_ch12", 25, 0, 500);
        
        self.h_ch1_vs_ch2 = TH2D("h_ch1_vs_ch2", "h_ch1_vs_ch2", 22, 0, 1100, 22, 0, 1100)
        self.h_ch3_vs_ch4 = TH2D("h_ch3_vs_ch4", "h_ch3_vs_ch4", 22, 0, 1100, 22, 0, 1100)
        self.h_ch5_vs_ch6 = TH2D("h_ch5_vs_ch6", "h_ch5_vs_ch6", 22, 0, 1100, 22, 0, 1100)
        self.h_ch10_vs_ch11 = TH2D("h_ch10_vs_ch11", "h_ch10_vs_ch11", 22, 0, 1100, 22, 0, 1100)
        self.h_ch12_vs_ch13 = TH2D("h_ch12_vs_ch13", "h_ch12_vs_ch13", 22, 0, 1100, 22, 0, 1100)
        
        self.h_ch1_ch2_x_vs_ratio = TProfile("h_ch1_ch2_x_vs_ratio", "h_ch1_ch2_x_vs_ratio", 40, -5, 35, 0, 2)
        self.h_ch3_ch4_x_vs_ratio = TProfile("h_ch3_ch4_x_vs_ratio", "h_ch3_ch4_x_vs_ratio", 40, -5, 35, 0, 2)
        self.h_ch5_ch6_x_vs_ratio = TProfile("h_ch5_ch6_x_vs_ratio", "h_ch5_ch6_x_vs_ratio", 40, -5, 35, 0, 2)
        self.h_ch10_ch11_x_vs_ratio = TProfile("h_ch10_ch11_x_vs_ratio", "h_ch10_ch11_x_vs_ratio", 40, -5, 35, 0, 2)
        self.h_ch12_ch13_x_vs_ratio = TProfile("h_ch12_ch13_x_vs_ratio", "h_ch12_ch13_x_vs_ratio", 40, -5, 35, 0, 2)

        self.h_ch1_x_vs_amp = TProfile("h_ch1_x_vs_amp", "h_ch1_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch2_x_vs_amp = TProfile("h_ch2_x_vs_amp", "h_ch2_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch3_x_vs_amp = TProfile("h_ch3_x_vs_amp", "h_ch3_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch4_x_vs_amp = TProfile("h_ch4_x_vs_amp", "h_ch4_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch5_x_vs_amp = TProfile("h_ch5_x_vs_amp", "h_ch5_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch6_x_vs_amp = TProfile("h_ch6_x_vs_amp", "h_ch6_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch10_x_vs_amp = TProfile("h_ch10_x_vs_amp", "h_ch10_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch11_x_vs_amp = TProfile("h_ch11_x_vs_amp", "h_ch11_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch12_x_vs_amp = TProfile("h_ch12_x_vs_amp", "h_ch12_x_vs_amp", 40, -5, 35, 0, 1000)
        self.h_ch13_x_vs_amp = TProfile("h_ch13_x_vs_amp", "h_ch13_x_vs_amp", 40, -5, 35, 0, 1000)

        self.h_ch1_x_vs_time = TProfile("h_ch1_x_vs_time", "h_ch1_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch2_x_vs_time = TProfile("h_ch2_x_vs_time", "h_ch2_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch3_x_vs_time = TProfile("h_ch3_x_vs_time", "h_ch3_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch4_x_vs_time = TProfile("h_ch4_x_vs_time", "h_ch4_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch5_x_vs_time = TProfile("h_ch5_x_vs_time", "h_ch5_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch6_x_vs_time = TProfile("h_ch6_x_vs_time", "h_ch6_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch10_x_vs_time = TProfile("h_ch10_x_vs_time", "h_ch10_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch11_x_vs_time = TProfile("h_ch11_x_vs_time", "h_ch11_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch12_x_vs_time = TProfile("h_ch12_x_vs_time", "h_ch12_x_vs_time", 40, -5, 35, 0, 100)
        self.h_ch13_x_vs_time = TProfile("h_ch13_x_vs_time", "h_ch13_x_vs_time", 40, -5, 35, 0, 100)

        self.h_allChannel_timing = TH1D("h_allChannel_timing", "h_allChannel_timing", 50, 0, 50)
        self.h_allChannel_timingRes = TH1D("h_allChannel_timingRes", "h_allChannel_timingRes", 20, -2000, 2000)

        self.c1 = TCanvas("c1", "c1", 800, 800)
        self.c2 = TCanvas("c2", "c2", 800, 800)
        self.c3 = TCanvas("c3", "c3", 800, 800)
        self.c4 = TCanvas("c4", "c4", 800, 800)

        gStyle.SetOptStat(0000)

        # make some directories if not already existent
        if not os.path.isdir(self.topDir):
            os.system( 'mkdir {0}'.format(self.topDir) )
        self.topDir = '{0}/{1}'.format(topDir, runType)
        if not os.path.isdir( '{0}/{1}'.format(topDir, runType) ):
            os.system( 'mkdir {0}'.format(self.topDir) )

        # run analysis
        self.loopEvents()

    # =============================

    def setVarsByRunType(self, runType):
        """ set various constants as function of runType"""

        if runType == "all5exposure":
            self.signalThreshold = 100
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

    def fillChannelPlots(self, event, barNum, h_b, h_mcp, h_r_vs_l, h_lr_x_vs_ratio, h_b_test, h_r_x_vs_amp, h_l_x_vs_amp, h_r_x_vs_time, h_l_x_vs_time, h_all_time, h_all_timeRes):
        """ function to fill bar-specific plots"""
        
        # calculate channel numbers given bar number --> there is probably a smarter way to automate this with fewer lines
        leftSiPMchannel = rightSiPMchannel = mcpChannel = timeChannel = -1
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
            
        # logic to see if event should be vetoed based on signals in other bars
        doVetoEvent = self.returnVetoDecision(event, barNum, self.vetoOpt) # vetoOpt = none, singleAdj, doubleAdj, allAdj, all

        if (event.amp[rightSiPMchannel] > self.signalThreshold and not doVetoEvent and abs(event.xSlope) < 0.0004 and abs(event.ySlope) < 0.0004):
            # leakage histograms and profiles
            h_b.Fill(event.x_dut[2], event.y_dut[2])
            h_mcp.Fill( event.amp[mcpChannel] )
            h_r_vs_l.Fill( event.amp[rightSiPMchannel], event.amp[leftSiPMchannel] )
            h_lr_x_vs_ratio.Fill(event.x_dut[2], event.amp[rightSiPMchannel] / event.amp[leftSiPMchannel] )
            h_r_x_vs_amp.Fill(event.x_dut[2], event.amp[rightSiPMchannel] )
            h_l_x_vs_amp.Fill(event.x_dut[2], event.amp[leftSiPMchannel] )

            # test area for hit integral defintion
            if (abs(event.y_dut[2] - self.yBoundaries[barNum - 1]) <= self.yIntegralOffset and event.x_dut[2]>=self.xBoundaries[0] and event.x_dut[2]<=self.xBoundaries[ len(self.xBoundaries)-1 ]):
                h_b_test.Fill(event.x_dut[2], event.y_dut[2])
                    
            # timing stuff
            #print len(event.time), len(event.channel)
            mipTime_R = self.getTimingForChannel(event.time, event.channel, timeChannel, rightSiPMchannel, event.i_evt)
            mipTime_L = self.getTimingForChannel(event.time, event.channel, timeChannel, leftSiPMchannel, event.i_evt)

            if mipTime_R != 0 and mipTime_L != 0:
                h_r_x_vs_time.Fill(event.x_dut[2], mipTime_R)
                h_l_x_vs_time.Fill(event.x_dut[2], mipTime_L)
                h_all_time.Fill(mipTime_L)
                h_all_time.Fill(mipTime_R)
                h_all_timeRes.Fill( 1000*(mipTime_L - mipTime_R) ) # multiple by 1000 to transfer from ns to ps

        return h_b, h_mcp, h_r_vs_l, h_lr_x_vs_ratio, h_b_test, h_r_x_vs_amp, h_l_x_vs_amp, h_r_x_vs_time, h_l_x_vs_time, h_all_time, h_all_timeRes

    # =============================

    def getTimingForChannel(self, time, channel, drs_time, drs_channel, i_evt):
        """ function to calculate and return information about waveform for fitting"""
        voltageFromFunction = 0
        timeStep = 29.999 # timestamp to start at --> rough guess from looking at waveforms
        
        i = 0
        l_time = []
        l_channel = []
        # *** 1. Make arrays of trace information
        while i < 1024:
            l_time.append(time[drs_time*1024 + i])
            l_channel.append(channel[drs_channel*1024 + i])
            i+=1            

        # *** 2. Then store a TGraph --> why not in same step? because it doesn't work for unknown reasons
        i=0
        g = TGraph()       
        while i < 1024:
            g.SetPoint(i, l_time[i], l_channel[i])        
            i=i+1

        fn1 = TF1("fn1", self.fitFunction) # first degree polynomial --> this is just a choice atm
        #fn1 = TF1('fn1', 'gaus(0)') # first degree polynomial --> this is just a choice atm
        startFit = self.getWaveformInfo_TOFPET(l_time, l_channel)
        if startFit > 1 and (startFit + self.fitTimeWindow) < 1024: # protection against weird waveforms
            fn1.SetRange(l_time[startFit], l_time[startFit + self.fitTimeWindow])
            #print 'N_points={0}, fit point {1}: ({2}, {3}), fit point {4}: ({5}, {6})'.format(g.GetN(), startFit, l_time[startFit], l_channel[startFit], startFit+ self.fitTimeWindow, l_time[startFit + self.fitTimeWindow], l_channel[startFit + self.fitTimeWindow])
            g.Fit("fn1", "QR")

            fnR = g.GetFunction("benFnc1")
            if i_evt%5000 == 0:
                c5 = TCanvas("c5", "c5", 800, 800)
                g.Draw()
                c5.Print( "{0}/waveformPlusPol1Fit_Ch{1}_Evt{2}.png".format(self.topDir, drs_channel, i_evt) )
            
            # numerical solution for timing info
            while abs(voltageFromFunction) < abs(self.fitVoltageForTiming):
                store = voltageFromFunction
                voltageFromFunction = fn1.Eval(timeStep)
                #print 'old: {0}, new: {1}'.format(store, voltageFromFunction)
                timeStep += 0.001
                if timeStep > 50:
                    break
               
        
        if timeStep <= 30 or timeStep > 50: # something wonky --> not good
            return 0
        else: # good timing result!
            #print timeStep
            return timeStep
        
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
    
    def loopEvents(self):
        """ function looping over all events in file"""

        print self.tree.GetEntries()
        nTotal=0

        for event in self.tree:        
            if nTotal > 10000 and self.isTest:
                break

            nTotal += 1
            if (nTotal % 10000 == 0):
                print nTotal, "processed"

            # bar 1
            self.h_b1, self.h_mcp0_ch1, self.h_ch1_vs_ch2, self.h_ch1_ch2_x_vs_ratio, self.h_b1_t, self.h_ch1_x_vs_amp, self.h_ch2_x_vs_amp, self.h_ch1_x_vs_time, self.h_ch2_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes = self.fillChannelPlots(event, 1, self.h_b1, self.h_mcp0_ch1, self.h_ch1_vs_ch2, self.h_ch1_ch2_x_vs_ratio, self.h_b1_t, self.h_ch1_x_vs_amp, self.h_ch2_x_vs_amp, self.h_ch1_x_vs_time, self.h_ch2_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes)
            # bar 2
            self.h_b2, self.h_mcp0_ch3, self.h_ch3_vs_ch4, self.h_ch3_ch4_x_vs_ratio, self.h_b2_t, self.h_ch3_x_vs_amp, self.h_ch4_x_vs_amp, self.h_ch3_x_vs_time, self.h_ch4_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes = self.fillChannelPlots(event, 2, self.h_b2, self.h_mcp0_ch3, self.h_ch3_vs_ch4, self.h_ch3_ch4_x_vs_ratio, self.h_b2_t, self.h_ch3_x_vs_amp, self.h_ch4_x_vs_amp, self.h_ch3_x_vs_time, self.h_ch4_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes)
            # bar 3
            self.h_b3, self.h_mcp0_ch5, self.h_ch5_vs_ch6, self.h_ch5_ch6_x_vs_ratio, self.h_b3_t, self.h_ch5_x_vs_amp, self.h_ch6_x_vs_amp, self.h_ch5_x_vs_time, self.h_ch6_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes = self.fillChannelPlots(event, 3, self.h_b3, self.h_mcp0_ch5, self.h_ch5_vs_ch6, self.h_ch5_ch6_x_vs_ratio, self.h_b3_t, self.h_ch5_x_vs_amp, self.h_ch6_x_vs_amp, self.h_ch5_x_vs_time, self.h_ch6_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes)
            # bar 4
            self.h_b4, self.h_mcp1_ch10, self.h_ch10_vs_ch11, self.h_ch10_ch11_x_vs_ratio, self.h_b4_t, self.h_ch10_x_vs_amp, self.h_ch11_x_vs_amp, self.h_ch10_x_vs_time, self.h_ch11_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes = self.fillChannelPlots(event, 4, self.h_b4, self.h_mcp1_ch10, self.h_ch10_vs_ch11, self.h_ch10_ch11_x_vs_ratio, self.h_b4_t, self.h_ch10_x_vs_amp, self.h_ch11_x_vs_amp, self.h_ch10_x_vs_time, self.h_ch11_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes)
            # bar 5
            self.h_b5, self.h_mcp1_ch12, self.h_ch12_vs_ch13, self.h_ch12_ch13_x_vs_ratio, self.h_b5_t, self.h_ch12_x_vs_amp, self.h_ch13_x_vs_amp, self.h_ch12_x_vs_time, self.h_ch13_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes = self.fillChannelPlots(event, 5, self.h_b5, self.h_mcp1_ch12, self.h_ch12_vs_ch13, self.h_ch12_ch13_x_vs_ratio, self.h_b5_t, self.h_ch12_x_vs_amp, self.h_ch13_x_vs_amp, self.h_ch12_x_vs_time, self.h_ch13_x_vs_time, self.h_allChannel_timing, self.h_allChannel_timingRes)

       # end filling loop     

        self.draw2Dbar(self.c1, self.h_b1, 1)
        self.draw2Dbar(self.c1, self.h_b2, 2)
        self.draw2Dbar(self.c1, self.h_b3, 3)
        self.draw2Dbar(self.c1, self.h_b4, 4)
        self.draw2Dbar(self.c1, self.h_b5, 5)
        self.draw2Dbar(self.c1, self.h_b1_t, 1, "test")
        self.draw2Dbar(self.c1, self.h_b2_t, 2, "test")
        self.draw2Dbar(self.c1, self.h_b3_t, 3, "test")
        self.draw2Dbar(self.c1, self.h_b4_t, 4, "test")
        self.draw2Dbar(self.c1, self.h_b5_t, 5, "test")
        
        self.c2.cd()
        self.c2.SetLeftMargin(0.15);
        self.c2.SetRightMargin(0.05);
        self.c2.SetBottomMargin(0.10);
        self.c2.SetTopMargin(0.05);
        # kBlack == 1, kRed == 632, kBlue == 600, kGreen == 416, kMagenta == 616
        self.h_mcp0_ch1.SetLineColor(1) # kBlack
        self.h_mcp0_ch3.SetLineColor(600) # kBlue
        self.h_mcp0_ch5.SetLineColor(632) # kRed
        self.h_mcp1_ch10.SetLineColor(416+2) # kGreen+2
        self.h_mcp1_ch12.SetLineColor(616-3) # kMagenta-3
        
        self.h_mcp0_ch1.SetLineWidth(3) 
        self.h_mcp0_ch3.SetLineWidth(3) 
        self.h_mcp0_ch5.SetLineWidth(3) 
        self.h_mcp1_ch10.SetLineWidth(3) 
        self.h_mcp1_ch12.SetLineWidth(3) 
        
        self.h_mcp1_ch12.SetYTitle("Noramlized Entries / 20 mV")
        self.h_mcp1_ch12.SetXTitle("MCP Amplitude [mV]")
        self.h_mcp1_ch12.SetTitle("")
        self.h_mcp1_ch12.DrawNormalized()
        self.h_mcp1_ch12.GetYaxis().SetRangeUser(0,1.6)
        
        self.h_mcp0_ch3.DrawNormalized("same")
        self.h_mcp0_ch5.DrawNormalized("same")
        self.h_mcp1_ch10.DrawNormalized("same")
        self.h_mcp0_ch1.DrawNormalized("same")
        
        leg = TLegend(0.5, 0.4, .85, .7);
        leg.AddEntry(self.h_mcp0_ch1, "MCP Amplitude: Bar 1 Signal", "l");
        leg.AddEntry(self.h_mcp0_ch3, "MCP Amplitude: Bar 2 Signal", "l");
        leg.AddEntry(self.h_mcp0_ch5, "MCP Amplitude: Bar 3 Signal", "l");
        leg.AddEntry(self.h_mcp1_ch10, "MCP Amplitude: Bar 4 Signal", "l");
        leg.AddEntry(self.h_mcp1_ch12, "MCP Amplitude: Bar 5 Signal", "l");
        leg.Draw("same");
        
        self.c2.Print("{0}/mcp_amplitudes.png".format(self.topDir) )

        
        self.drawLvsRinBar(self.c3, self.h_ch1_vs_ch2, 1)
        self.drawLvsRinBar(self.c3, self.h_ch3_vs_ch4, 2)
        self.drawLvsRinBar(self.c3, self.h_ch5_vs_ch6, 3)
        self.drawLvsRinBar(self.c3, self.h_ch10_vs_ch11, 4)
        self.drawLvsRinBar(self.c3, self.h_ch12_vs_ch13, 5)
        
        
        self.drawSingleProfile(self.c4, self.h_ch1_ch2_x_vs_ratio, 1, 'Right/Left')
        self.drawSingleProfile(self.c4, self.h_ch3_ch4_x_vs_ratio, 2, 'Right/Left')
        self.drawSingleProfile(self.c4, self.h_ch5_ch6_x_vs_ratio, 3, 'Right/Left')
        self.drawSingleProfile(self.c4, self.h_ch10_ch11_x_vs_ratio, 4, 'Right/Left')
        self.drawSingleProfile(self.c4, self.h_ch12_ch13_x_vs_ratio, 5, 'Right/Left')

        self.drawSingleProfile(self.c4, self.h_ch1_x_vs_amp, 1, 'Bar 1 [Right SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch2_x_vs_amp, 1, 'Bar 1 [Left SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch3_x_vs_amp, 2, 'Bar 2 [Right SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch4_x_vs_amp, 2, 'Bar 2 [Left SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch5_x_vs_amp, 3, 'Bar 3 [Right SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch6_x_vs_amp, 3, 'Bar 3 [Left SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch10_x_vs_amp, 4, 'Bar 4 [Right SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch11_x_vs_amp, 4, 'Bar 4 [Left SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch12_x_vs_amp, 5, 'Bar 5 [Right SiPM]')
        self.drawSingleProfile(self.c4, self.h_ch13_x_vs_amp, 5, 'Bar 5 [Left SiPM]')

        self.drawSingleProfile(self.c4, self.h_ch1_x_vs_time, 1, 'Bar 1 [Right SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch2_x_vs_time, 1, 'Bar 1 [Left SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch3_x_vs_time, 2, 'Bar 2 [Right SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch4_x_vs_time, 2, 'Bar 2 [Left SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch5_x_vs_time, 3, 'Bar 3 [Right SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch6_x_vs_time, 3, 'Bar 3 [Left SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch10_x_vs_time, 4, 'Bar 4 [Right SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch11_x_vs_time, 4, 'Bar 4 [Left SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch12_x_vs_time, 5, 'Bar 5 [Right SiPM] Time')
        self.drawSingleProfile(self.c4, self.h_ch13_x_vs_time, 5, 'Bar 5 [Left SiPM] Time')

        self.c4.cd()
        self.h_allChannel_timing.Draw()
        self.c4.Print( "{0}/h_allChannel_timing.png".format(self.topDir) )
        self.c4.cd()
        self.h_allChannel_timingRes.Draw()
        f_res = TF1("f_res", "gaus") # first degree polynomial --> this is just a choice atm
        f_res.SetRange(-500, 500)
        self.h_allChannel_timingRes.Fit("f_res", "")
        self.c4.Print( "{0}/h_allChannel_timingRes.png".format(self.topDir) )

        """self.drawXquadrants(self.c4, self.h_ch1_ch2_ratio_x1, self.h_ch1_ch2_ratio_x2, self.h_ch1_ch2_ratio_x3, self.h_ch1_ch2_ratio_x4, 1, 'Right/Left')
        self.drawXquadrants(self.c4, self.h_ch3_ch4_ratio_x1, self.h_ch3_ch4_ratio_x2, self.h_ch3_ch4_ratio_x3, self.h_ch3_ch4_ratio_x4, 2, 'Right/Left')
        self.drawXquadrants(self.c4, self.h_ch5_ch6_ratio_x1, self.h_ch5_ch6_ratio_x2, self.h_ch5_ch6_ratio_x3, self.h_ch5_ch6_ratio_x4, 3, 'Right/Left')
        self.drawXquadrants(self.c4, self.h_ch10_ch11_ratio_x1, self.h_ch10_ch11_ratio_x2, self.h_ch10_ch11_ratio_x3, self.h_ch10_ch11_ratio_x4, 4, 'Right/Left')
        self.drawXquadrants(self.c4, self.h_ch12_ch13_ratio_x1, self.h_ch12_ch13_ratio_x2, self.h_ch12_ch13_ratio_x3, self.h_ch12_ch13_ratio_x4, 5, 'Right/Left')"""

        """self.drawXquadrants(self.c4, self.h_ch1_x1, self.h_ch1_x2, self.h_ch1_x3, self.h_ch1_x4, 1, 'Bar 1 [Right SiPM]')
        self.drawXquadrants(self.c4, self.h_ch2_x1, self.h_ch2_x2, self.h_ch2_x3, self.h_ch2_x4, 1, 'Bar 1 [Left SiPM]')
        self.drawXquadrants(self.c4, self.h_ch3_x1, self.h_ch3_x2, self.h_ch3_x3, self.h_ch3_x4, 2, 'Bar 2 [Right SiPM]')
        self.drawXquadrants(self.c4, self.h_ch4_x1, self.h_ch4_x2, self.h_ch4_x3, self.h_ch4_x4, 2, 'Bar 2 [Left SiPM]')
        self.drawXquadrants(self.c4, self.h_ch5_x1, self.h_ch5_x2, self.h_ch5_x3, self.h_ch5_x4, 3, 'Bar 3 [Right SiPM]')
        self.drawXquadrants(self.c4, self.h_ch6_x1, self.h_ch6_x2, self.h_ch6_x3, self.h_ch6_x4, 3, 'Bar 3 [Left SiPM]')
        self.drawXquadrants(self.c4, self.h_ch10_x1, self.h_ch10_x2, self.h_ch10_x3, self.h_ch10_x4, 4, 'Bar 4 [Right SiPM]')
        self.drawXquadrants(self.c4, self.h_ch11_x1, self.h_ch11_x2, self.h_ch11_x3, self.h_ch11_x4, 4, 'Bar 4 [Left SiPM]')
        self.drawXquadrants(self.c4, self.h_ch12_x1, self.h_ch12_x2, self.h_ch12_x3, self.h_ch12_x4, 5, 'Bar 5 [Right SiPM]')
        self.drawXquadrants(self.c4, self.h_ch13_x1, self.h_ch13_x2, self.h_ch13_x3, self.h_ch13_x4, 5, 'Bar 5 [Left SiPM]')"""
        
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

    def drawSingleProfile(self, c0, h_p, barNum, varName):
        """ function to recieve canvas (c0), TProfile (h_p), and bar number"""
        c0.cd()
        c0.SetLeftMargin(0.15);
        c0.SetRightMargin(0.05);
        c0.SetBottomMargin(0.10);
        c0.SetTopMargin(0.10);
        
        if varName == 'Right/Left':
            xtitle = "Right SiPM Amplitude / Left SiPM Amplitude"
        elif 'Time' in varName:
            xtitle = '{0} Time'.format(varName)
        else:
            xtitle = '{0} Amplitude'.format(varName)
        
        h_p.SetYTitle("Noramlized Entries / Bin")
        h_p.SetXTitle(xtitle)
        h_p.SetTitle("Bar {0}".format(barNum))

        h_p.Draw()

        if varName == 'Right/Left':
            filename = '{0}/bar{1}_RL_ratio_x.png'.format(self.topDir, barNum)
        else:
            filename = '{0}/{1}_x.png'.format(self.topDir, varName.replace(' ', '_').replace('[','').replace(']',''))

        c0.Print(filename)

    # =============================
