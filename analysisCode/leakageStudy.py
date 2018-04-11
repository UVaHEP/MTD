# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 5, 2018
#Purpose: Script to combine executables for 1) converting testbeam .dat files to .root, and 2) analyzing root file and producing timing-relevant information by channel

import os,sys, argparse
from ROOT import gROOT, TH1D, TFile, TTree, TChain, TCanvas, TH2D, gStyle, TLegend

def draw2Dbar(c0, h0, barNum, test=""):
    """ function to recieve canvas (c0), histogram (h0), and name for 2D bar plot"""
    c0.cd()
    h0.SetTitle( "Bar {0}".format(barNum) )
    h0.SetXTitle("X [mm]")
    h0.SetYTitle("Y [mm]")
    h0.Draw("colz")

    #h0.Integral(x1, x2, y1, y2)
    #7-10, 12-14, 16
    ymin = 4*barNum + 3
    ymax = 4*barNum + 7
    window = h0.Integral(4, 38, ymin, ymax)

    if h0.GetEntries() > 0:
        print "In bar {0}: {1}\t Total: {2}\t % in Bar: {3}".format(barNum, window, h0.GetEntries(), window/h0.GetEntries())
    else:
        print "no entries"

    c0.Print( "bar{0}_python{1}.png".format(barNum, test) )


def drawLvsRinBar(c0, h0, barNum):
    """ function to recieve canvas (c0), histogram (h0), and name for 2D bar plot"""
    c0.cd()
    c0.SetLeftMargin(0.15);
    h0.SetTitle( "Bar {0}: R vs L Amplitudes".format(barNum) )
    h0.SetXTitle("Right SiPM [mV]")
    h0.SetYTitle("Left SiPM [mV]")
    h0.Draw("colz")

    c0.Print( "bar{0}_rightVleft.png".format(barNum) )


def drawXquadrants(c0, h_x1, h_x2, h_x3, h_x4, barNum):
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
    
    h_x4.SetYTitle("Noramlized Entries / Bin")
    h_x4.SetXTitle("Right SiPM Amplitude / Left SiPM Amplitude")
    h_x4.SetTitle("Bar {0}".format(barNum))
    h_x4.DrawNormalized()
    h_x4.GetYaxis().SetRangeUser(0,1.6)
    
    h_x3.DrawNormalized("same")
    h_x2.DrawNormalized("same")
    h_x1.DrawNormalized("same")
    
    leg = TLegend(0.6, 0.5, .85, .7);
    leg.AddEntry(h_x1, "Right/Left Amplitude: Q1", "l");
    leg.AddEntry(h_x2, "Right/Left Amplitude: Q2", "l");
    leg.AddEntry(h_x3, "Right/Left Amplitude: Q3", "l");
    leg.AddEntry(h_x4, "Right/Left Amplitude: Q4", "l");
    leg.Draw("same");
    
    c0.Print("bar{0}_RL_ratio_xq.png".format(barNum))



###########################################################


f0 = TFile('/home/btannenw/Desktop/MTD/testbeam_FNAL_03-2018/all5exposure/DataCMSVMETiming_5barExposure.root', 'READ')
t0 = f0.pulse


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

h_mcp0_ch1 = TH1D("h_mcp0_bar1", "h_mcp0_ch1", 25, 0, 500);
h_mcp0_ch3 = TH1D("h_mcp0_bar3", "h_mcp0_ch3", 25, 0, 500);
h_mcp0_ch5 = TH1D("h_mcp0_bar5", "h_mcp0_ch5", 25, 0, 500);
h_mcp1_ch10 = TH1D("h_mcp1_ch10", "h_mcp1_ch10", 25, 0, 500);
h_mcp1_ch12 = TH1D("h_mcp1_ch12", "h_mcp1_ch12", 25, 0, 500);

h_ch1_vs_ch2 = TH2D("h_ch1_vs_ch2", "h_ch1_vs_ch2", 22, 0, 1100, 22, 0, 1100)
h_ch3_vs_ch4 = TH2D("h_ch3_vs_ch4", "h_ch3_vs_ch4", 22, 0, 1100, 22, 0, 1100)
h_ch5_vs_ch6 = TH2D("h_ch5_vs_ch6", "h_ch5_vs_ch6", 22, 0, 1100, 22, 0, 1100)
h_ch10_vs_ch11 = TH2D("h_ch10_vs_ch11", "h_ch10_vs_ch11", 22, 0, 1100, 22, 0, 1100)
h_ch12_vs_ch13 = TH2D("h_ch12_vs_ch13", "h_ch12_vs_ch13", 22, 0, 1100, 22, 0, 1100)

h_ch1_ch2_ratio_xq1 = TH1D("h_ch1_ch2_ratio_xq1", "h_ch1_ch2_ratio_xq1", 100, 0.0, 2.0)
h_ch1_ch2_ratio_xq2 = TH1D("h_ch1_ch2_ratio_xq2", "h_ch1_ch2_ratio_xq2", 100, 0.0, 2.0)
h_ch1_ch2_ratio_xq3 = TH1D("h_ch1_ch2_ratio_xq3", "h_ch1_ch2_ratio_xq3", 100, 0.0, 2.0)
h_ch1_ch2_ratio_xq4 = TH1D("h_ch1_ch2_ratio_xq4", "h_ch1_ch2_ratio_xq4", 100, 0.0, 2.0)
h_ch3_ch4_ratio_xq1 = TH1D("h_ch3_ch4_ratio_xq1", "h_ch3_ch4_ratio_xq1", 100, 0.0, 2.0)
h_ch3_ch4_ratio_xq2 = TH1D("h_ch3_ch4_ratio_xq2", "h_ch3_ch4_ratio_xq2", 100, 0.0, 2.0)
h_ch3_ch4_ratio_xq3 = TH1D("h_ch3_ch4_ratio_xq3", "h_ch3_ch4_ratio_xq3", 100, 0.0, 2.0)
h_ch3_ch4_ratio_xq4 = TH1D("h_ch3_ch4_ratio_xq4", "h_ch3_ch4_ratio_xq4", 100, 0.0, 2.0)
h_ch5_ch6_ratio_xq1 = TH1D("h_ch5_ch6_ratio_xq1", "h_ch5_ch6_ratio_xq1", 100, 0.0, 2.0)
h_ch5_ch6_ratio_xq2 = TH1D("h_ch5_ch6_ratio_xq2", "h_ch5_ch6_ratio_xq2", 100, 0.0, 2.0)
h_ch5_ch6_ratio_xq3 = TH1D("h_ch5_ch6_ratio_xq3", "h_ch5_ch6_ratio_xq3", 100, 0.0, 2.0)
h_ch5_ch6_ratio_xq4 = TH1D("h_ch5_ch6_ratio_xq4", "h_ch5_ch6_ratio_xq4", 100, 0.0, 2.0)
h_ch10_ch11_ratio_xq1 = TH1D("h_ch10_ch11_ratio_xq1", "h_ch10_ch11_ratio_xq1", 100, 0.0, 2.0)
h_ch10_ch11_ratio_xq2 = TH1D("h_ch10_ch11_ratio_xq2", "h_ch10_ch11_ratio_xq2", 100, 0.0, 2.0)
h_ch10_ch11_ratio_xq3 = TH1D("h_ch10_ch11_ratio_xq3", "h_ch10_ch11_ratio_xq3", 100, 0.0, 2.0)
h_ch10_ch11_ratio_xq4 = TH1D("h_ch10_ch11_ratio_xq4", "h_ch10_ch11_ratio_xq4", 100, 0.0, 2.0)
h_ch12_ch13_ratio_xq1 = TH1D("h_ch12_ch13_ratio_xq1", "h_ch12_ch13_ratio_xq1", 100, 0.0, 2.0)
h_ch12_ch13_ratio_xq2 = TH1D("h_ch12_ch13_ratio_xq2", "h_ch12_ch13_ratio_xq2", 100, 0.0, 2.0)
h_ch12_ch13_ratio_xq3 = TH1D("h_ch12_ch13_ratio_xq3", "h_ch12_ch13_ratio_xq3", 100, 0.0, 2.0)
h_ch12_ch13_ratio_xq4 = TH1D("h_ch12_ch13_ratio_xq4", "h_ch12_ch13_ratio_xq4", 100, 0.0, 2.0)


c1 = TCanvas("c1", "c1", 800, 800)
c2 = TCanvas("c2", "c2", 800, 800)
c3 = TCanvas("c3", "c3", 800, 800)
c4 = TCanvas("c4", "c4", 800, 800)

gStyle.SetOptStat(0000)

print t0.GetEntries(), f0.pulse.GetEntries()
nTotal=0
for event in t0:

    #if (nTotal > 10000):
    #    break

    nTotal += 1
    if (nTotal % 10000 == 0):
        print nTotal, "processed"

    if (event.amp[1] > 100 and event.amp[3] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) ):
        h_b1.Fill(event.x_dut[2], event.y_dut[2])
        h_mcp0_ch1.Fill( event.amp[0] )
        h_ch1_vs_ch2.Fill( event.amp[1], event.amp[2] )

        if( event.x_dut[2]>=-2 and event.x_dut[2] < 6):
            h_ch1_ch2_ratio_xq1.Fill( event.amp[1] / event.amp[2] )
        if( event.x_dut[2]>= 6 and event.x_dut[2] < 15):
            h_ch1_ch2_ratio_xq2.Fill( event.amp[1] / event.amp[2] )
        if( event.x_dut[2]>=14 and event.x_dut[2] < 24):
            h_ch1_ch2_ratio_xq3.Fill( event.amp[1] / event.amp[2] )
        if( event.x_dut[2]>=22 and event.x_dut[2] < 33):
            h_ch1_ch2_ratio_xq4.Fill( event.amp[1] / event.amp[2] )

    if (event.amp[3] > 100 and event.amp[5] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) ):
        h_b2.Fill(event.x_dut[2], event.y_dut[2])
        h_mcp0_ch3.Fill( event.amp[0] )
        h_ch3_vs_ch4.Fill( event.amp[3], event.amp[4] )
        
        if( event.x_dut[2]>=-2 and event.x_dut[2] < 6):
            h_ch3_ch4_ratio_xq1.Fill( event.amp[3] / event.amp[4] )
        if( event.x_dut[2]>= 6 and event.x_dut[2] < 15):
            h_ch3_ch4_ratio_xq2.Fill( event.amp[3] / event.amp[4] )
        if( event.x_dut[2]>=14 and event.x_dut[2] < 24):
            h_ch3_ch4_ratio_xq3.Fill( event.amp[3] / event.amp[4] )
        if( event.x_dut[2]>=22 and event.x_dut[2] < 33):
            h_ch3_ch4_ratio_xq4.Fill( event.amp[3] / event.amp[4] )

    if (event.amp[5] > 100 and event.amp[10] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) ):
        h_b3.Fill(event.x_dut[2], event.y_dut[2])
        h_mcp0_ch5.Fill( event.amp[0] )
        h_ch5_vs_ch6.Fill( event.amp[5], event.amp[6] )
        
        if( event.x_dut[2]>=-2 and event.x_dut[2] < 6):
            h_ch5_ch6_ratio_xq1.Fill( event.amp[5] / event.amp[6] )
        if( event.x_dut[2]>= 6 and event.x_dut[2] < 15):
            h_ch5_ch6_ratio_xq2.Fill( event.amp[5] / event.amp[6] )
        if( event.x_dut[2]>=14 and event.x_dut[2] < 24):
            h_ch5_ch6_ratio_xq3.Fill( event.amp[5] / event.amp[6] )
        if( event.x_dut[2]>=22 and event.x_dut[2] < 33):
            h_ch5_ch6_ratio_xq4.Fill( event.amp[5] / event.amp[6] )
        
    if (event.amp[10] > 100 and event.amp[5] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) ):
        h_b4.Fill(event.x_dut[2], event.y_dut[2])
        h_mcp1_ch10.Fill( event.amp[9] )
        h_ch10_vs_ch11.Fill( event.amp[10], event.amp[11] )

        if( event.x_dut[2]>=-2 and event.x_dut[2] < 6):
            h_ch10_ch11_ratio_xq1.Fill( event.amp[10] / event.amp[11] )
        if( event.x_dut[2]>= 6 and event.x_dut[2] < 15):
            h_ch10_ch11_ratio_xq2.Fill( event.amp[10] / event.amp[11] )
        if( event.x_dut[2]>=14 and event.x_dut[2] < 24):
            h_ch10_ch11_ratio_xq3.Fill( event.amp[10] / event.amp[11] )
        if( event.x_dut[2]>=22 and event.x_dut[2] < 33):
            h_ch10_ch11_ratio_xq4.Fill( event.amp[10] / event.amp[11] )
        
    if (event.amp[12] > 100 and event.amp[10] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) ):
        h_b5.Fill(event.x_dut[2], event.y_dut[2])
        h_mcp1_ch12.Fill( event.amp[9] )
        h_ch12_vs_ch13.Fill( event.amp[12], event.amp[13] )

        if( event.x_dut[2]>=-2 and event.x_dut[2] < 6):
            h_ch12_ch13_ratio_xq1.Fill( event.amp[12] / event.amp[13] )
        if( event.x_dut[2]>= 6 and event.x_dut[2] < 15):
            h_ch12_ch13_ratio_xq2.Fill( event.amp[12] / event.amp[13] )
        if( event.x_dut[2]>=14 and event.x_dut[2] < 24):
            h_ch12_ch13_ratio_xq3.Fill( event.amp[12] / event.amp[13] )
        if( event.x_dut[2]>=22 and event.x_dut[2] < 33):
            h_ch12_ch13_ratio_xq4.Fill( event.amp[12] / event.amp[13] )
        
    if (event.amp[1] > 100 and event.amp[3] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) and abs(event.y_dut[2] -8.5) < 2 and event.x_dut[2]>=-2 and event.x_dut[2]<=33):
        h_b1_t.Fill(event.x_dut[2], event.y_dut[2])
    if (event.amp[3] > 100 and event.amp[5] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) and abs(event.y_dut[2] -12.5) < 2 and event.x_dut[2]>=-2 and event.x_dut[2]<=33):
        h_b2_t.Fill(event.x_dut[2], event.y_dut[2])
    if (event.amp[5] > 100 and event.amp[10] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) and abs(event.y_dut[2] -16.5) < 2 and event.x_dut[2]>=-2 and event.x_dut[2]<=33):
        h_b3_t.Fill(event.x_dut[2], event.y_dut[2])
    if (event.amp[10] > 100 and event.amp[5] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) and abs(event.y_dut[2] -20.5) < 2 and event.x_dut[2]>=-2 and event.x_dut[2]<=33):
        h_b4_t.Fill(event.x_dut[2], event.y_dut[2])
    if (event.amp[12] > 100 and event.amp[10] < 100 and abs(event.xSlope) < 0.0004 and abs(event.ySlope) and abs(event.y_dut[2] -24.5) < 2 and event.x_dut[2]>=-2 and event.x_dut[2]<=33):
        h_b5_t.Fill(event.x_dut[2], event.y_dut[2])

draw2Dbar(c1, h_b1, 1)
draw2Dbar(c1, h_b2, 2)
draw2Dbar(c1, h_b3, 3)
draw2Dbar(c1, h_b4, 4)
draw2Dbar(c1, h_b5, 5)
draw2Dbar(c1, h_b1_t, 1, "test")
draw2Dbar(c1, h_b2_t, 2, "test")
draw2Dbar(c1, h_b3_t, 3, "test")
draw2Dbar(c1, h_b4_t, 4, "test")
draw2Dbar(c1, h_b5_t, 5, "test")

c2.cd()
c2.SetLeftMargin(0.15);
c2.SetRightMargin(0.05);
c2.SetBottomMargin(0.10);
c2.SetTopMargin(0.05);
# kBlack == 1, kRed == 632, kBlue == 600, kGreen == 416, kMagenta == 616
h_mcp0_ch1.SetLineColor(1) # kBlack
h_mcp0_ch3.SetLineColor(600) # kBlue
h_mcp0_ch5.SetLineColor(632) # kRed
h_mcp1_ch10.SetLineColor(416+2) # kGreen+2
h_mcp1_ch12.SetLineColor(616-3) # kMagenta-3

h_mcp0_ch1.SetLineWidth(3) 
h_mcp0_ch3.SetLineWidth(3) 
h_mcp0_ch5.SetLineWidth(3) 
h_mcp1_ch10.SetLineWidth(3) 
h_mcp1_ch12.SetLineWidth(3) 

h_mcp1_ch12.SetYTitle("Noramlized Entries / 20 mV")
h_mcp1_ch12.SetXTitle("MCP Amplitude [mV]")
h_mcp1_ch12.SetTitle("")
h_mcp1_ch12.DrawNormalized()
h_mcp1_ch12.GetYaxis().SetRangeUser(0,1.6)

h_mcp0_ch3.DrawNormalized("same")
h_mcp0_ch5.DrawNormalized("same")
h_mcp1_ch10.DrawNormalized("same")
h_mcp0_ch1.DrawNormalized("same")

leg = TLegend(0.5, 0.4, .85, .7);
leg.AddEntry(h_mcp0_ch1, "MCP Amplitude: Bar 1 Signal", "l");
leg.AddEntry(h_mcp0_ch3, "MCP Amplitude: Bar 2 Signal", "l");
leg.AddEntry(h_mcp0_ch5, "MCP Amplitude: Bar 3 Signal", "l");
leg.AddEntry(h_mcp1_ch10, "MCP Amplitude: Bar 4 Signal", "l");
leg.AddEntry(h_mcp1_ch12, "MCP Amplitude: Bar 5 Signal", "l");
leg.Draw("same");

c2.Print("mcp_amplitudes.png")


drawLvsRinBar(c3, h_ch1_vs_ch2, 1)
drawLvsRinBar(c3, h_ch3_vs_ch4, 2)
drawLvsRinBar(c3, h_ch5_vs_ch6, 3)
drawLvsRinBar(c3, h_ch10_vs_ch11, 4)
drawLvsRinBar(c3, h_ch12_vs_ch13, 5)


drawXquadrants(c4, h_ch1_ch2_ratio_xq1, h_ch1_ch2_ratio_xq2, h_ch1_ch2_ratio_xq3, h_ch1_ch2_ratio_xq4, 1)
drawXquadrants(c4, h_ch3_ch4_ratio_xq1, h_ch3_ch4_ratio_xq2, h_ch3_ch4_ratio_xq3, h_ch3_ch4_ratio_xq4, 2)
drawXquadrants(c4, h_ch5_ch6_ratio_xq1, h_ch5_ch6_ratio_xq2, h_ch5_ch6_ratio_xq3, h_ch5_ch6_ratio_xq4, 3)
drawXquadrants(c4, h_ch10_ch11_ratio_xq1, h_ch10_ch11_ratio_xq2, h_ch10_ch11_ratio_xq3, h_ch10_ch11_ratio_xq4, 4)
drawXquadrants(c4, h_ch12_ch13_ratio_xq1, h_ch12_ch13_ratio_xq2, h_ch12_ch13_ratio_xq3, h_ch12_ch13_ratio_xq4, 5)


