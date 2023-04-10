#!/cvmfs/cms.cern.ch/slc7_amd64_gcc10/cms/cmssw/CMSSW_12_6_3/external/slc7_amd64_gcc10/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "histos.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();
'''histo2d = TH2F("h_2dquss","y,x,#entries", 100, -4., 4, 100, -4., 4.)

c1 = TCanvas('cHisto','cHisto',600,600)
histo = gROOT.FindObject('histo')
histo1 = gROOT.FindObject('histo1')
histo2d.Fill(histo.GetBinContent,histo1.GetBinContent)
histo2d.Draw()'''
c1 = TCanvas('cHisto','cHisto',600,600,700,500)

histo = gROOT.FindObject('histo')
histo1 = gROOT.FindObject('histo1')
histo2 = gROOT.FindObject('histo2')
histo3 = gROOT.FindObject('histo3')
#histo = gROOT.FindObject('histo2D')
#histo2 = gROOT.FindObject('histo2D1')
histo.SetLineColor(2)
#histo.SetFillColor(46)
histo1.SetLineColor(3)
#histo1.SetFillColor(30)
histo2.SetLineColor(4)
histo3.SetLineColor(5)
histo.GetYaxis().SetTitleOffset(1.5)
histo.GetYaxis().SetNdivisions(503)
histo.GetYaxis().SetTitle('Entries')
histo.GetXaxis().SetTitle('\eta')
histo.SetTitle('jet \eta')
histo.SetStats(0)
histo.Draw()
histo1.Draw("SAME")
histo2.Draw("SAME")
histo3.Draw("SAME")
c1.SetLeftMargin(0.2)
#c1.SetLogy()
#c1.SetLogx()
#c1.Update()

legend = TLegend(0.2,0.8,0.4,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("histo","\eta - pt>10GeV","f")
legend.AddEntry("histo1","\eta - pt>50GeV","f")
legend.AddEntry("histo2","\eta - pt>100GeV","f")
legend.AddEntry("histo3","\eta - pt>200GeV","f")
#legend.AddEntry("histo1","l1phi","f")
#legend.AddEntry("histo2","p","f")
legend.Draw()
c1.Print("histo_prez_jet_obceta.png")
input('press enter to exit')
