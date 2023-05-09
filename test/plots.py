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

peta = gROOT.FindObject('peta')
peta2 = gROOT.FindObject('peta2')
peta3 = gROOT.FindObject('peta3')
peta4 = gROOT.FindObject('peta4')
peta5 = gROOT.FindObject('peta5')
peta6 = gROOT.FindObject('peta6')
'''
histo = gROOT.FindObject('histo3')
histo1 = gROOT.FindObject('histo4')
histo2 = gROOT.FindObject('histo5')
#histo3 = gROOT.FindObject('histo3')
#histo = gROOT.FindObject('histo2D')
#histo2 = gROOT.FindObject('histo2D1')
#histo.SetLineColor(2)
#histo.SetFillColor(46)
#histo1.SetLineColor(3)
#histo1.SetFillColor(30)
#histo2.SetLineColor(4)
#histo3.SetLineColor(5)
histo.GetYaxis().SetTitleOffset(1.5)
histo.GetYaxis().SetNdivisions(503)
histo.GetYaxis().SetTitle('Entries')
histo.GetXaxis().SetTitle('number of muons')
#histo.GetXaxis().SetTitle('vertexz')
histo.SetTitle('number of muons')
histo.SetStats(0)
'''
peta.Draw("AP")
peta.SetTitle('Final countdown Jet; jet \eta; Probability')
peta.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
peta.SetLineColor(3)
peta3.Draw("SAME")
peta3.SetMarkerColor(2)
peta3.SetLineColor(2)
peta4.Draw("SAME")
'''
histo.Draw()
histo.SetLineColor(2)
histo1.Draw("SAME")
histo1.SetLineColor(3)
histo2.Draw("SAME")
histo2.SetLineColor(4)
#histo3.Draw("SAME")
'''
c1.SetLeftMargin(0.2)
c1.SetGridy()
c1.SetGridx()
'''
#c1.SetLogy()
#c1.SetLogx()
#c1.Update()
'''
legend = TLegend(0.7,0.8,0.9,0.9)
legend.SetTextSize(0.03)
#legend.AddEntry("histo","\sqrt{s}","f")
#legend.AddEntry("histo","vertex z diffrence","f")
#legend.AddEntry("peta5","L1T to PAT muon pt>5","f")
#legend.AddEntry("peta6","PAT to L1T muon pt>5","f")
legend.AddEntry("peta","PAT muon pt>5","f")
legend.AddEntry("peta3","PAT muon pt>10","f")
legend.AddEntry("peta4","PAT muon pt>20","f")
#legend.AddEntry("histo","P_{t}>5GeV","f")
#legend.AddEntry("histo1","P_{t}>10GeV","f")
#legend.AddEntry("histo2","P_{t}>20GeV","f")
#legend.AddEntry("histo3","\eta - pt>200GeV","f")
#legend.AddEntry("histo1","l1phi","f")
#legend.AddEntry("histo2","p","f")
legend.Draw()
c1.Print("histo_glownyrysunek8test.png")
input('press enter to exit')
