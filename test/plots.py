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

pjeta5 = gROOT.FindObject('pjeta5')
pjeta10 = gROOT.FindObject('pjeta10')
pjeta20 = gROOT.FindObject('pjeta20')
pjeta1_5 = gROOT.FindObject('pjeta1_5')
pjeta1_10 = gROOT.FindObject('pjeta1_10')
pjeta1_20 = gROOT.FindObject('pjeta1_20')
pmeta50 = gROOT.FindObject('pmeta50')
pmeta100 = gROOT.FindObject('pmeta100')
pmeta200 = gROOT.FindObject('pmeta200')
pmeta1_50 = gROOT.FindObject('pmeta1_50')
pmeta1_100 = gROOT.FindObject('pmeta1_100')
pmeta1_200 = gROOT.FindObject('pmeta1_200')
pml1t5 = gROOT.FindObject('pml1t5')
pl1tm5 = gROOT.FindObject('pl1tm5')

'''
histo = gROOT.FindObject('vertex')
#histo1 = gROOT.FindObject('jeteta50')
#histo2 = gROOT.FindObject('jeteta100')
#histo3 = gROOT.FindObject('jeteta200')
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
histo.GetXaxis().SetTitle('z [cm]')
#histo.GetXaxis().SetTitle('vertexz')
histo.SetTitle('\Delta~vertex~of~two~muons')
histo.SetStats(0)
'''
pmeta50.Draw("AP")
pmeta50.SetTitle('Probabilty; muon \eta; Probability')
pmeta50.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pmeta50.SetLineColor(3)
pmeta100.Draw("SAME")
pmeta100.SetMarkerColor(2)
pmeta100.SetLineColor(2)
pmeta200.Draw("SAME")
'''
histo.Draw()
#histo.SetLineColor(2)
#histo1.Draw("SAME")
#histo1.SetLineColor(3)
#histo2.Draw("SAME")
#histo2.SetLineColor(4)
#histo3.Draw("SAME")
#histo3.SetLineColor(5)
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
#legend.AddEntry("vertex","\Delta vertex","f")
#legend.AddEntry("jeteta50","\eta - p_{t}>50GeV","f")
#legend.AddEntry("jeteta100","\eta - p_{t}>100GeV","f")
#legend.AddEntry("jeteta200","\eta - p_{t}>200GeV","f")
legend.AddEntry("pmeta50","jet pt>50","f")
legend.AddEntry("pmeta100","jet pt>100","f")
legend.AddEntry("pmeta200","jet pt>200","f")
#legend.AddEntry("peta3","PAT muon pt>10","f")
#legend.AddEntry("peta4","PAT muon pt>20","f")
#legend.AddEntry("histo","P_{t}>5GeV","f")
#legend.AddEntry("histo1","P_{t}>10GeV","f")
#legend.AddEntry("histo2","P_{t}>20GeV","f")
#legend.AddEntry("histo3","\eta - pt>200GeV","f")
#legend.AddEntry("histo1","l1phi","f")
#legend.AddEntry("histo2","p","f")
legend.Draw()
c1.Print("probability_muon_jet_cats.png")
input('press enter to exit')
