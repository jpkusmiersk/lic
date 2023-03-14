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

c1 = TCanvas('cHisto','cHisto',600,600)
histo = gROOT.FindObject('histo')
histo1 = gROOT.FindObject('histo1')
histo2 = gROOT.FindObject('histo2')
histo.SetLineColor(2)
#histo.SetFillColor(46)
histo1.SetLineColor(3)
#histo1.SetFillColor(30)
histo2.SetLineColor(4)
histo.Draw()
histo1.Draw("SAME")
histo2.Draw("SAME")
legend = TLegend(0.1,0.7,0.48,0.9)
legend.AddEntry("histo","vx","f")
legend.AddEntry("histo1","vy","f")
legend.AddEntry("histo2","vz","f")
legend.Draw()
c1.Print("histo_v.png")
input('press enter to exit')
