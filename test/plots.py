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
jetpt = gROOT.FindObject('jetpt')
jetpt.GetYaxis().SetTitleOffset(1.5)
jetpt.GetYaxis().SetNdivisions(503)
jetpt.GetYaxis().SetTitle('Entries')
jetpt.GetXaxis().SetTitle('p_{T}[GeV]')
#histo.GetXaxis().SetTitle('vertexz')
jetpt.SetTitle('Transverse momentum of jets')
jetpt.SetStats(0)
jetpt.SetLineColor(2)
jetpt.Draw()
c1.SetLogy()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetpt","p_{T}","f")
#legend.Draw()
c1.Print("jet_transvers_momentum.png")

#/////////////////////////////////////////////////////

c2 = TCanvas('cHisto','cHisto',600,600,700,500)
muonpt = gROOT.FindObject('muonpt')
muonpt.GetYaxis().SetTitleOffset(1.5)
muonpt.GetYaxis().SetNdivisions(503)
muonpt.GetYaxis().SetTitle('Entries')
muonpt.GetXaxis().SetTitle('p_{T}[GeV]')
#histo.GetXaxis().SetTitle('vertexz')
muonpt.SetTitle('Transverse momentum of muons')
muonpt.SetStats(0)
muonpt.Draw()
c2.SetLogy()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("muonpt","p_{T}","f")
#legend.Draw()
c2.Print("muon_transvers_momentum.png")

#///////////////////////////////////////////////////

c3 = TCanvas('cHisto','cHisto',600,600,700,500)
jetdphi = gROOT.FindObject('jetdphi')
jetdphi.GetYaxis().SetTitleOffset(1.5)
jetdphi.GetYaxis().SetNdivisions(503)
jetdphi.GetYaxis().SetTitle('Entries')
jetdphi.GetXaxis().SetTitle('\Delta\phi')
#histo.GetXaxis().SetTitle('vertexz')
jetdphi.SetTitle('\Delta\phi~of~two~jets')
jetdphi.SetStats(0)
jetdphi.SetLineColor(3)
jetdphi.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c3.Print("jets_delta_phi.png")

#///////////////////////////////////////////////////

c4 = TCanvas('cHisto','cHisto',600,600,700,500)
jeteta10 = gROOT.FindObject('jeteta10')
jeteta50 = gROOT.FindObject('jeteta50')
jeteta100 = gROOT.FindObject('jeteta100')
jeteta200 = gROOT.FindObject('jeteta200')
jeteta10.GetYaxis().SetTitleOffset(1.5)
jeteta10.GetYaxis().SetNdivisions(503)
jeteta10.GetYaxis().SetTitle('Entries')
jeteta10.GetXaxis().SetTitle('\Delta\phi')
#histo.GetXaxis().SetTitle('vertexz')
jeteta10.SetTitle('\Delta\phi~of~two~jets')
jeteta10.SetStats(0)
jeteta10.SetLineColor(2)
jeteta10.Draw()
jeteta50.SetLineColor(3)
jeteta50.Draw("SAME")
jeteta100.SetLineColor(4)
jeteta100.Draw('SAME')
jeteta200.SetLineColor(5)
jeteta200.Draw("SAME")
legend = TLegend(0.4,0.7,0.6,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jeteta10","p_{T}>10GeV","f")
legend.AddEntry("jeteta50","p_{T}>50GeV","f")
legend.AddEntry("jeteta100","p_{T}>100GeV","f")
legend.AddEntry("jeteta200","p_{T}>200GeV","f")
legend.Draw()
c4.Print("jet_eta_ceica.png")

#/////////////////////////////////////////////////

c5 = TCanvas('cHisto','cHisto',600,600,700,500)
deltar = gROOT.FindObject('deltar')
deltar.GetYaxis().SetTitleOffset(1.5)
deltar.GetYaxis().SetNdivisions(503)
deltar.GetYaxis().SetTitle('Entries')
deltar.GetXaxis().SetTitle('\Delta R')
#histo.GetXaxis().SetTitle('vertexz')
deltar.SetTitle('\Delta R~between~muons~and~jets')
deltar.SetStats(0)
deltar.SetLineColor(2)
deltar.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c5.Print("muon_jet_deltaR.png")

#///////////////////////////////////////////////

c6 = TCanvas('cHisto','cHisto',600,600,700,500)
deltarpik = gROOT.FindObject('deltarpik')
deltarpik.GetYaxis().SetTitleOffset(1.5)
deltarpik.GetYaxis().SetNdivisions(503)
deltarpik.GetYaxis().SetTitle('Entries')
deltarpik.GetXaxis().SetTitle('\Delta R')
#histo.GetXaxis().SetTitle('vertexz')
deltarpik.SetTitle('\Delta R_{min}~between~muons~and~jets')
deltarpik.SetStats(0)
deltarpik.SetLineColor(2)
deltarpik.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c6.Print("muon_jet_deltaRpik.png")

#///////////////////////////////////////////////

c7 = TCanvas('cHisto','cHisto',600,600,700,500)
masan2m5_15_50 = gROOT.FindObject('masan2m5_15_50')
masan2m5_15_50.GetYaxis().SetTitleOffset(1.5)
masan2m5_15_50.GetYaxis().SetNdivisions(503)
masan2m5_15_50.GetYaxis().SetTitle('Entries')
masan2m5_15_50.GetXaxis().SetTitle('M [GeV]')
#histo.GetXaxis().SetTitle('vertexz')
masan2m5_15_50.SetTitle('Invariant mass of two muons')
masan2m5_15_50.SetStats(0)
masan2m5_15_50.SetLineColor(4)
masan2m5_15_50.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c7.Print("muon_2muon_masa_n051550.png")

#////////////////////////////////////////////

c8 = TCanvas('cHisto','cHisto',600,600,700,500)
masan2m5_15_100 = gROOT.FindObject('masan2m5_15_100')
masan2m5_15_100.GetYaxis().SetTitleOffset(1.5)
masan2m5_15_100.GetYaxis().SetNdivisions(503)
masan2m5_15_100.GetYaxis().SetTitle('Entries')
masan2m5_15_100.GetXaxis().SetTitle('M [GeV]')
#histo.GetXaxis().SetTitle('vertexz')
masan2m5_15_100.SetTitle('Invariant mass of two muons')
masan2m5_15_100.SetStats(0)
masan2m5_15_100.SetLineColor(4)
masan2m5_15_100.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c8.Print("muon_2muon_masa_n0515100.png")

#//////////////////////////////////////////

c9 = TCanvas('cHisto','cHisto',600,600,700,500)
masan2m5_15_200 = gROOT.FindObject('masan2m5_15_200')
masan2m5_15_200.GetYaxis().SetTitleOffset(1.5)
masan2m5_15_200.GetYaxis().SetNdivisions(503)
masan2m5_15_200.GetYaxis().SetTitle('Entries')
masan2m5_15_200.GetXaxis().SetTitle('M [GeV]')
#histo.GetXaxis().SetTitle('vertexz')
masan2m5_15_200.SetTitle('Invariant mass of two muons')
masan2m5_15_200.SetStats(0)
masan2m5_15_200.SetLineColor(4)
masan2m5_15_200.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c9.Print("muon_2muon_masa_n0515200.png")

#////////////////////////////////////////

c10 = TCanvas('cHisto','cHisto',600,600,700,500)
masan2m0_12 = gROOT.FindObject('masan2m0_12')
masan2m0_12.GetYaxis().SetTitleOffset(1.5)
masan2m0_12.GetYaxis().SetNdivisions(503)
masan2m0_12.GetYaxis().SetTitle('Entries')
masan2m0_12.GetXaxis().SetTitle('M [GeV]')
#histo.GetXaxis().SetTitle('vertexz')
masan2m0_12.SetTitle('Invariant mass of two muons')
masan2m0_12.SetStats(0)
masan2m0_12.SetLineColor(4)
masan2m0_12.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c10.Print("muon_2muon_masa_012.png")

#///////////////////////////////////////////

c11 = TCanvas('cHisto','cHisto',600,600,700,500)
masan2m0_120 = gROOT.FindObject('masan2m0_120')
masan2m0_120.GetYaxis().SetTitleOffset(1.5)
masan2m0_120.GetYaxis().SetNdivisions(503)
masan2m0_120.GetYaxis().SetTitle('Entries')
masan2m0_120.GetXaxis().SetTitle('M [GeV]')
masan2m0_120.SetTitle('Invariant mass of two muons')
masan2m0_120.SetStats(0)
masan2m0_120.SetLineColor(4)
masan2m0_120.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c11.Print("muon_2muon_masa_0120.png")

#///////////////////////////////////////////

c12 = TCanvas('cHisto','cHisto',600,600,700,500)
vertex = gROOT.FindObject('vertex')
vertex.GetYaxis().SetTitleOffset(1.5)
vertex.GetYaxis().SetNdivisions(503)
vertex.GetYaxis().SetTitle('Entries')
vertex.GetXaxis().SetTitle('z [cm]')
vertex.SetTitle('Difference of vertex in Oz')
vertex.SetStats(0)
vertex.SetLineColor(3)
vertex.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c12.Print("muon2_vertex.png")

#/////////////////////////////////////////

c13 = TCanvas('cHisto','cHisto',600,600,700,500)
muonvsjeteta = gROOT.FindObject('muonvsjeteta')
muonvsjeteta.GetYaxis().SetTitleOffset(1.5)
muonvsjeteta.GetYaxis().SetNdivisions(503)
muonvsjeteta.GetYaxis().SetTitle('muons\eta')
muonvsjeteta.GetXaxis().SetTitle('jets\eta')
muonvsjeteta.SetTitle('Muons~\eta~vs~jets~\eta')
muonvsjeteta.SetStats(0)
muonvsjeteta.SetLineColor(3)
muonvsjeteta.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c13.Print("2D_muon_vs_jet_eta.png")

#/////////////////////////////////////////////

c14 = TCanvas('cHisto','cHisto',600,600,700,500)
muonvsjetetamin = gROOT.FindObject('muonvsjetetamin')
muonvsjetetamin.GetYaxis().SetTitleOffset(1.5)
muonvsjetetamin.GetYaxis().SetNdivisions(503)
muonvsjetetamin.GetYaxis().SetTitle('muons\eta')
muonvsjetetamin.GetXaxis().SetTitle('jets\eta')
muonvsjetetamin.SetTitle('Muons~\eta~vs~jets~\eta')
muonvsjetetamin.SetStats(0)
muonvsjetetamin.SetLineColor(3)
muonvsjetetamin.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c14.Print("2D_muon_vs_jet_eta_min.png")

#////////////////////////////////////////////

c15 = TCanvas('cHisto','cHisto',600,600,700,500)
deltaphivspt = gROOT.FindObject('deltaphivspt')
deltaphivspt.GetYaxis().SetTitleOffset(1.5)
deltaphivspt.GetYaxis().SetNdivisions(503)
deltaphivspt.GetYaxis().SetTitle('\Delta\phi')
deltaphivspt.GetXaxis().SetTitle('jets p_{T}[GeV]')
deltaphivspt.SetTitle('Two~jets:~\Delta\phi~vs~p_{T} ')
deltaphivspt.SetStats(0)
deltaphivspt.SetLineColor(3)
deltaphivspt.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c15.Print("2D_jet_deltaphi_vs_pt.png")

#///////////////////////////////////////////////

c16 = TCanvas('cHisto','cHisto',600,600,700,500)
deltarvspt = gROOT.FindObject('deltarvspt')
deltarvspt.GetYaxis().SetTitleOffset(1.5)
deltarvspt.GetYaxis().SetNdivisions(503)
deltarvspt.GetYaxis().SetTitle('\Delta R')
deltarvspt.GetXaxis().SetTitle('jets p_{T}[GeV]')
deltarvspt.SetTitle('Jets:~\Delta R~vs~p_{T} ')
deltarvspt.SetStats(0)
deltarvspt.SetLineColor(3)
deltarvspt.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c16.Print("2D_jet_deltar_vs_pt.png")

#////////////////////////////////////////////

c17 = TCanvas('cHisto','cHisto',600,600,700,500)
deltarvsptmin = gROOT.FindObject('deltarvsptmin')
deltarvsptmin.GetYaxis().SetTitleOffset(1.5)
deltarvsptmin.GetYaxis().SetNdivisions(503)
deltarvsptmin.GetYaxis().SetTitle('\Delta R')
deltarvsptmin.GetXaxis().SetTitle('jets p_{T}[GeV]')
deltarvsptmin.SetTitle('Jets:~\Delta R~vs~p_{T} ')
deltarvsptmin.SetStats(0)
deltarvsptmin.SetLineColor(3)
deltarvsptmin.Draw()
legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("jetdphi","p_{T}","f")
#legend.Draw()
c17.Print("2D_jet_deltar_vs_pt_min.png")

#//////////////////////////////////////////////////

#Probability
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

#/////////////////////////////////////////////////////

c18 = TCanvas('cHisto','cHisto',600,600,700,500)
pjeta5.Draw("AP")
pjeta5.SetTitle('Probabilty; jets \eta; Probability')
pjeta5.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pjeta5.SetLineColor(3)
pjeta1_5.Draw("SAME")
pjeta1_5.SetMarkerColor(2)
pjeta1_5.SetLineColor(2)
legend = TLegend(0.4,0.8,0.7,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("pjeta5","muon PAT p_{T}>5GeV","f")
legend.AddEntry("pjeta1_5","muon L1T p_{T}>5GeV","f")
legend.Draw()
c18.Print("probability_pat_vs_l1t.png")

#/////////////////////////////////////////////////////

c19 = TCanvas('cHisto','cHisto',600,600,700,500)
pjeta5.Draw("AP")
pjeta5.SetTitle('Probabilty; jets \eta; Probability')
pjeta5.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pjeta5.SetLineColor(3)
pjeta10.Draw("SAME")
pjeta10.SetMarkerColor(2)
pjeta10.SetLineColor(2)
pjeta20.Draw("SAME")
legend = TLegend(0.4,0.8,0.7,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("pjeta5","muon p_{T}>5GeV","f")
legend.AddEntry("pjeta10","muon p_{T}>10GeV","f")
legend.AddEntry("pjeta20","muon p_{T}>24GeV","f")
legend.Draw()
c19.Print("probability_pat_cats.png")

#//////////////////////////////////////////////////////

c20 = TCanvas('cHisto','cHisto',600,600,700,500)
pjeta1_5.Draw("AP")
pjeta1_5.SetTitle('Probabilty; jets \eta; Probability')
pjeta1_5.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pjeta1_5.SetLineColor(3)
pjeta1_10.Draw("SAME")
pjeta1_10.SetMarkerColor(2)
pjeta1_10.SetLineColor(2)
pjeta1_20.Draw("SAME")
legend = TLegend(0.4,0.8,0.7,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("pjeta1_5","L1T muon p_{T}>5GeV","f")
legend.AddEntry("pjeta1_10","L1T muon p_{T}>10GeV","f")
legend.AddEntry("pjeta1_20","L1T muon p_{T}>24GeV","f")
legend.Draw()
c20.Print("probability_l1t_cats.png")

#//////////////////////////////////////////////////////

c21 = TCanvas('cHisto','cHisto',600,600,700,500)
pmeta50.Draw("AP")
pmeta50.SetTitle('Probabilty; PAT muon \eta; Probability')
pmeta50.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pmeta50.SetLineColor(3)
pmeta100.Draw("SAME")
pmeta100.SetMarkerColor(2)
pmeta100.SetLineColor(2)
pmeta200.Draw("SAME")
legend = TLegend(0.4,0.8,0.7,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("pmeta50","jets p_{T}>50GeV","f")
legend.AddEntry("pmeta100","jets p_{T}>100GeV","f")
legend.AddEntry("pmeta200","jets p_{T}>200GeV","f")
legend.Draw()
c21.Print("probability_muon_jet_cats.png")

#/////////////////////////////////////////////////////

c22 = TCanvas('cHisto','cHisto',600,600,700,500)
pmeta1_50.Draw("AP")
pmeta1_50.SetTitle('Probabilty; L1T muon \eta; Probability')
pmeta1_50.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pmeta1_50.SetLineColor(3)
pmeta1_100.Draw("SAME")
pmeta1_100.SetMarkerColor(2)
pmeta1_100.SetLineColor(2)
pmeta1_200.Draw("SAME")
legend = TLegend(0.4,0.8,0.7,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("pmeta1_50","jets p_{T}>50GeV","f")
legend.AddEntry("pmeta1_100","jets p_{T}>100GeV","f")
legend.AddEntry("pmeta1_200","jets p_{T}>200GeV","f")
legend.Draw()
c22.Print("probability_l1t_jet_cats.png")

#//////////////////////////////////////////////////////

c18 = TCanvas('cHisto','cHisto',600,600,700,500)
pml1t5.Draw("AP")
pml1t5.SetTitle('Probabilty; muon \eta; Probability')
pml1t5.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pml1t5.SetLineColor(3)
pl1tm5.Draw("SAME")
pl1tm5.SetMarkerColor(2)
pl1tm5.SetLineColor(2)
legend = TLegend(0.4,0.8,0.7,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("pml1t5","L1T to PAT p_{T}>5GeV","f")
legend.AddEntry("pl1tm5","PAT to L1T p_{T}>5GeV","f")
legend.Draw()
c18.Print("probability_pat_l1t.png")
'''
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


histo = gROOT.FindObject('masan2m0_120')
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
histo.GetXaxis().SetTitle('\sqrt{s}[GeV]')
#histo.GetXaxis().SetTitle('vertexz')
histo.SetTitle('\sqrt{s}~of~two~muon')
histo.SetStats(0)

pmeta50.Draw("AP")
pmeta50.SetTitle('Probabilty; muon \eta; Probability')
pmeta50.SetMarkerColor(3)
#peta.SetMarkerStyle(4)
pmeta50.SetLineColor(3)
pmeta100.Draw("SAME")
pmeta100.SetMarkerColor(2)
pmeta100.SetLineColor(2)
pmeta200.Draw("SAME")

histo.Draw()
#histo.SetLineColor(2)
#histo1.Draw("SAME")
#histo1.SetLineColor(3)
#histo2.Draw("SAME")
#histo2.SetLineColor(4)
#histo3.Draw("SAME")
#histo3.SetLineColor(5)

c1.SetLeftMargin(0.2)
c1.SetGridy()
c1.SetGridx()

#c1.SetLogy()
#c1.SetLogx()
#c1.Update()

legend = TLegend(0.8,0.8,0.9,0.9)
legend.SetTextSize(0.03)
legend.AddEntry("masan2m0_120","\sqrt{s}","f")
#legend.AddEntry("jeteta50","\eta - p_{t}>50GeV","f")
#legend.AddEntry("jeteta100","\eta - p_{t}>100GeV","f")
#legend.AddEntry("jeteta200","\eta - p_{t}>200GeV","f")
#legend.AddEntry("pmeta50","jet pt>50","f")
#legend.AddEntry("pmeta100","jet pt>100","f")
#legend.AddEntry("pmeta200","jet pt>200","f")
#legend.AddEntry("peta3","PAT muon pt>10","f")
#legend.AddEntry("peta4","PAT muon pt>20","f")
#legend.AddEntry("histo","P_{t}>5GeV","f")
#legend.AddEntry("histo1","P_{t}>10GeV","f")
#legend.AddEntry("histo2","P_{t}>20GeV","f")
#legend.AddEntry("histo3","\eta - pt>200GeV","f")
#legend.AddEntry("histo1","l1phi","f")
#legend.AddEntry("histo2","p","f")
legend.Draw()
c1.Print("muon_2muon_masa_n0120.png")
'''
input('press enter to exit')
