import sys
import math
from ROOT import *

print ("Hello ROOT")
fileName = "histos.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();

def calculate_histogram_mean(histogram):
    num_bins = histogram.GetNbinsX()
    sum_value = 0.0
    total_weight = 0.0

    for i in range(1, num_bins + 1):
        bin_content = histogram.GetBinContent(i)
        bin_center = histogram.GetBinCenter(i)
        bin_width = histogram.GetBinWidth(i)

        sum_value += bin_content * bin_center
        total_weight += bin_content * bin_width

    if total_weight > 0.0:
        return sum_value / total_weight

    return 0.0


jetcount50 = gROOT.FindObject('jetcount50')
jetcount100 = gROOT.FindObject('jetcount100')
jetcount200 = gROOT.FindObject('jetcount200')


jet50 = calculate_histogram_mean(jetcount50)
jet100 = calculate_histogram_mean(jetcount100)
jet200 = calculate_histogram_mean(jetcount200)

print('Srednia liczba jetow 50: ', jet50)
print('Srednia liczba jetow 100: ', jet100)
print('Srednia liczba jetow 200: ', jet200)

muoncount5 = gROOT.FindObject('muoncount5')
muoncount10 = gROOT.FindObject('muoncount10')
muoncount20 = gROOT.FindObject('muoncount20')

muon5 = calculate_histogram_mean(muoncount5)
muon10 = calculate_histogram_mean(muoncount10)
muon20 = calculate_histogram_mean(muoncount20)

print('Srednia liczba muonow 5: ', muon5)
print('Srednia liczba muonow 10: ',muon10)
print('Srednia liczba muonow 200: ', muon20)

