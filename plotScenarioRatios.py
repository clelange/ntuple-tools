import ROOT
# import os
# import sys
# import optparse
import math
# import hgcalHistHelpers
import pickle
from array import array

ROOT.gROOT.SetBatch(1)


def getData(path, pkl_file):
    dataDict = None
    with open(path + "/" + pkl_file, "rb") as inFile:
        dataDict = pickle.load(inFile)
    return dataDict


def main():
    # pkl_file = "211_PF_manual_noPU_PF_manual_PU200.pkl"
    # nominal = "manual"
    # scenarios = {}
    # scenarios["manual_18_9_9"] = "18-9-9"
    # scenarios["manual_24_11_11"] = "24-11-11"
    # scenarios["manual_skipFH"] = "32-6-12"
    pkl_file = "211_Mega_noPU_Mega_PU200_Mega_noPU_PU200.pkl"
    nominal = "nominal"
    scenarios = {}
    scenarios["skip_18_9_9"] = "18-9-9"
    scenarios["skip_24_11_11"] = "24-11-11"
    scenarios["dropFH"] = "32-6-12"

    # first get nominal
    nominalDict = getData(nominal, pkl_file)
    # then loop over scenarios
    for scen, scenName in scenarios.iteritems():
        # create plain ratio and squared difference
        thisDict = getData(scen, pkl_file)
        ratioGraphs = {}
        diffGraphs = {}
        maxRatio = 1.5
        minRatio = 1.
        maxDiff = 40
        minDiff = 0
        for puScen, values in thisDict.iteritems():
            print puScen
            counter = 0
            x_axis = []
            ratios = []
            squaredDiff = []
            for x, y in values:
                x_axis.append(x)
                print x, y, nominalDict[puScen][counter][1]
                ratios.append(y/nominalDict[puScen][counter][1])
                squaredDiff.append(math.sqrt(abs(y*y - nominalDict[puScen][counter][1]*nominalDict[puScen][counter][1])))
                counter += 1
            print x_axis
            print ratios
            print squaredDiff
            ratioGraphs[puScen] = ROOT.TGraph(len(x_axis), array('d', x_axis), array('d', ratios))
            diffGraphs[puScen] = ROOT.TGraph(len(x_axis), array('d', x_axis), array('d', squaredDiff))
            # ratioGraphs[puScen].SetMarkerSize(20)
            ratioGraphs[puScen].SetMaximum(maxRatio)
            ratioGraphs[puScen].SetMinimum(minRatio)
            # diffGraphs[puScen].SetMarkerSize(20)
            diffGraphs[puScen].SetMaximum(maxDiff)
            diffGraphs[puScen].SetMinimum(minDiff)
            ratioGraphs[puScen].GetXaxis().SetTitle("p_{T} [GeV]")
            diffGraphs[puScen].GetXaxis().SetTitle("p_{T} [GeV]")
            ratioGraphs[puScen].GetYaxis().SetTitle("#sigma_{modified} / #sigma_{nominal}")
            diffGraphs[puScen].GetYaxis().SetTitle("#sqrt{ #Delta_{#sigma}^{2} } [%]")
            ratioGraphs[puScen].SetTitle(scenName)
            diffGraphs[puScen].SetTitle(scenName)
            if puScen.find("200") >= 0:
                ratioGraphs[puScen].SetLineColor(ROOT.kRed)
                diffGraphs[puScen].SetLineColor(ROOT.kRed)
                ratioGraphs[puScen].SetMarkerColor(ROOT.kRed)
                diffGraphs[puScen].SetMarkerColor(ROOT.kRed)
                ratioGraphs[puScen].SetMarkerStyle(21)
                diffGraphs[puScen].SetMarkerStyle(21)
            else:
                ratioGraphs[puScen].SetMarkerStyle(22)
                diffGraphs[puScen].SetMarkerStyle(22)
        c = ROOT.TCanvas("c", "", 800, 600)
        c.Draw()
        leg = ROOT.TLegend(0.15, 0.70, 0.82, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        firstGraph = True
        for puScen, values in thisDict.iteritems():
            leg.AddEntry(ratioGraphs[puScen], puScen, "pl")
            if firstGraph:
                ratioGraphs[puScen].Draw("apl")
                firstGraph = False
            else:
                ratioGraphs[puScen].Draw("pl same")
        leg.Draw()
        c.SaveAs(scen + "_ratio.pdf")
        c.SaveAs(scen + "_ratio.png")
        firstGraph = True
        for puScen, values in thisDict.iteritems():
            if firstGraph:
                diffGraphs[puScen].Draw("apl")
                firstGraph = False
            else:
                diffGraphs[puScen].Draw("pl same")
        leg.Draw()
        c.SaveAs(scen + "_diff.pdf")
        c.SaveAs(scen + "_diff.png")







if __name__ == '__main__':
    main()
