#!/bin/env python
import optparse
from NtupleDataFormat import HGCalNtuple
import hgcalHelpers
import numpy as np
import pandas as pd
from RecHitCalibration import RecHitCalibration as rhc
import RealisticSimClusterCalibration
maxlayer = 52


def getPFClusters(pfClusters, recHits, skipLayers=[]):
    """
    get the skimmed pfClusters.
    returns a dataframe containing 4-vectors
    """

    newPFClusters = []

    for idx, pfClus in pfClusters.iterrows():
        print pfClus
        newPFCluster = pfClus
        if abs(pfClus.eta) < 3:
            newPFCluster = getSinglePFCluster(pfClus, recHits, skipLayers)
        newPFClusters.append(newPFCluster)

    newPFClustersDF = pd.DataFrame(newPFClusters, columns=['pt', 'eta', 'phi', 'energy', 'correctedEnergy'])

    return newPFClustersDF


def getSinglePFCluster(pfClus, recHits, skipLayers):

    energySum = 0
    pTSum = 0
    skippedPreviousLayer = False
    layerEnergyCorrection = 1.
    clusterCalib = RealisticSimClusterCalibration.RealisticSimClusterCalibration()
    clusterEnergyCorrection = clusterCalib.getEnergyCorrection(pfClus['eta'])
    if skipLayers:
        dEdX_weights = rhc().dEdX_weights
    # find associated RecHits
    recHitIndices = hgcalHelpers.getRecHitsAssociatedToSimCluster(pfClus, recHits)
    associatedRecHits = recHits.iloc[recHitIndices]
    print type(associatedRecHits)
    print pfClus[["fractions"]]
    print pfClus[["fractions"]].values
    # df.at['C', 'x'] = 10
    associatedRecHits['fracs'] = associatedRecHits['index']
    for i, val in enumerate(associatedRecHits):
        print pfClus[["fractions"]][i]
        associatedRecHits.at[val["index"], "fracs"] = pfClus[["fractions"]][i]
    associatedRecHits["frac"] = pfClus[["fractions"]].values
    print associatedRecHits
    associatedRecHits = associatedRecHits.assign(fraction=pd.Series(pfClus[["fractions"]].values))
    print associatedRecHits
    associatedRecHits.loc[:, 'fraction'] = pd.Series(pfClus[["fractions"]], index=associatedRecHits.index)
    print associatedRecHits
    associatedRecHits['fraction'] = pd.Series(pfClus[["fractions"]], index=associatedRecHits.index)
    print associatedRecHits
    sys.exit(1)
    associatedRecHits['fractions']
    # now find layer clusters within the multicluster/track cone
    # maybe it's good to do this per layer to save some computing time
    for layer in range(1, maxlayer+1):
        if skippedPreviousLayer:
            if layer > 1:
                layerEnergyCorrection = (dEdX_weights[layer-1] + dEdX_weights[layer]) / dEdX_weights[layer]
            skippedPreviousLayer = False
        else:
            layerEnergyCorrection = 1.
        if layer in skipLayers:
            skippedPreviousLayer = True
            continue
        # require same layer
        selectedRecHits = associatedRecHits[(associatedRecHits.layer == layer)]
        print 'pfClus', type(pfClus[["fractions"]]), pfClus[["fractions"]]
        print 'selectedRecHits', type(selectedRecHits[["energy"]]), selectedRecHits[["energy"]]
        energySum += (selectedRecHits[["energy"]]*pfClus[["fractions"]]).sum()[0]*layerEnergyCorrection
        pTSum += selectedRecHits[["pt"]].sum()[0]*layerEnergyCorrection

    correctedEnergy = energySum*clusterEnergyCorrection
    newPFCluster = [pTSum, pfClus['eta'], pfClus['phi'], energySum, correctedEnergy]
    print newPFCluster
    return newPFCluster


def getCollections(event):
    """
    get the collections to be used
    need pfClusters, recHits
    """
    pfClusters = event.getDataFrame(prefix="pfcluster")
    recHits = event.getDataFrame(prefix="rechit")
    return pfClusters, recHits


def main():
    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-f', '--files', dest='fileString', type='string',  default='root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-PU200FEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-PU200FEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_2.root', help='comma-separated file list')
    parser.add_option('-s', '--skipLayers', dest='skipLayers', type='string',  default='', help='comma-separated list of layers to skip')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    print "files:", opt.fileString
    print "skipLayers:", opt.skipLayers

    fileList = opt.fileString.split(",")
    skipLayers = []
    if opt.skipLayers:
        skipLayers = [int(i) for i in opt.skipLayers.split(",")]

    for fileName in fileList:
        ntuple = HGCalNtuple(opt.fileString)

        for event in ntuple:
            if (event.entry() > 11):
                break
            # get collections
            pfClusters, recHits = getCollections(event)
            newPFClusters = getPFClusters(pfClusters, recHits, skipLayers=skipLayers)
            print newPFClusters


if __name__ == '__main__':
    main()
