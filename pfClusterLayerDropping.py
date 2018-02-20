#!/bin/env python
import optparse
from NtupleDataFormat import HGCalNtuple
import hgcalHelpers
# import numpy as np
import pandas as pd
from RecHitCalibration import RecHitCalibration as rhc
import RealisticSimClusterCalibration
maxlayer = 52


def getPFClusters(genParticles, pfClusters, recHits, gun_type, GEN_engpt, pidSelected, skipLayers=[]):
    """
    get the skimmed pfClusters.
    returns a dataframe containing 4-vectors
    """

    newPFClusters = []

    # use genParticles with generated pT/energy that reach EE before converting
    selectedGen = genParticles[(abs(genParticles.pid) == pidSelected) & (genParticles.reachedEE > 0)]
    if gun_type == "pt":
        selectedGen = selectedGen[(selectedGen.pt >= GEN_engpt*.999)]
    else:
        selectedGen = selectedGen[(selectedGen.energy >= GEN_engpt*.999)]

    if pfClusters.shape[0] == 0:
        return pd.DataFrame(columns=['pt', 'eta', 'phi', 'energy', 'correctedEnergy'])

    # for the axis, take highest energy/pt object within dR = 0.1; 0.2 for <= 7 GeV guns
    bestTruthMatchedIndices = None
    if GEN_engpt <= 7.5:
        bestTruthMatchedIndices = hgcalHelpers.getHighestEnergyObjectIndex(selectedGen[['eta', 'phi']], pfClusters[['eta', 'phi']], pfClusters['energy'], 0.2)
    else:
        bestTruthMatchedIndices = hgcalHelpers.getHighestEnergyObjectIndex(selectedGen[['eta', 'phi']], pfClusters[['eta', 'phi']], pfClusters['energy'], 0.1)

    for idx, genPart in selectedGen.iterrows():
        if idx not in bestTruthMatchedIndices:
            # this is the case if there's no match between gen and the reconstructed object
            continue
        pfClus = pfClusters.iloc[[bestTruthMatchedIndices[idx]]]
        # print "### matched pfClus"
        # print pfClus
        # print "### end matched pfClus"

        newPFCluster = pfClus
        # import pdb
        # pdb.set_trace()
        if abs(pfClus['eta'].item()) < 3:
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
    # print type(associatedRecHits)
    # print pfClus[["fractions"]]
    # print pfClus[["fractions"]].values
    # df.at['C', 'x'] = 10
    # import pdb
    # pdb.set_trace()
    # associatedRecHits['fracs'] = associatedRecHits['index']
    # for i, val in enumerate(associatedRecHits):
    #     print pfClus[["fractions"]][i]
    #     associatedRecHits.at[val["index"], "fracs"] = pfClus[["fractions"]][i]
    # associatedRecHits["frac"] = pfClus[["fractions"]].values
    # print associatedRecHits
    # associatedRecHits = associatedRecHits.assign(fraction=pd.Series(pfClus[["fractions"]].values))
    # print associatedRecHits
    # associatedRecHits.loc[:, 'fraction'] = pd.Series(pfClus[["fractions"]], index=associatedRecHits.index)
    # print associatedRecHits
    # associatedRecHits['fraction'] = pd.Series(pfClus[["fractions"]], index=associatedRecHits.index)
    """ all of the following commands give an error/warning, but seem to work:
    /cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/py2-pippkgs_depscipy/3.0-ghjeda6/lib/python2.7/site-packages/pandas/core/indexing.py:337: SettingWithCopyWarning:
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead

    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      self.obj[key] = _infer_fill_value(value)
    /cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/py2-pippkgs_depscipy/3.0-ghjeda6/lib/python2.7/site-packages/pandas/core/indexing.py:517: SettingWithCopyWarning:
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead

    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      self.obj[item] = s"""
    # disable warning
    pd.set_option('mode.chained_assignment', None)
    # associatedRecHits.loc[:, ('fraction')] = pd.Series(pfClus["fractions"], index=associatedRecHits.index)  # this worked
    # associatedRecHits.loc[:, 'fraction'] = pd.Series(pfClus["fractions"], index=associatedRecHits.index)  # this also works apparently
    assert(len(associatedRecHits) == len(*pfClus["fractions"])), "Cannot assign fractions!"
    associatedRecHits.loc[:, 'fraction'] = tuple(*pfClus["fractions"])  # also works
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
            # skip this layer
            skippedPreviousLayer = True
            continue
        # require same layer
        selectedRecHits = associatedRecHits[(associatedRecHits.layer == layer)]
        energySum += (selectedRecHits["energy"]*selectedRecHits["fraction"]).sum()*layerEnergyCorrection
        pTSum += (selectedRecHits["pt"]*selectedRecHits["fraction"]).sum()*layerEnergyCorrection

    correctedEnergy = energySum*clusterEnergyCorrection
    newPFCluster = [pTSum, pfClus['eta'].values[0], pfClus['phi'].values[0], energySum, correctedEnergy[0]]
    # print "############### newPFCluster ###############"
    # print newPFCluster
    # print "############### end newPFCluster ###############"
    return newPFCluster


def getCollections(event):
    """
    get the collections to be used
    need genParticles, pfClusters, recHits
    """
    genParticles = event.getDataFrame(prefix="genpart")
    pfClusters = event.getDataFrame(prefix="pfcluster")
    recHits = event.getDataFrame(prefix="rechit")
    return genParticles, pfClusters, recHits


def main():
    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-f', '--files', dest='fileString', type='string',  default='root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-PU200FEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-PU200FEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_2.root', help='comma-separated file list')
    parser.add_option('-t', '--gunType', dest='gunType', type='string',  default='pt', help='pt or e')
    parser.add_option('-p', '--pid', dest='pid', type='int',  default=211, help='pdgId int')
    parser.add_option('-g', '--genValue', dest='genValue', type='float',  default=50, help='generated pT or energy')
    parser.add_option('-s', '--skipLayers', dest='skipLayers', type='string',  default='', help='comma-separated list of layers to skip')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    print "files:", opt.fileString
    print "gunType:", opt.gunType
    print "pid:", opt.pid
    print "GEN_engpt:", opt.genValue
    print "skipLayers:", opt.skipLayers

    fileList = opt.fileString.split(",")
    gun_type = opt.gunType
    pidSelected = opt.pid
    GEN_engpt = opt.genValue
    skipLayers = []
    if opt.skipLayers:
        skipLayers = [int(i) for i in opt.skipLayers.split(",")]

    for fileName in fileList:
        ntuple = HGCalNtuple(opt.fileString)

        for event in ntuple:
            if (event.entry() > 11):
                break
            # get collections
            genParticles, pfClusters, recHits = getCollections(event)
            newPFClusters = getPFClusters(genParticles, pfClusters, recHits, gun_type, GEN_engpt, pidSelected, skipLayers=skipLayers)
            print newPFClusters


if __name__ == '__main__':
    main()
