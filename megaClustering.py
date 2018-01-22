# investigate shower development based on RecHits and SimClusters
# import ROOT
# import os
import optparse
# from array import array
# from HGCalImagingAlgo import recHitAboveThreshold
from NtupleDataFormat import HGCalNtuple
# from GeoUtils import GeoUtil
# import math
import hgcalHelpers
# import hgcalHistHelpers
import numpy as np
import pandas as pd
from itertools import repeat
maxlayer = 52
energyWeights = list(repeat(1.02, 28)) + list(repeat(0.86, 12)) + list(repeat(1.12, 12))


def getConeRadius(frontRadius, backRadius, z, maxval=9999.):
    depthTerm = backRadius * (abs(z)-320.7)/(407.8-320.7)
    val = frontRadius + depthTerm
    if val > maxval:
        return maxval
    return val


def pileupSubtraction(matchedObject, selectedLayerClusters, recHits, layer, energyRadius, frontRadius, backRadius):
    """For now, the code is the same as in getMegaClusters, but the phi coordinate is changed by pi."""

    pTSum = 0
    energySum = 0
    # take first layer cluster z value
    layer_z = selectedLayerClusters.head(1).z.item()
    # get multi cluster x and y coordinates
    matchedObject_phi = float(matchedObject.phi) - np.pi
    if (matchedObject_phi < -np.pi):
        matchedObject_phi += 2*np.pi
    multiClusPosDF = hgcalHelpers.convertToXY(matchedObject.eta, matchedObject_phi, layer_z)
    # calculate radius based on current layer's z position
    coneRadius = getConeRadius(frontRadius, backRadius, layer_z)
    # mind that we need only the first index since there is only one multiCluster
    layerClusterIndices = hgcalHelpers.getIndicesWithinRadius(multiClusPosDF[['x', 'y']], selectedLayerClusters[['x', 'y']], coneRadius)
    # now we need to recalculate the layer cluster energies using associated RecHits
    for layerClusterIndex in layerClusterIndices[0]:
        associatedRecHits = recHits.iloc[selectedLayerClusters.iloc[layerClusterIndex].rechits]
        # find maximum energy RecHit
        maxEnergyRecHitIndex = associatedRecHits['energy'].argmax()
        # considering only associated RecHits within a radius of energyRadius (6 cm)
        matchedRecHitIndices = hgcalHelpers.getIndicesWithinRadius(associatedRecHits.loc[[maxEnergyRecHitIndex]][['x', 'y']], associatedRecHits[['x', 'y']], energyRadius)[maxEnergyRecHitIndex]
        # sum up energies and pT
        selectedRecHits = associatedRecHits.iloc[matchedRecHitIndices]
        # correct energy by subdetector weights
        energySum += selectedRecHits[["energy"]].sum()[0]*energyWeights[layer-1]*1.38
        pTSum += selectedRecHits[["pt"]].sum()[0]*energyWeights[layer-1]*1.38

    return (energySum, pTSum)


def getMegaClusters(genParticles, axisCollection, axisCollectionName, layerClusters, recHits, gun_type, GEN_engpt, pidSelected, energyRadius=6, frontRadius=3, backRadius=8, doPileupSubtraction=True):
    """
    get the actual mega clusters.
    frontRadius: cone at front of EE
    backRadius: cone at back of FH (frontRadius to be added to it)
    returns a dataframe containing 4-vectors
    """

    useTracks = False
    if "track" in axisCollectionName:
        useTracks = True

    # use genParticles with generated pT/energy that reach EE before converting
    selectedGen = genParticles[(abs(genParticles.pid) == pidSelected) & (genParticles.reachedEE > 0)]
    if gun_type == "pt":
        selectedGen = selectedGen[(selectedGen.pt >= GEN_engpt*.999)]
    else:
        selectedGen = selectedGen[(selectedGen.energy >= GEN_engpt*.999)]
    # print selectedGen

    if axisCollection.shape[0] == 0:
        return pd.DataFrame(columns=['pt', 'eta', 'phi', 'energy'])
    bestTruthMatchedIndices = None
    maxPtOrEnergy = 'energy'
    if useTracks:
        maxPtOrEnergy = 'pt'
    # for the axis, take highest energy/pt object within dR = 0.1; 0.2 for <= 7 GeV guns
    if GEN_engpt <= 7.5:
        bestTruthMatchedIndices = hgcalHelpers.getHighestEnergyObjectIndex(selectedGen[['eta', 'phi']], axisCollection[['eta', 'phi']], axisCollection[maxPtOrEnergy], 0.2)
    else:
        bestTruthMatchedIndices = hgcalHelpers.getHighestEnergyObjectIndex(selectedGen[['eta', 'phi']], axisCollection[['eta', 'phi']], axisCollection[maxPtOrEnergy], 0.1)
    # print bestTruthMatchedIndices

    megaClusters = []

    for idx, genPart in selectedGen.iterrows():
        if idx not in bestTruthMatchedIndices:
            # this is the case if there's no match between gen and the reconstructed object
            continue
        matchedObject = None
        if useTracks:
            matchedObject = axisCollection.iloc[[bestTruthMatchedIndices[idx]]]
        else:
            matchedObject = axisCollection.iloc[[bestTruthMatchedIndices[idx]]]
        megaCluster = getSingleMegaCluster(matchedObject, layerClusters, recHits, useTracks, energyRadius, frontRadius, backRadius, doPileupSubtraction)
        megaClusters.append(megaCluster)

    megaClustersDF = pd.DataFrame(megaClusters, columns=['pt', 'eta', 'phi', 'energy'])

    return megaClustersDF


def getSingleMegaCluster(matchedObject, layerClusters, recHits, useTracks, energyRadius, frontRadius, backRadius, doPileupSubtraction):

    energySum = 0
    pTSum = 0
    # now find layer clusters within the multicluster/track cone
    # maybe it's good to do this per layer to save some computing time
    for layer in range(1, maxlayer+1):
        # match only in same detector side
        etaSide = float(matchedObject.eta)  # not sure why this needs to be cast
        selectedLayerClusters = layerClusters[(layerClusters.layer == layer) & (layerClusters.eta*etaSide > 0)]
        if selectedLayerClusters.shape[0] == 0:
            # continue if no layer clusters selected
            print "layer", layer, "no layer clusters selected"
            continue
        layer_z = 0
        posDF = None
        if useTracks:
            layer_z = matchedObject["posz"].item()[layer-1]
            # get x and y coordinates
            d = {'x': [matchedObject["posx"].item()[layer-1]], 'y': [matchedObject["posy"].item()[layer-1]]}
            posDF = pd.DataFrame(data=d)
        else:
            # take first layer cluster z value
            layer_z = selectedLayerClusters.head(1).z.item()
            # get multi cluster x and y coordinates
            posDF = hgcalHelpers.convertToXY(matchedObject.eta, matchedObject.phi, layer_z)
        # calculate radius based on current layer's z position
        coneRadius = getConeRadius(frontRadius, backRadius, layer_z)
        # mind that we need only the first index since there is only one multiCluster
        layerClusterIndices = hgcalHelpers.getIndicesWithinRadius(posDF[['x', 'y']], selectedLayerClusters[['x', 'y']], coneRadius)
        # now we need to recalculate the layer cluster energies using associated RecHits
        for layerClusterIndex in layerClusterIndices[0]:
            associatedRecHits = recHits.iloc[selectedLayerClusters.iloc[layerClusterIndex].rechits]
            # find maximum energy RecHit
            maxEnergyRecHitIndex = associatedRecHits['energy'].argmax()
            # considering only associated RecHits within a radius of energyRadius (6 cm)
            matchedRecHitIndices = hgcalHelpers.getIndicesWithinRadius(associatedRecHits.loc[[maxEnergyRecHitIndex]][['x', 'y']], associatedRecHits[['x', 'y']], energyRadius)[maxEnergyRecHitIndex]
            # sum up energies and pT
            selectedRecHits = associatedRecHits.iloc[matchedRecHitIndices]
            # correct energy by subdetector weights
            energySum += selectedRecHits[["energy"]].sum()[0]*energyWeights[layer-1]*1.38
            pTSum += selectedRecHits[["pt"]].sum()[0]*energyWeights[layer-1]*1.38
        if (doPileupSubtraction):
            (pu_energySum, pu_pTSum) = pileupSubtraction(matchedObject, selectedLayerClusters, recHits, layer, energyRadius, frontRadius, backRadius)
            energySum -= pu_energySum
            pTSum -= pu_pTSum

    # use as coordinates eta and phi of matched multi cluster
    megaCluster = [pTSum, matchedObject['eta'].item(), matchedObject['phi'].item(), energySum]
    return megaCluster


def getCollections(event, axisCollectionName):
    """
    get the collections to be fed to the mega clustering.
    need genParticles, axisCollectionName, layerClusters, recHits
    """
    genParticles = event.getDataFrame(prefix="genpart")
    axisCollection = event.getDataFrame(prefix=axisCollectionName)
    layerClusters = event.getDataFrame(prefix="cluster2d")
    recHits = event.getDataFrame(prefix="rechit")
    return genParticles, axisCollection, layerClusters, recHits


def main():
    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    # parser.add_option('', '--files', dest='fileString', type='string',  default='root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_1.root', help='comma-separated file list')
    parser.add_option('-f', '--files', dest='fileString', type='string',  default='root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-PU200FEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-PU200FEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_2.root', help='comma-separated file list')
    parser.add_option('-t', '--gunType', dest='gunType', type='string',  default='pt', help='pt or e')
    parser.add_option('-p', '--pid', dest='pid', type='int',  default=211, help='pdgId int')
    parser.add_option('-g', '--genValue', dest='genValue', type='float',  default=50, help='generated pT or energy')
    parser.add_option('-a', '--axisCollectionName', dest='axisCollectionName', type='string',  default='multiclus', help='object to use as clustering axis (multiclus/track)')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    print "files:", opt.fileString
    print "gunType:", opt.gunType
    print "pid:", opt.pid
    print "GEN_engpt:", opt.genValue
    print "axisCollectionName:", opt.axisCollectionName

    # set sample/tree - for photons
    gun_type = opt.gunType
    pidSelected = opt.pid
    GEN_engpt = opt.genValue
    axisCollectionName = opt.axisCollectionName

    fileList = opt.fileString.split(",")

    for fileName in fileList:
        ntuple = HGCalNtuple(opt.fileString)

        for event in ntuple:
            if (event.entry() > 11):
                break
            # get collections

            genParticles, axisCollection, layerClusters, recHits = getCollections(event, axisCollectionName)
            megaClusters = getMegaClusters(genParticles, axisCollection, axisCollectionName, layerClusters, recHits, gun_type, GEN_engpt, pidSelected)
            print megaClusters


if __name__ == '__main__':
    main()
