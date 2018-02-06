import numpy as np


class RealisticSimClusterCalibration:

    def __init__(self):
        # source https://github.com/cms-sw/cmssw/blob/master/RecoParticleFlow/PFClusterProducer/python/particleFlowRealisticSimClusterHGCCalibrations_cfi.py
        self.minEtaCorrection = 1.4
        self.maxEtaCorrection = 3.0
        self.hadronCorrections = np.array([1.24, 1.24, 1.24, 1.23, 1.24, 1.25, 1.29, 1.29])
        self.egammaCorrections = np.array([1.00, 1.00, 1.01, 1.01, 1.02, 1.03, 1.04, 1.04])
        self.hadronBins = np.linspace(self.minEtaCorrection, self.maxEtaCorrection, len(self.hadronCorrections)+1)
        self.egammaBins = np.linspace(self.minEtaCorrection, self.maxEtaCorrection, len(self.hadronCorrections)+1)

    def getEnergyCorrection(self, eta, isHadron=True):
        if isHadron:
            return self.hadronCorrections[np.digitize(abs(eta), self.hadronBins)-1]
        else:
            return self.egammaCorrections[np.digitize(abs(eta), self.egammaBins)-1]
