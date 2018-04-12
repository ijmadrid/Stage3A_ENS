# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 15:09:09 2018

@author: ignacio
"""

############################################################################
# PACKAGES AND MODULES

## To import the Polymer and Simulation modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from modules.Experiment import Experiment
from modules.Polymers import RCLPolymer

import numpy as np
from time import time

from pathos.multiprocessing import Pool


p0 = RCLPolymer(100, 3, 0.2, 3)
params = dict(dt_relax = 0.5,
              diffusionConstant = 0.008,
#              numSteps = 100,
              dt = 0.1,
              genomicDistance = 2,
              Nb = 2,
              waitingSteps = 200,
              numMaxSteps = 600,
              encounterDistance = 0.1)
results = {}


def makeExperiment(numrealisations):
    params['numRealisations'] = numrealisations
    np.random.seed()
    return Experiment(p0, results, params,"EncounterSimulation") 

print("Something is in here")

###########################################################################

if __name__ == '__main__':
    numrealisations = 1000
    print(":::::::::: PARALLEL EXPERIMENTS :::::::::::::::::::")
    start = time()
    workers = 10
    pool = Pool(workers)
    realisationsPerWorker = int(numrealisations/workers)
#    experiment = Experiment(p0, results, params)
    parallelExperiments = pool.map(makeExperiment, [realisationsPerWorker]*workers)
    print("Experiment duration : ", time()-start)
    
    print(":::::::::::: SERIAL EXPERIMENT ::::::::::::::::::::")
    start = time()
    serialExperiment = makeExperiment(numrealisations)
    print("Experiment duration : ", time()-start)