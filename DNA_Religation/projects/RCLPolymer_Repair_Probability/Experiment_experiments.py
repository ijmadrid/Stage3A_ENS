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
import matplotlib.pyplot as plt
import scipy.stats as sts
from time import strftime, time

import multiprocess as mp

###########################################################################

if __name__ == '__main__':
    
    start = time()
    print("NAIVE----------------")
    P = []
    for i in range(1000):
        p0 = RCLPolymer(100, 3, 0.2, 5)
        p0.step(1000,0.01,0.008)
        P.append(p0)
    print("For time : ",time()-start)
    
    start = time()
    print("PARALEL--------------")
    def sample():
        p0 = RCLPolymer(100, 3, 0.2, 5)
        p0.step(1000,0.01,0.008)
        return p0
    
    pool = mp.Pool(8)
    
    res = [pool.apply_async(sample) for _ in range(1000)]
    res = [h.get() for h in res]
    print("Paralel time : ",time()-start)
    
    
else:
    p0 = RCLPolymer(100, 3, 0.2, 5)
    params = dict(numRealisations = 500,
                  dt_relax = 0.5,
                  diffusionConstant = 0.008,
                  numSteps = 10,
                  dt = 0.1)
    results = {}
    experiment = Experiment(p0, results, params)
    
#    # Test with custom experiment
#    def customExperiment(exp):
#        print(exp.stringTest)
#        exp.addResults('polymer_same',True)
#    
#    results = {}
#    params = dict(stringTest = 'HelloWorld!')
#    exp = Experiment(p0, results, params, customExperiment)
#    
#    exp = Experiment(p0, results, params, 10)