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
#from modules.Experiment import Experiment
from modules.Polymers import RCLPolymer

#import numpy as np
#import matplotlib.pyplot as plt
#import scipy.stats as sts
from time import time

#from joblib import Parallel, delayed
#import multiprocessing as mp
import ipyparallel as ipp

def sample(HowManySteps):
    p0 = RCLPolymer(100, 3, 0.2, 5)
    p0.step(HowManySteps,0.01,0.008)
    return p0

print("Something is in here")

###########################################################################

if __name__ == '__main__':
    
    steps = range(900,1000)
    start = time()
    print("NAIVE----------------")
    P = []
    for i in steps:
        p0 = RCLPolymer(100, 3, 0.2, 5)
        p0.step(i,0.01,0.008)
        P.append(p0)
    print("For time : ",time()-start)
    
    start = time()
    print("PARALLEL--------------")
#    res = Parallel(n_jobs=2)(delayed(sample)(i) for i in steps)
    pool = mp.Pool(4) 
    res = pool.map(sample, steps)
    print("Paralel time : ",time()-start)
    
    
else:
    
    print(__name__)
    print("Something is happening here")
#    p0 = RCLPolymer(100, 3, 0.2, 5)
#    params = dict(numRealisations = 500,
#                  dt_relax = 0.5,
#                  diffusionConstant = 0.008,
#                  numSteps = 10,
#                  dt = 0.1)
#    results = {}
#    experiment = Experiment(p0, results, params)
    
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