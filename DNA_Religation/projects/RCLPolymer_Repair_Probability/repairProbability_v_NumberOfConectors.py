# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:48:01 2018

@author: ignacio
"""

#############################################################
############ MODULES AND PACKAGES ###########################
#############################################################

## To import the Polymer and Simulation modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt

from modules.Polymers import RCLPolymer
from modules.Experiment import Experiment

from time import strftime, time

from pathos.multiprocessing import Pool

#############################################################
############ PARAMETERS #####################################
#############################################################

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

numRealisations = 500
params['numRealisations'] = numRealisations


test_genomic_distances = [3,5,10]
x_Nc = np.arange(3,20,2)

#############################################################
############ FUNCTIONS ######################################
#############################################################    

def makeExperiment(numrealisations):
    params['numRealisations'] = numrealisations
    np.random.seed()
    return Experiment(p0, results, params,"EncounterSimulation") 


def proba_vs_conectorsNumber(test_genomic_distances,x_Nc,keepCL,errorbars=False):
    probas = []
    demiCI = []

    output = np.zeros((len(test_genomic_distances),3,len(x_Nc)))
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    plt.figure()
    
    for i, genomicDistance in enumerate(test_genomic_distances):
        probas = []
        demiCI = []
        
        for Nc in x_Nc:
            results = {}
            p0 = RCLPolymer(100, 3, 0.2, Nc, keepCL)
            params['genomicDistance'] = genomicDistance
            mc = Experiment(p0, results, params,"EncounterSimulation") 
            probas.append(mc.results['repair_probability_CI'][0])
            demiCI.append(mc.results['repair_probability_CI'][1])
        
        probas = np.array(probas)
        demiCI = np.array(demiCI)
        
        output[i][0] = x_Nc
        output[i][1] = probas
        output[i][2] = demiCI

        np.save('results/proba_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(100)+'monomers_'+
            str(numRealisations)+'iterations'+
            date+'.npy',output)
        
        if errorbars:
            plt.errorbar(x=x_Nc, y=probas, yerr=demiCI,
                     fmt='-o', label=r'$g = $ '+str(genomicDistance), capsize = 4)
        else:
            plt.plot(x_Nc,probas,'-o',label=r'$g = $ '+str(genomicDistance))

    np.save('results/proba_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(100)+'monomers_'+
            str(numRealisations)+'iterations'+
            date+'.npy',output)
    
    plt.legend()
    plt.show()


def measureProba(g):
    probas = []
    params['genomicDistance'] = g
    
    for Nc in x_Nc:            
        results = {}
        np.random.seed()
        p0 = RCLPolymer(100, 3, 0.2, Nc, False)
        mc = Experiment(p0, results, params,"EncounterSimulation") 
        probas.append(mc.results['repair_probability_CI'][0])
   
    return probas


if __name__ == '__main__':
    
    print("::::::: SERIAL EXPERIMENTS :::::::::::")
    
    verystart = time()
    
    proba_vs_conectorsNumber(test_genomic_distances,x_Nc,False,errorbars=True)

    print("Total duration : ", time()-verystart)
    
    print("::::::: PARALLEL EXPERIMENTS :::::::::")
    
    plt.figure()
    
    verystart = time()
    workers = len(test_genomic_distances)
    realisationsPerWorker = numRealisations
    pool = Pool(workers)
    
    res = pool.map(measureProba, test_genomic_distances)
    
    for i in range(len(res)):
        plt.plot(x_Nc,res[i],'-o',label=r'$g = $ '+str(test_genomic_distances[i]))
    
    plt.legend()
    plt.show()    
    
    print("Total duration : ", time()-verystart)