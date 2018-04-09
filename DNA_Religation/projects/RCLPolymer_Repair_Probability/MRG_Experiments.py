# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:26:42 2018

@author: ignacio
"""

"""
RADIUS OF GYRATION TEST
Tests the mean radius of gyration against the number of random conectors
in polymers cleaved by two DSBs separed by different genomic distances
"""

############################################################################
# PACKAGES AND MODULES

## To import the Polymer and Simulation modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from modules.Simulation import EncounterSimulation
from modules.Polymers import RCLPolymer

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts
from time import strftime

###########################################################################


############################################################################
# PARAMETERS

params = dict(
# Polymer parameters
dimension                   = 3,
monomerNumber               = 100,
b                           = 0.2,

# Simulation parameters
# Physicial parameters
diffusionConstant           = 0.08,
encounterDistance           = 0.1,
keepCL                      = True, # keep CL in cleaved monomers?
# Numerical parameters
numRealisations             = 500,  # num of realisations per set 
maxIterationsPerExperiment  = 700,  # if there is no encounter before
dt                          = 0.01, # time step after relaxation time
dt_relax                    = 0.05, # time step until relaxation time
waitingSteps                = 1000, # steps after the DBSs

# Test parameters
test_connectorsNumber       = [],
test_genomic_distances      = [1,2,3,4,5]
)

############################################################################


############################################################################
# FUNCTIONS

def MRG_vs_CLNumber(**params):
    """
    Main function
    """
    
    # Output file
    output = np.zeros((len(test_genomic_distances),3,len(test_connectorsNumber)))
    
    date = strftime("%Y-%m-%d-%H:%M")
    
    plt.figure()
    
    for i, genomicDistance in enumerate(test_genomic_distances):
        MRG = []
        demiCI = []
        
        for Nc in test_connectorsNumber:
            p0 = RCLPolymer(monomerNumber, dimension, b, Nc, keepCL)
            mc = EncounterSimulation(dt, diffusionConstant, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, genomicDistance, encounterDistance, waitingSteps)
            mc.run()
            MRG.append(mc.get_msrg())
            demiCI.append(1.96*np.std(mc.msrg)/np.sqrt(len(mc.msrg)))
        
        MRG = np.array(MRG)
        demiCI = np.array(demiCI)
        
        output[i][0] = test_connectorsNumber
        output[i][1] = MRG
        output[i][2] = demiCI

        np.save('results/MRG_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(monomerNumber)+'monomers_'+
            str(numRealisations)+'iterations'+
            date+'.npy',output)
        
        if errorbars:
            plt.errorbar(x=test_connectorsNumber, y=MRG, yerr=demiCI,
                     fmt='-o', label=r'$g = $ '+str(genomicDistance), capsize = 4)
        else:
            plt.plot(test_connectorsNumber,MRG,'-o',label=r'$g = $ '+str(genomicDistance))

    np.save('results/MRG_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(monomerNumber)+'monomers_'+
            str(numRealisations)+'iterations'+
            date+'.npy',output)
    
    plt.legend()
    plt.show()


############################################################################

if __name__ == '__main__':
    MRG_vs_CLNumber(**params)