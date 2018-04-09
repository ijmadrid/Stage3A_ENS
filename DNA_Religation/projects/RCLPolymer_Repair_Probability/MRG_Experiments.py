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
numRealisations             = 1000,  # num of realisations per set 
maxIterationsPerExperiment  = 400,  # if there is no encounter before
dt                          = 0.01, # time step after relaxation time
dt_relax                    = 0.05, # time step until relaxation time
waitingSteps                = 100, # steps after the DBSs

# Test parameters
test_connectorsNumber       = np.arange(2,20,dtype=int),
test_genomic_distances      = 3,

# Plot options
errorbars                   = False 
)

############################################################################


############################################################################
# FUNCTIONS

def MRG_vs_CLNumber(dimension,monomerNumber,b,diffusionConstant,encounterDistance,keepCL,
                    numRealisations,maxIterationsPerExperiment,dt,dt_relax,waitingSteps,
                    test_connectorsNumber,test_genomic_distances, errorbars):
    """
    Main function
    """
    
    # Output file
    output = np.zeros((len(test_genomic_distances),3,len(test_connectorsNumber)))
    
    date = strftime("%Y_%m_%d_%H_%M")
    
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

def MSRG_vs_CLremoval(dimension,monomerNumber,b,diffusionConstant,encounterDistance,
                    numRealisations,maxIterationsPerExperiment,dt,dt_relax,waitingSteps,
                    test_genomic_distances, errorbars,**kwargs):
    """
    Main function
    """
    
    # Output file
    genomic_distance = test_genomic_distances
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    plt.figure()
    plt.xlabel('Number of random cross-links')
    plt.ylabel(r'Mean square radius of gyration ($\mu^2$)')
       
    for keepCL in [True, False]:
        MSRG = []
        demiCI = []
        msrg = 100
        Nc_max = 1
        while(msrg > b + 0.1*b):
            Nc_max += 1
            p0 = RCLPolymer(monomerNumber, dimension, b, Nc_max, keepCL)
            mc = EncounterSimulation(dt, diffusionConstant, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, genomic_distance, encounterDistance, waitingSteps)
            mc.run()
            msrg = mc.get_msrg()
            MSRG.append(msrg)
            demiCI.append(1.96*np.std(mc.msrg)/np.sqrt(len(mc.msrg)))
        
        MSRG = np.array(MSRG)
        demiCI = np.array(demiCI)
        
        output = np.zeros((3,Nc_max-1))
        output[0] = np.arange(2,Nc_max+1)
        output[1] = MSRG
        output[2] = demiCI

        np.save('results/MSRG_vs_RCLremoval__'+
            'keepCL_'+str(keepCL)+
            str(monomerNumber)+'monomers_'+
            str(genomic_distance)+'genomicDistance_'+
            str(numRealisations)+'iterations'+
            date+'.npy',output)
        
        if keepCL:
            labelExp = "Keeping RCLs"
        else:
            labelExp = "Removing RCLs"
            
        if errorbars:
            plt.errorbar(x=np.arange(2,Nc_max+1), y=MSRG, yerr=demiCI,
                     fmt='-o', label=labelExp, capsize = 4)
        else:
            plt.plot(np.arange(2,Nc_max+1),MSRG,'-o',label=labelExp)
            
    plt.legend()
    plt.show()
    

############################################################################

if __name__ == '__main__':
    MSRG_vs_CLremoval(**params)