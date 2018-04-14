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
from modules.Experiment import Experiment

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import strftime, time

###########################################################################


############################################################################
# PARAMETERS

params = dict(
# Polymer parameters
dimension                   = 3,
monomerNumber               = 100,
b                           = 0.2,
keepCL                      = True, # keep CL in cleaved monomers?

# Simulation parameters
simulationParams = dict(
# Physicial parameters
Nb                          = 2,
diffusionConstant           = 0.08,
encounterDistance           = 0.1,
# Numerical parameters
numRealisations             = 10,  # num of realisations per set 
dt                          = 0.01, # time step after relaxation time
dt_relax                    = 0.01, # time step until relaxation time
waitingSteps                = 10, # steps after the DBSs
),

# Test parameters
test_connectorsNumber       = np.arange(2,100,5,dtype=int),
test_genomic_distances      = [10],  # Only one if the test is MSRG_vs_CLremoval
        
# Plot options
errorbars                   = True 
)


############################################################################


############################################################################
# FUNCTIONS

def MRG_vs_CLNumber(dimension,monomerNumber,b,keepCL,simulationParams,
                    test_connectorsNumber,test_genomic_distances, errorbars):
    """
    Main function
    """
    
    # Output file
    output = np.zeros((len(test_genomic_distances),3,len(test_connectorsNumber)))
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    plt.figure()
    rcParams.update({'axes.labelsize': 'xx-large'})
    plt.xlabel("Number of cross-links")
    plt.ylabel("Mean radius of gyration ($\mu$m)")
    
    for i, genomicDistance in enumerate(test_genomic_distances):
        MRG = []
        demiCI = []
        
        for Nc in test_connectorsNumber:
            p0 = RCLPolymer(monomerNumber, dimension, b, Nc, keepCL)
            results = {}
            simulationParams['genomicDistance'] = genomicDistance
            mc = Experiment(p0, results, simulationParams, "twoDSB") 
            MRG.append(mc.results['MSRG'])
            demiCI.append(mc.results['MSRG_95CI'])
        
        MRG = np.array(MRG)
        demiCI = np.array(demiCI)
        
        output[i][0] = test_connectorsNumber
        output[i][1] = np.sqrt(MRG)
        output[i][2] = demiCI

        np.save('results/MRG_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(monomerNumber)+'monomers_'+
            str(simulationParams['numRealisations'])+'iterations'+
            date+'.npy',output)
        
        if errorbars:
            plt.errorbar(x=test_connectorsNumber, y=np.sqrt(MRG), yerr=demiCI,
                     fmt='-o', label=r'$g = $ '+str(genomicDistance), capsize = 4)
        else:
            plt.plot(test_connectorsNumber,MRG,'-o',label=r'$g = $ '+str(genomicDistance))

    np.save('results/MRG_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(monomerNumber)+'monomers_'+
            str(simulationParams['numRealisations'])+'iterations'+
            date+'.npy',output)
    
    plt.legend()
    plt.show()
    
    return output

def MSRG_vs_CLremoval(dimension,monomerNumber,b,keepCL,simulationParams,
                    test_connectorsNumber, test_genomic_distances, errorbars):
    """
    Main function
    """
    
    # Output file
    genomic_distance = test_genomic_distances[-1]
    simulationParams['genomicDistance'] = genomic_distance
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    output = np.zeros((2,3,len(test_connectorsNumber)))
    
    plt.figure()
    rcParams.update({'axes.labelsize': 'xx-large'})
    plt.xlabel('Number of random cross-links')
    plt.ylabel(r'Mean radius of gyration ($\mu$m)')
       
    for i, keepCL in enumerate([True, False]):

        MRG = []
        demiCI = []
        
        for Nc in test_connectorsNumber:
            p0 = RCLPolymer(monomerNumber, dimension, b, Nc, keepCL)
            results = {}
            mc = Experiment(p0, results, simulationParams, "twoDSB") 
            MRG.append(mc.results['MSRG'])
            demiCI.append(mc.results['MSRG_95CI'])
        
        MRG = np.array(MRG)
        demiCI = np.array(demiCI)
        
        output[i][0] = test_connectorsNumber
        output[i][1] = np.sqrt(MRG)
        output[i][2] = demiCI
        
        if keepCL:
            labelExp = "Keeping RCLs"
        else:
            labelExp = "Removing RCLs"
            
        if errorbars:
            plt.errorbar(x=test_connectorsNumber, y=MRG, yerr=demiCI,
                     fmt='-o', label=labelExp, capsize = 4)
        else:
            plt.plot(test_connectorsNumber, MRG,'-o',label=labelExp)

    np.save('results/MSRG_vs_RCLremoval__'+
        str(monomerNumber)+'monomers_'+
        str(genomic_distance)+'genomicDistance_'+
        str(simulationParams['numRealisations'])+'iterations'+
        date+'.npy',output)
            
    plt.legend(fontsize='xx-large')
    plt.show()
    

############################################################################

if __name__ == '__main__':
    start = time()
    z = MSRG_vs_CLremoval(**params)
    print("Experiment duration : ", time()-start)   