# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 17:44:03 2018

@author: ignacio
"""


############################################################################
# PACKAGES AND MODULES #####################################################
############################################################################

## To import the Polymer and Experiments modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from modules.Polymers import RCLPolymer
from modules.Experiment import Experiment
from modules.Forces import ExcludedVolume
#from modules.Forces import ExcludedVolume, LocalExcludedVolume

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import strftime
import scipy.stats as sts
#import pickle
import csv

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

NcMatrix = np.ones((3,3),dtype=int)*0
NcMatrix[0,0] = 10
NcMatrix[1,1] = 10
NcMatrix[2,2] = 10
                     
polymerParams = dict(numMonomers = np.array([100,100,100]), #np.array([40,40]), # TODO (.., ..., ...)
                     dim         = 3,
                     b           = 0.2,
                     Nc          = NcMatrix,
                     keepCL      = False
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 500, 
                        dt                = 0.01,
                        dt_relax          = 0.01,
                        numSteps          = 5000,
                        excludedVolumeCutOff = 0.3,
                        waitingSteps = 250,
#                        numMaxSteps = 10000,
                        encounterDistance = 0.1,
#                        genomicDistance = 4,
#                        Nb = 2
#                        selectedSubDomain = 0
                        )

#test_distances = np.arange(1,30,3,dtype = int)
#
#gmax = 10
#gStep = 2
#test_epsilons = [0.1]
#
x_Nc = np.arange(2,7,1)
#
errorbars = True

############################################################################
############################################################################



def DSBclustering(polymerParams,simulationParams):
    """ 
    Tracks the breaks positions assuming they are persistent 
    """
    
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_interbreakdistances' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        rcParams.update({'legend.fontsize': 'large'})
        plt.xlabel('Time (sec)')
        plt.ylabel('Distance between A1 and A2') #(r'$\mathbb{P}$(Repair)')
        
        first_time = True
        for i, VolumeExclusion in enumerate([False]):
            
            if VolumeExclusion:
                labelkeep = 'Excluding volume with a cutoff of ' + str(simulationParams['excludedVolumeCutOff']) + ' μm'
            else:
                labelkeep = 'Without excluded volume'       
                simulationParams['excludedVolumeCutOff'] = 0

            p0 = RCLPolymer(**polymerParams)
            results = {**polymerParams, **simulationParams}
            mc = Experiment(p0, results, simulationParams, 'persistentDSB')

            if first_time:
                fieldnames = ['experimentSetID']+list(mc.results)
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                first_time = False
            writer.writerow({**{'experimentSetID' : str(i)}, **mc.results})
            
            plt.title(labelkeep)
            time = np.arange(simulationParams['numSteps']+1)*simulationParams['dt']
            
            from itertools import combinations
            z = combinations([chr(97 + i//2) + str(1 + i%2) for i in range(2*len(p0.subDomainnumMonomers))], r = 2)
            names = [i for i in z]

            for i in range(15):
                distanceA2B1 = mc.results['MeanInterBreakDistance'][:,i]
                plt.plot(time,distanceA2B1,'-',label='d(%s, %s)' % names[i])
            
        plt.legend()        
        plt.show()    
    

def potentialEnergyplot():
    q0 = RCLPolymer(**polymerParams)
    dt = 0.01
    T = 120
    numSteps = int(T/dt)
    cutoffs = np.array([0.5, 1., 1.5, 2., 2.5, 3.])*0.2
    numRealisations = 100
    
    for cutoff in cutoffs:
        
        p0 = q0.copy()
        kappa = 3*0.008/(p0.b**2)
        repulsionForce = lambda polymer : - kappa * ExcludedVolume(polymer, cutoff)
        p0.addnewForce(repulsionForce)
        
        energy = np.zeros((numRealisations,numSteps))
        for i in range(numRealisations):
            for j in range(numSteps):
                p0.step(1,dt,0.008)
                energy[i][j] = p0.potentialEnergy()
        meanEnergy = np.mean(energy, axis = 0)
        plt.plot(np.linspace(0,T,numSteps),meanEnergy,label='Cut-off: %s $\mu$m' % cutoff)
    
    plt.axvline(x=p0.relaxTime(0.008),color='r')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    
    DSBclustering(polymerParams,simulationParams)