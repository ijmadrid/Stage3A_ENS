# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 17:24:05 2018

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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import strftime
#import pickle
import csv

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

                     
polymerParams = dict(numMonomers = 100, # np.array([100,100]), # TODO (.., ..., ...)
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 25, #NcMatrix,
                     keepCL      = False
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 200, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
                        numSteps          = 24000,
                        excludedVolumeCutOff = 0,
#                        excludedVolumeSpringConstant = 0.6,
                        waitingSteps = 250,
#                        numMaxSteps = 12000,
                        encounterDistance = 0.05,
                        genomicDistance = 59 - 30 - 1,
                        Nb = 2,
                        Nc_inDamageFoci = 0,
                        distanceThreshold = 0.05,
                        A1 = 30,
                        B1 = 59
#                        selectedSubDomain = 0
                        )

#x_Nb = np.arange(2,10)
x_Nc = np.arange(2,20,4)
#x_g = [21] #[3,12,21,30]

############################################################################
# SCRIPT FUNCTIONS  ########################################################
############################################################################

def watchOneSimulation(polymerParams, simulationParams):
    p0 = RCLPolymer(**polymerParams)
    results = {}
    mc = Experiment(p0, results, simulationParams, "watchEncounter")
    print(mc.results['FET'])
    return mc

def adaptiveEpsilon(xi, N, b):
    y = 1 + (N*xi)/(2*(1 -xi))
    return 2*np.sqrt(6 * b**2 / (N * 2*(1-xi) * np.sqrt(y**2-1)))


def trackDefinedDSBtrajectory(polymerParams, simulationParams, x_Nc):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_trackDSB' + '.csv'

    first_time = True
    N = polymerParams['numMonomers']
    with open('results/'+filename, 'w') as csvfile:

        for i, Nc in enumerate(x_Nc):
  
            print("\n Simulation for Nc = %s CLs " % Nc)
            polymerParams['Nc'] = Nc
            
            ### ADAPTIVE ENCOUNTER DISTANCE
            xi = 2*(Nc)/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
            simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
            simulationParams['distanceThreshold'] = simulationParams['encounterDistance']
                        
            p0 = RCLPolymer(**polymerParams)
                    
            results = {**polymerParams, **simulationParams}
        
            mc = Experiment(p0, results, simulationParams,"trackDSB")

            if first_time:
                fieldnames = ['experimentSetID']+list(mc.results)
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                first_time = False
            writer.writerow({**{'experimentSetID' : str(i)}, **mc.results})



def trackDistanceDistribution(polymerParams, simulationParams, x_Nc, x_g):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_trackDSB' + '.csv'

    first_time = True
  
    with open('results/'+filename, 'w') as csvfile:

        for i, Nc in enumerate(x_Nc):
  
            print("Simulation for Nc = %s CLs " % Nc)
            polymerParams['Nc'] = Nc
            
            for j, g in enumerate(x_g):
                
                simulationParams['genomicDistance'] = g
                print("and g = %s . " % g)
                
                p0 = RCLPolymer(**polymerParams)
                
#                simulationParams['waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt_relax']).astype(int)
            
                results = {**polymerParams, **simulationParams}
            
                mc = Experiment(p0, results, simulationParams,"trackDSB")

                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+"_"+str(j)}, **mc.results})
    
############################################################################
# EXECUTION  ###############################################################
############################################################################

if __name__ == "__main__":
    
#    trackDistanceDistribution(polymerParams, simulationParams, x_Nc, x_g)
    
    trackDefinedDSBtrajectory(polymerParams, simulationParams, x_Nc)
    
#    p0 = RCLPolymer(**polymerParams)
#    
#    simulationParams['waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt_relax']).astype(int)
#
#    results = {**polymerParams, **simulationParams}
#
#    mc = Experiment(p0, results, simulationParams,"trackDSB")
#                    
#    above = mc.results['aboveTimes']
#    below = mc.results['belowTimes']
#    
#    plt.figure()
#    plt.hist(above['a1-a2'])