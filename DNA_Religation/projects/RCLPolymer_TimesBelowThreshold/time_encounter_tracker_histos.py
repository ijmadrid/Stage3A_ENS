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
from time import strftime
#import pickle
import csv

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

                     
polymerParams = dict(numMonomers = 100, # np.array([100,100]),
                     dim         = 3,
                     b           = 0.6,
                     Nc          = 50, #NcMatrix,
                     keepCL      = False
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 100, 
                        dt                = 0.01,
                        dt_relax          = 0.01,
                        numSteps          = 6000, #15000,
                        excludedVolumeCutOff = 0,
#                        excludedVolumeSpringConstant = 0.6,
                        waitingSteps = 500,
#                        numMaxSteps = 12000,
                        encounterDistance = 0.05,
#                        genomicDistance = 59 - 30 - 1,
                        Nb = 2,
                        Nc_inDamageFoci = 0,
                        distanceThreshold = 0.01
#                        A1 = 10,
#                        B1 = 80
#                        selectedSubDomain = 0
                        )

#x_Nb = np.arange(2,10)
x_Nc = np.arange(2,40)
x_g = [3,12,21,30]

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
    filename = date + '_trackDSB_fixedDSBs' + '.csv'

    first_time = True
    with open('results/'+filename, 'w') as csvfile:

        for i, Nc in enumerate(x_Nc):
  
            print("\n Simulation for Nc = %s CLs " % Nc)
            polymerParams['Nc'] = Nc
            
#            ### ADAPTIVE ENCOUNTER DISTANCE
#            xi = 2*(Nc)/((N-1)*(N-2))
##                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
#            simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
#                
#            simulationParams['distanceThreshold'] = simulationParams['encounterDistance']
                        
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
    
##    trackDistanceDistribution(polymerParams, simulationParams, x_Nc, x_g)
#
#    x_Nc = np.arange(5,90,5)
#    x_g = [5,40]
#    trackDistanceDistribution(polymerParams, simulationParams, x_Nc, x_g)
    
#    x_Nc = np.arange(20,60,5)
    simulationParams['A1'] = 30
    simulationParams['B1'] = 59
#    trackDefinedDSBtrajectory(polymerParams, simulationParams, x_Nc)
#    
    
    p0 = RCLPolymer(**polymerParams)
    
#    simulationParams['waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt_relax']).astype(int)

    results = {**polymerParams, **simulationParams}

    mc = Experiment(p0, results, simulationParams,"trackDSB")
                    
    above = mc.results['aboveTimes']
    below = mc.results['belowTimes']
    
    plt.figure()
    plt.hist(below['a1-a2'])
    
    
    aboveTimes = {'a1-a2': np.array([1.0050e+00, 4.2000e-01, 9.5000e-02, 5.0000e-03, 1.0000e-02, 4.5000e-02, 1.6500e-01, 5.0000e-03, 5.5000e-02, 5.0000e-03, 1.0000e-02, 4.7500e-01, 8.5000e-02, 1.5000e-02, 8.0000e-02, 5.0000e-03, 5.0000e-03, 5.5000e-02, 4.8500e-01, 7.0000e-02, 9.5000e-02, 5.0000e-03, 1.0000e-02, 1.7000e-01, 1.1500e-01, 5.0000e-03, 8.9000e-01, 4.0000e-02, 2.0000e-02, 5.0000e-03, 2.5000e-01, 5.0000e-03, 5.0000e-03, 2.0000e-02, 5.0000e-03, 5.0000e-03, 2.0000e-02, 1.5000e-02, 1.6555e+01, 1.0000e-02, 1.5000e-01, 5.0000e-03, 5.0000e-03, 5.0000e-03, 2.1760e+01, 2.5000e-02, 5.0000e-03, 2.5000e-02, 4.5000e-02, 1.5000e-02, 5.0000e-03, 2.5000e-02, 5.0000e-03, 3.5500e+00, 4.5000e-02, 7.5000e-02, 4.5000e-02, 5.0000e-03, 4.0000e-02, 5.0000e-03, 2.5000e-02, 3.5000e-02, 1.0300e+00, 5.0000e-03, 5.0000e-02, 5.0000e-03, 1.3230e+01, 1.1000e-01, 4.4700e+00, 5.0000e-03]), 
                  'a1-b1': np.array([], dtype=np.float64), 
                  'a1-b2': np.array([0.03 , 0.005, 0.025, 0.07 , 0.005]), 
                  'a2-b1': np.array([0.005, 0.005, 0.005]), 
                  'a2-b2': np.array([0.01 , 0.05 , 0.005, 0.015]), 
                  'b1-b2': np.array([5.0000e-03, 5.0000e-03, 1.0000e-02, 5.0000e-03, 2.0000e-02, 7.5000e-02, 1.5000e-02, 1.0000e-02, 1.3765e+01, 2.6500e-01, 1.0000e-02, 5.0000e-03, 5.0000e-03, 1.0000e-02, 2.0000e-02, 2.2050e+00, 3.0000e-02, 3.0000e-02, 4.5000e-02, 3.0000e-02, 2.0000e-02, 1.0000e-02, 1.6000e+00, 5.0000e-03, 4.0000e-02, 1.0000e-02, 5.0000e-03, 2.5000e-02, 1.5000e-02, 6.6000e-01, 4.3500e-01, 5.0000e-03, 1.0500e-01, 1.0500e-01, 9.9500e-01, 1.1750e+00, 1.5000e-02, 5.0000e-03, 2.5000e-02, 5.0400e+00, 1.5000e-02, 1.9000e-01, 5.0000e-03, 1.0000e-02, 5.0000e-03, 5.0000e-03, 5.0000e-03, 2.0000e-02, 5.0000e-03, 5.0000e-03, 2.3500e-01, 5.0000e-03, 6.0000e-02, 5.0000e-03, 5.0000e-03, 1.2000e-01, 2.0000e-02, 5.0000e-03, 2.1000e-01, 5.0000e-03, 1.2500e-01])}

    
    
    belowTimes = {'a1-a2': np.array([0.01 , 0.005, 0.015, 0.005, 0.005, 0.025, 0.015, 0.045, 0.005, 0.015, 0.015, 0.005, 0.015, 0.01 , 0.005, 0.015, 0.01 , 0.03 , 0.005, 0.015, 0.005, 0.055, 0.01 , 0.025, 0.015, 0.005, 0.03 , 0.02 , 0.005, 0.015, 0.005, 0.005, 0.005, 0.005, 0.01 , 0.025, 0.005, 0.025, 0.01 , 0.005, 0.02 , 0.005, 0.025, 0.04 , 0.03 , 0.005, 0.035, 0.05 , 0.01 , 0.05 , 0.02 , 0.005, 0.005, 0.025, 0.005, 0.005, 0.005, 0.01 , 0.005, 0.01 , 0.02 , 0.01 , 0.03 , 0.055, 0.02 , 0.005, 0.01 , 0.005, 0.075, 0.06 , 0.04 , 0.005, 0.015, 0.01 , 0.035, 0.005, 0.01 , 0.02 , 0.01 , 0.01 , 0.03 , 0.005, 0.005, 0.015, 0.01 , 0.005, 0.015, 0.005, 0.045, 0.01 , 0.01 , 0.03 , 0.005, 0.005, 0.005, 0.005, 0.005, 0.015, 0.01 , 0.005]), 
                  'a1-b1': np.array([0.075]), 
                  'a1-b2': np.array([0.005, 0.015, 0.005, 0.01 , 0.005, 0.035, 0.04 , 0.015, 0.005]), 
                  'a2-b1': np.array([0.01 , 0.08 , 0.01 , 0.015]), 
                  'a2-b2': np.array([0.05 , 0.01 , 0.015, 0.005, 0.01 , 0.01 , 0.015, 0.01 , 0.005]), 
                  'b1-b2': np.array([0.005, 0.035, 0.005, 0.01 , 0.07 , 0.045, 0.015, 0.005, 0.02 , 0.005, 0.02 , 0.015, 0.01 , 0.005, 0.015, 0.005, 0.01 , 0.015, 0.015, 0.015, 0.01 , 0.04 , 0.005, 0.005, 0.02 , 0.005, 0.08 , 0.05 , 0.01 , 0.01 , 0.01 , 0.025, 0.015, 0.01 , 0.06 , 0.005, 0.04 , 0.005, 0.03 , 0.02 , 0.005, 0.005, 0.01 , 0.055, 0.005, 0.005, 0.01 , 0.01 , 0.01 , 0.01 , 0.015, 0.03 , 0.005, 0.03 , 0.005, 0.02 , 0.005, 0.015, 0.01 , 0.005, 0.015, 0.01 , 0.005, 0.005, 0.005, 0.005, 0.02 , 0.02 , 0.005, 0.01 , 0.005, 0.01 , 0.01 , 0.015, 0.085, 0.015, 0.005, 0.035, 0.005, 0.005, 0.01 , 0.005, 0.01 , 0.005, 0.005])}
