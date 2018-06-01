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
                     keepCL      = True
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 250, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
#                        numSteps          = 12000,
                        excludedVolumeCutOff = 0,
#                        excludedVolumeSpringConstant = 0.6,
                        waitingSteps = 250,
                        numMaxSteps = 12000,
                        encounterDistance = 0.05,
                        genomicDistance = 5,
                        Nb = 5,
                        Nc_inDamageFoci = 0
#                        selectedSubDomain = 0
                        )

x_Nb = np.arange(2,10)
x_g = [3,5,7,9]

############################################################################
# SCRIPT FUNCTIONS  ########################################################
############################################################################

def watchOneSimulation(polymerParams, simulationParams):
    p0 = RCLPolymer(**polymerParams)
    results = {}
    mc = Experiment(p0, results, simulationParams, "watchEncounter")
    print(mc.results['FET'])
    return mc

def trackDistanceDistribution(polymerParams, simulationParams, x_Nb, x_g):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_StatsVsNb' + '.csv'

    first_time = True
  
    with open('results/'+filename, 'w') as csvfile:

        for i, Nb in enumerate(x_Nb):
  
            print("Simulation for Nb = %s breaks " % Nb)
            simulationParams['Nb'] = Nb
            
            for j, g in enumerate(x_g):
                
                simulationParams['genomicDistance'] = g
                print("and g = %s . " % g)
                
                p0 = RCLPolymer(**polymerParams)
                
                simulationParams['waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt_relax']).astype(int)
            
                results = {**polymerParams, **simulationParams}
            
                mc = Experiment(p0, results, simulationParams,"EncounterSimulation")

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
    
    stats_vs_Nb(polymerParams, simulationParams, x_Nb, x_g)
#    mc = watchOneSimulation(polymerParams, simulationParams)
#    ani = mc.plot_trajectoire(show=True)
#    
#    plt.figure()
#    plt.plot(mc.results['realtime'],mc.results['a1MSD'])
#    plt.plot(mc.results['realtime'],mc.results['a2MSD'])
#    plt.plot(mc.results['realtime'],mc.results['b1MSD'])
#    plt.plot(mc.results['realtime'],mc.results['b2MSD'])
#    plt.plot(mc.results['realtime'],mc.results['polymerMSD'])
#    plt.show()