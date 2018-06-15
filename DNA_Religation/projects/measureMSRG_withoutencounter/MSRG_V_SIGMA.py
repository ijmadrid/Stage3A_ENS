# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 17:56:05 2018

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
import csv

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

                     
polymerParams = dict(numMonomers = 100,
                     dim         = 3,
                     b           = 0.6,
                     Nc          = 50,
                     keepCL      = True
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 500, 
                        dt                = 0.005,
                        dt_relax          = 0.01,

                        excludedVolumeCutOff = 0,
#                        excludedVolumeSpringConstant = 0.6,
                        waitingSteps = 2000,
                        encounterDistance = 0.10,
                        genomicDistance = 10,
                        Nb = 2,
                        Nc_inDamageFoci = 0,
                        times2sample = [1,500,1000,5000,10000]
#                        A1 = 10,
#                        B1 = 80
                        )

x_sigma = np.linspace(0.10,0.30,num=5)


############################################################################
# SCRIPT FUNCTIONS  ########################################################
############################################################################

def proba_v_sigma(polymerParams,simulationParams,x_sigma):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_msrg-v_sigma' + '.csv'

    first_time = True
        
    with open('results/'+filename, 'w') as csvfile:

        for i, keepCL in enumerate([True]):
            
  
            polymerParams['keepCL'] = keepCL

            for j, sigma in enumerate(x_sigma):    
                print("Simulation for σ = %s " % sigma)
                simulationParams['excludedVolumeCutOff'] = sigma
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                if sigma == 0:
                    expName = "measureMSRG"
                else:
                    expName = "measureMSRG"
                mc = Experiment(p0, results, simulationParams,expName)

                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})

if __name__ == "__main__":
        
    proba_v_sigma(polymerParams,simulationParams,x_sigma)