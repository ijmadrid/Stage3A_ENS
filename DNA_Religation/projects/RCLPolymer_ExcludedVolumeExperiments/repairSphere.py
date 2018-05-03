# -*- coding: utf-8 -*-
"""
Created on Thu May  3 16:57:48 2018

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

NcMatrix = np.ones((2,2),dtype=int)*0
NcMatrix[0,0] = 10
NcMatrix[1,1] = 10
                     
polymerParams = dict(numMonomers = 100, # np.array([100,100]), # TODO (.., ..., ...)
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 10, #NcMatrix,
                     keepCL      = False
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 200, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
#                        numSteps          = 500,
                        excludedVolumeCutOff = 0.1,
                        waitingSteps = 200,
                        numMaxSteps = 6000,
                        encounterDistance = 0.05,
                        genomicDistance = 12,
                        Nb = 2
#                        selectedSubDomain = 0
                        )

#test_distances = np.arange(1,30,3,dtype = int)
#
#gmax = 50
#gStep = 4
#test_epsilons = [0.1]
#
#x_Nc = np.arange(25,45,5)
x_Nc = np.arange(2,33,10)
#TADsizes = [20,50,100,200,300]
#connectivityFraction = 0.002
errorbars = True

############################################################################
############################################################################

def proba_vs_Nc_andVE(polymerParams,simulationParams,x_Nc,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + 'proba_VE_vs_Nc' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        rcParams.update({'legend.fontsize': 'large'})
        plt.xlabel('Number of connectors')
        plt.ylabel(r'$\mathbb{P}$(Repair)') #('Mean first encounter time (sec)') #
        
        first_time = True
        for i, VolumeExclusion in enumerate([True, False]):
            
            if VolumeExclusion:
                labelkeep = 'Excluding volume with a cutoff of ' + str(simulationParams['excludedVolumeCutOff']) + ' μm'
                experimentName = "Encounter_withRepairSphere"
            else:
                labelkeep = 'Without excluded volume'
                experimentName = "EncounterSimulation"        
                simulationParams['excludedVolumeCutOff'] = 0
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            probas = np.zeros(len(x_Nc))
            demiCI = np.zeros(len(x_Nc))
            for j, Nc in enumerate(x_Nc):    
                print("Simulation for Nc =",Nc,"and",labelkeep)
                polymerParams['Nc'] = Nc
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,experimentName)
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
#                mfet[j] = mc.results['meanFET']
#                efet[j] = mc.results['halfCI_FET']
                
                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
            if errorbars:
                plt.errorbar(x=x_Nc, y=probas, yerr=demiCI,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(x_Nc,probas,'-o',label=labelkeep)
        
        plt.legend()        
        plt.show()

if __name__ == "__main__":
        
#    file = 'results/2018_04_20_12_55_VE_proba_vs_genomicDistance_.csv'
#    wanted = {'genomicDistance',
#              'excludedVolumeCutOff',
#              'repair_probability',
#              'repair_CIhalflength',
#              'FETs'}
#    res = opencsv(file,wanted)
#
#    a = [row.split(' ')[1:] for row in res['FETs'][0].split('\n')]
#    b = [s for r in a for s in r]

    
    
#    mc = proba_v_interDomainDistance(polymerParams, simulationParams, TAD_size, test_distances, errorbars)
#    mc = proba_v_connectivityFraction(polymerParams, simulationParams, test_xi, test_TADsize, errorbars)
#    mc = FET_Simulation(polymerParams,simulationParams)

#    proba_vs_genomicDistance(polymerParams,simulationParams,gmax,gStep,errorbars)
#    proba_vs_genomicDistance_andVE(polymerParams,simulationParams,gmax,gStep,errorbars)
    proba_vs_Nc_andVE(polymerParams,simulationParams,x_Nc,errorbars)   ###########
#    print(__name__)
#    DSBclustering(polymerParams,simulationParams)
    
#    proba_vs_interNc(polymerParams,simulationParams,x_Nc,errorbars)
#    TADsizes = [20,50,100,200,300]
#    proba_v_neighborSize(polymerParams,simulationParams,TADsizes,connectivityFraction,x_Nc,errorbars)
