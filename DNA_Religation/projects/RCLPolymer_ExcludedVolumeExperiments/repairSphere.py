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
                     keepCL      = True
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 5, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
#                        numSteps          = 500,
                        excludedVolumeCutOff = 0.1,
                        waitingSteps = 200,
                        numMaxSteps = 9000,
                        encounterDistance = 0.05,
                        genomicDistance = 25,
                        Nb = 2,
                        Nc_inDamageFoci = 1
#                        selectedSubDomain = 0
                        )

#test_distances = np.arange(1,30,3,dtype = int)
#
gmax = 10
gStep = 5
#test_epsilons = [0.1]
#
#x_Nc = np.arange(25,45,5)
#x_Nc = np.arange(5,40,5)
#TADsizes = [20,50,100,200,300]
#connectivityFraction = 0.002
errorbars = True
x_sigma = np.array([0, 0.1, 0.2])
############################################################################
############################################################################


def adaptiveEpsilon(xi, N, b):
    y = 1 + (N*xi)/(2*(1 -xi))
    return 2*np.sqrt(6 * b**2 / (N * 2*(1-xi) * np.sqrt(y**2-1)))

def proba_vs_Nc_andVE(polymerParams,simulationParams,x_Nc,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + 'proba_VE_vs_Nc_NcinDF_adptEPSILONandSIGMA' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize' : 'xx-large'})
        rcParams.update({'legend.fontsize': 'xx-large'})
        rcParams.update({'xtick.labelsize': 'large'}) 
        plt.xlabel('Number of connectors')
        plt.ylabel(r'$\mathbb{P}$(Repair)') #('Mean first encounter time (sec)') #
        
        first_time = True
        
        N = polymerParams['numMonomers']
#        Nc0 = x_Nc[0]
#        xi0 = 2*Nc0/((N-1)*(N-2))
#        eps0 = simulationParams['encounterDistance']
#        sigma0 = simulationParams['excludedVolumeCutOff']
                
        for i, VolumeExclusion in enumerate([True, False]):
            
            if VolumeExclusion:
                labelkeep = r"Excluding volume ($\xi$-adaptive $\epsilon$ and $\sigma$)"
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
                print("Simulation for Nc =",Nc)
                polymerParams['Nc'] = Nc

                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*Nc/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE CUTOFF RADIUS
                simulationParams['excludedVolume'] = 2 * adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s and σ = %s" % (simulationParams['encounterDistance'],simulationParams['excludedVolume']))
                
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

def proba_vs_genomicDistance_andVE(polymerParams,simulationParams,x_sigma,gmax,gStep,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_VE_proba_vs_genomicDistance_' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
        
        gmin = 2
        
        
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        plt.xlabel('genomic distance (in number of monomers)')
        plt.ylabel(r'$\mathbb{P}$(Repair)')
        
        xg = np.arange(gmin,gmax,gStep,dtype=int)
        
        for i, sigma in enumerate(x_sigma):
            
            if sigma == 0:
                labelkeep = 'Without excluded volume'
                experimentName = "EncounterSimulation"        
                simulationParams['excludedVolumeCutOff'] = 0
            else:
                labelkeep = r"Excluding volume ($\sigma = %s $)" % sigma
                experimentName = "Encounter_withRepairSphere"
                simulationParams['excludedVolumeCutOff'] = sigma
            probas = np.zeros(len(xg))
            demiCI = np.zeros(len(xg))
            for j, g in enumerate(xg):    
                print("Simulation for g =",g,"and ",labelkeep)
                p0 = RCLPolymer(**polymerParams)
                simulationParams['genomicDistance'] = g
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,experimentName)
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
                
                if i == 0 and g == gmin:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
            if errorbars:
                plt.errorbar(x=xg, y=probas, yerr=demiCI,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(xg,probas,'-o',label=labelkeep)
        
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
    proba_vs_genomicDistance_andVE(polymerParams,simulationParams,x_sigma,gmax,gStep,errorbars)
#    proba_vs_Nc_andVE(polymerParams,simulationParams,x_Nc,errorbars)   ###########
#    print(__name__)
#    DSBclustering(polymerParams,simulationParams)
    
#    proba_vs_interNc(polymerParams,simulationParams,x_Nc,errorbars)
#    TADsizes = [20,50,100,200,300]
#    proba_v_neighborSize(polymerParams,simulationParams,TADsizes,connectivityFraction,x_Nc,errorbars)
