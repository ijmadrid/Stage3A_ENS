# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:30:51 2018

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
                     Nc          = 12, #NcMatrix,
                     keepCL      = True
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 300, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
#                        numSteps          = 500,
                        #excludedVolumeCutOff = 0.1,
                        waitingSteps = 200,
                        numMaxSteps = 9000,
                        encounterDistance = 0.05,
#                        genomicDistance = 25,
                        Nb = 2,
                        Nc_inDamageFoci = 3
#                        selectedSubDomain = 0
                        )


#x_Nc = np.arange(0,5)

gmax = 30
gStep = 3

errorbars = True

############################################################################
############################################################################


def adaptiveEpsilon(xi, N, b):
    y = 1 + (N*xi)/(2*(1 -xi))
    return 2*np.sqrt(6 * b**2 / (N * 2*(1-xi) * np.sqrt(y**2-1)))

def manipulate(function, arguments):
    return 

def proba_vs_keepDFCL(polymerParams,simulationParams,x_Nc,errorbars=False):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + 'proba_VE_vs_Nc_NcinDF_adptEPSILONandSIGMA' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize' : 'xx-large'})
        rcParams.update({'legend.fontsize': 'xx-large'})
        rcParams.update({'xtick.labelsize': 'large'}) 
        plt.xlabel('Number of connectors in DF')
        plt.ylabel(r'$\mathbb{P}$(Repair)') #('Mean first encounter time (sec)') #
        
        first_time = True
        
        N = polymerParams['numMonomers']

        for i, keepCL in enumerate([True, False]):
            
            if keepCL:
                labelkeep = "Keeping CLs in DF"
                
            else:
                labelkeep = 'Removing CLs in DF'       
            polymerParams['keepCL'] = keepCL
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            probas = np.zeros(len(x_Nc))
            demiCI = np.zeros(len(x_Nc))
            for j, Nc in enumerate(x_Nc):    
                print("Simulation for %s CL in DF" % Nc)
                simulationParams['Nc_inDamageFoci'] = Nc

                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(polymerParams['Nc']+Nc)/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE CUTOFF RADIUS
                simulationParams['excludedVolume'] = 2 * adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s and σ = %s" % (simulationParams['encounterDistance'],simulationParams['excludedVolume']))
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,"EncounterSimulation")
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

def proba_vs_genomicDistance(polymerParams,simulationParams,gmax,gStep,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba_vs_genomicDistance_Ncs_' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
        
        gmin = 2
        
        
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        plt.xlabel('genomic distance (in number of monomers)')
        plt.ylabel(r'$\mathbb{P}$(Repair)')
        
        xg = np.arange(gmin,gmax,gStep,dtype=int)
        
        first_time = True
        
        for i, keepCL in enumerate([True, False]):
            
            if keepCL:
                labelkeep = "Keeping CLs in DF"
                
            else:
                labelkeep = 'Removing CLs in DF'       
            polymerParams['keepCL'] = keepCL
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            probas = np.zeros(len(xg))
            demiCI = np.zeros(len(xg))  

            for j, g in enumerate(xg):    
                print("Simulation for g = %s " % g)
                simulationParams['genomicDistance'] = g

                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,'EncounterSimulation')
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
                
                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
            if errorbars:
                plt.errorbar(x=xg, y=probas, yerr=demiCI,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(xg,probas,'-o',label=labelkeep)
        
        plt.legend()        
        plt.show()


if __name__ == "__main__":
        
#    proba_vs_keepDFCL(polymerParams,simulationParams,x_Nc,errorbars)   ###########
    proba_vs_genomicDistance(polymerParams,simulationParams,gmax,gStep,errorbars)