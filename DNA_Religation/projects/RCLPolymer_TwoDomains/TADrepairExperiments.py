# -*- coding: utf-8 -*-
"""
Created on Wed May 23 17:47:33 2018

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
NcMatrix[0,0] = 20
NcMatrix[1,1] = 20
NcMatrix[0,1] = NcMatrix[1,0] = 0

polymerParams = dict(numMonomers = np.array([50,50]),
                     dim         = 3,
                     b           = 0.2,
                     Nc          = NcMatrix,
                     keepCL      = True
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 500, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
#                        numSteps          = 12000,
#                        excludedVolumeCutOff = 0.15,
#                        excludedVolumeSpringConstant = 0.60,
                        waitingSteps = 250,
                        numMaxSteps = 12000,
                        encounterDistance = 0.05,
                        genomicDistance = 10,
                        Nb = 2,
#                        Nc_inDamageFoci = 5,
                        selectedSubDomain = 0
                        )


x_Nlr = np.arange(6)
x_TADsizes = np.array([100])
x_g = np.arange(2,40,4)

############################################################################
############################################################################

def adaptiveEpsilon(xi, N, b):
    y = 1 + (N*xi)/(2*(1 -xi))
    return 2*np.sqrt(6 * b**2 / (N * 2*(1-xi) * np.sqrt(y**2-1)))

def proba_v_NlrandSizes(polymerParams,simulationParams,x_Nlr,x_TADsizes):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba-v_Nlr_TADsize' + '.csv'

    first_time = True
    
    with open('results/'+filename, 'w') as csvfile:

        for i, nlr in enumerate(x_Nlr):
  
            print("Simulation for Nlr = %s " % nlr)
            polymerParams['Nc'][0,1] = polymerParams['Nc'][1,0] = nlr

            for j, tadsize in enumerate(x_TADsizes):    
                print("Simulation for a parasyte of size %s " % tadsize)
                polymerParams['numMonomers'][1] = tadsize
                
                N = polymerParams['numMonomers'].sum()
                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(polymerParams['Nc'].sum())/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s" % (simulationParams['encounterDistance']))
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams, **{'Nlr' : polymerParams['Nc'][0,1],
                                                                   'TADsize' : polymerParams['numMonomers'][1]}}
                
                mc = Experiment(p0, results, simulationParams,"oneTAD_Repair")

                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})


def proba_v_Nlr_g(polymerParams,simulationParams,x_Nlr,x_g):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba-v_Nlr_g' + '.csv'

    first_time = True
    
    with open('results/'+filename, 'w') as csvfile:

        for i, nlr in enumerate(x_Nlr):
  
            print("Simulation for Nlr = %s " % nlr)
            polymerParams['Nc'][0,1] = polymerParams['Nc'][1,0] = nlr

            for j, g in enumerate(x_g):    
                print("Simulation for g = %s " % g)
                simulationParams['genomicDistance'] = g
                
                N = polymerParams['numMonomers'].sum()
                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(polymerParams['Nc'].sum())/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s" % (simulationParams['encounterDistance']))
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams, **{'Nlr' : polymerParams['Nc'][0,1],
                                                                   'TADsize' : polymerParams['numMonomers'][1]}}
                
                mc = Experiment(p0, results, simulationParams,"oneTAD_Repair")

                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})


if __name__ == "__main__":
    
#    proba_v_NlrandSizes(polymerParams,simulationParams,x_Nlr,x_TADsizes)
    proba_v_Nlr_g(polymerParams,simulationParams,x_Nlr,x_g)
#    p0 = RCLPolymer(**polymerParams)
#    p0.step(500,0.05,0.008)
#    p0.colors[:100] = ['y']*100
#    p0.plot()