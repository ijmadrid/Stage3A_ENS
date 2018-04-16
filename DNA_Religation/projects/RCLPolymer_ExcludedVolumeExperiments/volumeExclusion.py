# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 12:11:44 2018

@author: ignacio
"""


############################################################################
# PACKAGES AND MODULES #####################################################
############################################################################

## To import the Polymer and Experiments modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from modules.Simulation import EncounterSimulation
from modules.Polymers import RCLPolymer
from modules.Experiment import Experiment
from modules.Forces import ExcludedVolume, LocalExcludedVolume

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import strftime, time

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

polymerParams = dict(numMonomers = 100,
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 5
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 1, 
                        dt                = 0.01,
                        dt_relax          = 0.01,
                        numSteps          = 500
                        )

excludedVolumeCutOff = 5*0.2

############################################################################
############################################################################


############################################################################
# TESTS ####################################################################
############################################################################

if __name__ == '__main__':
    
    ### SIMPLE SIMULATIONS WITH EXCLUDED VOLUME
    
    # Polymer parameters
    p0 = RCLPolymer(**polymerParams)
    
    # Add excluded volume
    # First define the potential gradient
    # It should depend on the polymer only !!    
    kappa = 3*0.008/(p0.b**2)
    
    # Loci dependant kappa
    breakPoints = [5,6,10,11,70,71]
#    amplitude = np.ones(p0.numMonomers)/10.
#    amplitude[breakPoints] = 1.
#    kappa = kappa * amplitude
    
#    repulsionForce = lambda polymer : (ExcludedVolume(polymer, excludedVolumeCutOff, method='spring-like').T * -kappa).T
    
    repulsionForce = lambda polymer : - kappa * LocalExcludedVolume(polymer, breakPoints, excludedVolumeCutOff)
    
    p0.addnewForce(repulsionForce)
    # Simulation
    print(len(p0.forces), 'added forces.')
    results = {}
    mc_ve = Experiment(p0, results, simulationParams)
    for monomer in breakPoints: mc_ve.addMonomer2follow(monomer)
    print('MRG with volume exclusion :', np.sqrt(mc_ve.results['MSRG']))
    
#    # Without excluded volume
#    p0 = RCLPolymer(**polymerParams)
#    print(len(p0.forces), 'added forces.')
#    results = {}
#    mc_normal = Experiment(p0, results, simulationParams)
#    print('MRG without volume exclusion :', np.sqrt(mc_normal.results['MSRG']))    
    
############################################################################
############################################################################