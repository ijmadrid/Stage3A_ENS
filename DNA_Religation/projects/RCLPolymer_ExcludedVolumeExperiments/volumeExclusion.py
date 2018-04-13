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
from modules.Forces import ExcludedVolume

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import strftime, time

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

polymerParams = dict(numMonomers = 500,
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 0
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.08,
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
    kappa = 3*0.08/(p0.b**2)
    repulsionForce = lambda polymer : -kappa * ExcludedVolume(polymer, excludedVolumeCutOff, method='spring-like')
    p0.addnewForce(repulsionForce)
    # Simulation
    print(len(p0.forces), 'added forces.')
    results = {}
    mc_ve = Experiment(p0, results, simulationParams)
    print('MRG with volume exclusion :', np.sqrt(mc_ve.results['MSRG']))
    
#    # Without excluded volume
#    p0 = RCLPolymer(**polymerParams)
#    print(len(p0.forces), 'added forces.')
#    results = {}
#    mc_normal = Experiment(p0, results, simulationParams)
#    print('MRG without volume exclusion :', np.sqrt(mc_normal.results['MSRG']))    
    
############################################################################
############################################################################