# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 16:33:08 2018

@author: ignacio
"""

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
                     b           = 0.2,
                     Nc          = 20,
                     keepCL      = True
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 1, 
                        dt                = 0.005,
                        dt_relax          = 0.01,

                        excludedVolumeCutOff = 0,
#                        excludedVolumeSpringConstant = 0.6,
                        waitingSteps = 2000,
                        encounterDistance = 0.10,
                        genomicDistance = 10,
                        Nb = 2,
                        Nc_inDamageFoci = 0,
#                        times2sample = [1,500,1000,5000,10000]
#                        A1 = 10,
#                        B1 = 80
                        )

numSteps = 24000
D = 0.008
dt = 0.005


p0 = RCLPolymer(**polymerParams)
d = np.zeros((numSteps,6))

mc = Experiment(p0, {}, simulationParams, "break")

for i in range(numSteps):
    p0.step(1, dt, D)
    d[i] = p0.interBreakDistance()[0]

realtime = np.arange(numSteps) * dt

import matplotlib.pyplot as plt
plt.figure()
plt.plot(realtime,d,lw=1)
plt.hlines(0.1,0,realtime[-1])