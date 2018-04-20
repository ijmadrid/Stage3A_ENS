# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 10:07:26 2018

@author: ignacio
"""

import csv
from modules.Polymers import RCLPolymer
from modules.Experiment import Experiment

polymerParams = dict(numMonomers = 100, # TODO (.., ..., ...)
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 25,
                     keepCL      = True
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 5, 
                        dt                = 0.01,
                        dt_relax          = 0.01,
#                        numSteps          = 500,
#                        excludedVolumeCutOff = 0.2,
                        waitingSteps = 200,
                        numMaxSteps = 500,
                        encounterDistance = 0.1,
                        genomicDistance = 10,
                        Nb = 2,
                        selectedSubDomain = 0
                        )

with open('results/test_results.csv', 'w') as csvfile:
    
    for i in range(5):
        p0 = RCLPolymer(**polymerParams)
        results = {**polymerParams, **simulationParams}
        mc = Experiment(p0, results, simulationParams,"EncounterSimulation")
        
        if i == 0:
            fieldnames = ['experimentSetID']+list(mc.results)
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
        writer.writerow({**{'experimentSetID' : str(i)}, **mc.results})
