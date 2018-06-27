# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 12:18:22 2018

@author: ignacio
"""

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from modules.Polymers import RCLPolymer
from modules.Experiment import Experiment
from modules.Forces import RepairSphere
import numpy as np


polymerParams = dict(numMonomers = 20,
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 10,
                     keepCL      = False
                     )


p0 = RCLPolymer(**polymerParams)