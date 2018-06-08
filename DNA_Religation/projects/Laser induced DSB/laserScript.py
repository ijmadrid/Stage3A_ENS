# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 13:14:59 2018

@author: ignacio
"""

######################################
######## LASER EXPERIMENT ############
######################################

import numpy as np

DestructionZone = CylindricLaser()
# TO DO , EQUATION OF CYLINDER (DIRECTION, RADIUS)


class Laser:
    
    def __init__(self, end1, end2, radius):
        """
        Creates a laser centered in end1, with direction end2 - end1 and radius radius
        """
        self.P1 = end1
        self.P2 = end2
        self.r  = radius

    def distance2(self, x0):
        """
        Point-line distance 
        """        
        d = np.linalg.norm(np.cross(x0 - self.P1, x0 - self.P2))/np.linalg.norm(self.P2 - self.P1)
        return d
        
    def contains(self, point):
        return self.distance2(point) < self.r


class Polymer:
    
    def laserInduceDSB(self, destructionRadius, destructionProbability):
        """
        
        """
        laser = Laser(np.zeros(3), self.centerOfMass(), destructionRadius)
        for m, rm in enumerate(self.get_r()):
            if laser.contains(rm):
                if np.random.rand(1) < destructionProbability:
                    try:
                        self.cut(m, m+1)
