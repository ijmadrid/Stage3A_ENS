# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 10:29:51 2018

@author: ignacio
"""

import numpy as np

    
def ExcludedVolume(polymer,cutoff,method='spring-like'):
    """
    Return the Exluded Volume Potential Gradient (a callable)
    $ \nabla \phi (R) $
    """
    
    if method == 'spring-like':
        """
        The Excluded volume potential is
        sum_{n != m} { (Rm - Rn)^2 * 1_{dist(Rm,Rn) < epsilon} }
        where epssilon is a cutoff that has to be given in the arguments
        """
        distanceMatrix = polymer.distanceMatrix()
        interactionMatrix = np.where(distanceMatrix < cutoff, -1, 0) + np.eye(polymer.numMonomers)
        interactionMatrix -= np.diag(np.sum(interactionMatrix, axis = 0))
        return np.dot(interactionMatrix,polymer.get_r())
    
    if method == 'lennard-jones':
        """
        The Excluded volume potential is
        sum { U_{LJ}(R_m - R_n) }
        where
        U_{LJ}(r) = 4 * [ (s/r)^12 - (s/r)^6 + 1/4 ] if |r| > 2^(1/6)*s
                  = 0                                otherwise
        and s is a cutoff that has to be given in the arguments
        """
        distanceMatrix = polymer.distanceMatrix() + 0.5*cutoff*np.eye(polymer.numMonomers)
        similarityMatrix = (cutoff/distanceMatrix)**6
        return np.dot(8*(similarityMatrix-1)*(similarityMatrix < 0.5), polymer.get_r())
        
    
    if method == 'exponential-like':
        raise NotImplementedError
    
    else:
        print("Available methods : spring-like, lennard-jones, exponential-like")
        raise NotImplementedError 