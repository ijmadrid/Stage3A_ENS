# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 10:29:51 2018

@author: ignacio
"""

import numpy as np
#from scipy.spatial.distance import pdist, squareform
    
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
        

def LocalExcludedVolume(polymer, exclusionLoci, cutoff):
    """
    Add excluded volume forces in the monomers indicated by exclusionLoci
    """
##    indexes = np.concatenate([np.arange(x-window,x+1) for x in exclusionLoci])
    localPositions = polymer.get_r()[exclusionLoci]
    localDistances = np.array([np.linalg.norm(polymer.get_r() - rn, axis = 1) for rn in localPositions])
    interactionMatrix = np.zeros((polymer.numMonomers,polymer.numMonomers))
    localInteractionMatrix = np.where((localDistances < cutoff) & (localDistances > 0), -1, 0)
    for i, monomer_i in enumerate(exclusionLoci):
#        localInteractionMatrix[i,monomer_i] = - np.sum(localInteractionMatrix[i]) - 1
        interactionMatrix[monomer_i] = localInteractionMatrix[i]
        interactionMatrix[:,monomer_i] = localInteractionMatrix[i]
        interactionMatrix -= np.diag(interactionMatrix.sum(axis=0))
    return np.dot(interactionMatrix,polymer.get_r())


def RepairSphere(polymer, exclusionLoci, cutoff):
    """
    Add excluded volume forces in the monomers indicated by exclusionLoci,
    however exclusionLoci monomers do not feel this repulsion between them
    """
##    indexes = np.concatenate([np.arange(x-window,x+1) for x in exclusionLoci])
    localPositions = polymer.get_r()[exclusionLoci]
    localDistances = np.array([np.linalg.norm(polymer.get_r() - rn, axis = 1) for rn in localPositions])
    interactionMatrix = np.zeros((polymer.numMonomers,polymer.numMonomers))
    localInteractionMatrix = np.where(localDistances < cutoff, -1, 0)
    localInteractionMatrix[:,exclusionLoci] = 0
    for i, monomer_i in enumerate(exclusionLoci):
#        localInteractionMatrix[i,monomer_i] = - np.sum(localInteractionMatrix[i]) - 1
        interactionMatrix[monomer_i] = localInteractionMatrix[i]
        interactionMatrix[:,monomer_i] = localInteractionMatrix[i]
        interactionMatrix += np.diag(interactionMatrix.sum(axis=0))
    
    return np.dot(interactionMatrix,polymer.get_r())    


def ExternalStatiticForce(polymer,source,cutoff):
    """
    Add a external harmonic potential created by an "obstacle"
    """
    r = polymer.get_r()
    distance2source = np.array([np.linalg.norm(r - rn, axis = 1) for rn in source])
    distance2source = np.min(distance2source, axis = 1)
    interactionMatrix = np.zeros((polymer.numMonomers,polymer.numMonomers))

    sourceInteractionMatrix = np.where(distance2source < cutoff, -1, 0)
    for i, monomer_i in enumerate(source):
        sourceInteractionMatrix[i,monomer_i] = - np.sum(sourceInteractionMatrix[i]) - 1
        interactionMatrix[monomer_i] = sourceInteractionMatrix[i]
    return np.dot(interactionMatrix,polymer.get_r())