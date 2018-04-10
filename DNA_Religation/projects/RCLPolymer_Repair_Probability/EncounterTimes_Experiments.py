# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:48:01 2018

@author: ignacio
"""

## To import the Polymer and Simulation modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts

from modules.Simulation import EncounterSimulation
from modules.Polymers import RCLPolymer

from time import strftime
import pickle

def plot_bar_from_counter(counter, ax=None):
    """"
    This function creates a bar plot from a counter.

    :param counter: This is a counter object, a dictionary with the item as the key
     and the frequency as the value
    :param ax: an axis of matplotlib
    :return: the axis wit the object in it
    """

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    frequencies = counter.values()
    names = counter.keys()

    x_coordinates = np.arange(len(counter))
    ax.bar(x_coordinates, frequencies, align='center')

    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(names)

    return ax


nb_monomers         = 100

numRealisations     = 500
numConectors        = 4  #number of added connectors

maxIterationsPerExperiment = 5000
dt                  = 0.01
dt_relax            = 0.05

dimension           = 3
diffusionCte        = 0.008
b                   = 0.2

genomicDistance = 10
encounterDistance = 0.1

waitingSteps = 200

def plotFET(mc):
    plt.figure()
    loc, scale = sts.expon.fit(mc.FETs)    
    plt.hist(mc.FETs,bins='auto',normed=True)
    x = np.linspace(mc.FETs.min(),mc.FETs.max(), 100)
    plt.plot(x, sts.expon.pdf(x,loc=loc,scale=scale),'r-', label="Fitted exponential ($\lambda \exp(-\lambda x)$)  :  $\lambda$=%5.3f" % (1/scale))
    plt.title("Distribution of the First Time Encounter")
    plt.xlabel('time (sec)')
    plt.ylabel('distribution')
    plt.legend()
    plt.show()
    
def FET_Simulation(encounterDistance,waitingSteps,numRealisations, keepCL):
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    p0 = RCLPolymer(nb_monomers, dimension, b, numConectors, keepCL)
         
    mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment,2,genomicDistance,encounterDistance,waitingSteps)
    mc.run()
        
    plotFET(mc)
    
    halfCI = 1.96*np.std(mc.FETs)/np.sqrt(numRealisations)
    print('Mean FTE : '+str(np.mean(mc.FETs))+' ± '+str(halfCI))
    
    #
    events = mc.events
    plot_bar_from_counter(events)
    plt.show()

    filepath = 'results/FET_Distribution__'+'keepCL_'+str(keepCL)+\
    str(nb_monomers)+'monomers_'+str(numRealisations)+'iterations'+\
    date+'.pkl'
    
    with open(filepath, 'wb') as output:
        pickle.dump(mc, output, pickle.HIGHEST_PROTOCOL)
    
#    proba = events.get('Repair')/sum(events.values())

def openFETtest(file):
    with open(file, 'rb') as input:
        return pickle.load(input, )


def proba_vs_genomicDistance(nb_monomers,gmax,gStep,test_epsilons,numRealisations,errorbars=False):
#    gmax = nb_monomers - 3
    gmin = 2
    
    probas = []
    demiCI = []
    
    plt.figure()
    plt.xlabel('genomic distance (in number of monomers)')
    plt.ylabel(r'$\mathbb{P}$(Repair)')
    
    xg = np.arange(gmin,gmax,gStep,dtype=int)
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    output = np.empty((len(test_epsilons),3,len(xg)))
    for i, eps in enumerate(test_epsilons):
        
        probas = []
        demiCI = []
        for g in xg:
            p0 = RCLPolymer(nb_monomers, dimension, b, numConectors)
            mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, g, eps, waitingSteps)
            mc.run()
            probas.append(mc.repairProba[0])
            demiCI.append(mc.repairProba[1])
        
        probas = np.array(probas)
        demiCI = np.array(demiCI)
        
        output[i][0] = xg
        output[i][1] = probas
        output[i][2] = demiCI
        
        if errorbars:
            plt.errorbar(x=xg, y=probas, yerr=demiCI,
                         fmt='-o', label=r'$\varepsilon = $ '+str(eps), capsize = 4)
        else:
            plt.plot(xg,probas,'-o',label=r'$\varepsilon = $ '+str(eps))

    np.save('results/proba_vs_genomicDistance__'+
            str(nb_monomers)+'monomers_'+
            str(numRealisations)+'iterations_'+date+
            '.npy',output)
    
    plt.legend()        
    plt.show()



def proba_vs_encounterDistance(numRealisations):
    probas = []
    demiCI = []
    
    testDistances = np.arange(0.004,0.006,0.00025)
    date = strftime("%Y_%m_%d_%H_%M")
    
    plt.figure()
    plt.xlabel('Encounter distance ($\mu$ m)')
    plt.ylabel('$\mathbb{P}$(Repair)')
    
    for epsilon in testDistances:
        p0 = RCLPolymer(nb_monomers, dimension, b, numConectors)
        mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, genomicDistance, epsilon,waitingSteps)
        mc.run()
        probas.append(mc.repairProba[0])
        demiCI.append(mc.repairProba[1])
    
    probas = np.array(probas)
    
    output = np.zeros((3,len(probas)))
    output[0] = testDistances
    output[1] = probas
    output[2] = demiCI
    
    np.save('results/proba_vs_encounterDistance__'+
            str(nb_monomers)+'monomers_'+
            str(numRealisations)+'iterations_'+date+
            '.npy',output)
    
    

    plt.plot(testDistances, probas, '-o')
    plt.show()


def proba_vs_conectorsNumber(numRealisations,nb_monomers,test_genomic_distances,x_Nc,keepCL,errorbars=False):
    probas = []
    demiCI = []

    output = np.zeros((len(test_genomic_distances),3,len(x_Nc)))
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    plt.figure()
    
    for i, genomicDistance in enumerate(test_genomic_distances):
        probas = []
        demiCI = []
        
        for Nc in x_Nc:
            p0 = RCLPolymer(nb_monomers, dimension, b, Nc, keepCL)
            mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, genomicDistance, encounterDistance, waitingSteps)
            mc.run()
            probas.append(mc.repairProba[0])
            demiCI.append(mc.repairProba[1])
        
        probas = np.array(probas)
        demiCI = np.array(demiCI)
        
        output[i][0] = x_Nc
        output[i][1] = probas
        output[i][2] = demiCI

        np.save('results/proba_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(nb_monomers)+'monomers_'+
            str(numRealisations)+'iterations'+
            date+'.npy',output)
        
        if errorbars:
            plt.errorbar(x=x_Nc, y=probas, yerr=demiCI,
                     fmt='-o', label=r'$g = $ '+str(genomicDistance), capsize = 4)
        else:
            plt.plot(x_Nc,probas,'-o',label=r'$g = $ '+str(genomicDistance))

    np.save('results/proba_vs_conectorsNumber__'+
            'keepCL_'+str(keepCL)+
            str(nb_monomers)+'monomers_'+
            str(numRealisations)+'iterations'+
            date+'.npy',output)
    
    plt.legend()
    plt.show()

def proba_vs_DSBNumber(numRealisations=500, nb_monomers=50, genomicDistance=4):
    probas = []
    demiCI = []
    NbMax = 1+int((nb_monomers-1)/(genomicDistance +1))
    for Nb in range(1,NbMax):
        p0 = RCLPolymer(nb_monomers, dimension, b, numConectors)
        mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment, Nb, genomicDistance, encounterDistance, waitingSteps)
        mc.run()
        probas.append(mc.repairProba[0])
        demiCI.append(mc.repairProba[1])
    
    probas = np.array(probas)
    
    plt.figure()
    plt.plot(np.arange(1,NbMax),probas)
    plt.plot(np.arange(1,NbMax),probas-demiCI,'r--',lw=0.4)
    plt.plot(np.arange(1,NbMax),probas+demiCI,'r--',lw=0.4)
    plt.show()