# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:48:01 2018

@author: ignacio
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts

from Simulation import EncounterSimulation
from Polymers import RCLPolymer

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

numRealisations     = 1000
numConectors        = 25  #number of added connectors

maxIterationsPerExperiment = 700
dt                  = 0.01
dt_relax            = 0.05

dimension           = 3
diffusionCte        = 0.008
b                   = 0.2

genomicDistance = 10
encounterDistance = 0.1

waitingSteps = 1000

def FET_Simulation(encounterDistance,waitingSteps):
    
    p0 = RCLPolymer(nb_monomers, dimension, b, numConectors)
         
    mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment,2,genomicDistance,encounterDistance,waitingSteps)
    mc.run()
    
#    msd = mc.get_msd_per_monomer()
#    mmsd = mc.get_avg_msd()
#    msrg = mc.get_msrg()
#    b2e = mc.get_b2e()
#    
    # Exponential fit for the FET
#    timeline = mc.timeline
#    f = lambda t, a, l : a*np.exp(-t*l)
#    popt, pcov = curve_fit(f, timeline, self.get_avg_msd())
    
#    
#    # Fitted parameters
#    fit_params = mc.MSD_fit_mean()
#    
#    print('Estimated MSRG : '+str(msrg))
#    print('Real MSRG      : '+str(nb_monomers*b*b/6))
#    print('Estimated b2   : '+str(b2e))
#    print('Real b2        : '+str(b*b))
    
    plt.rc('text', usetex=True)
    plt.figure()
    loc, scale = sts.expon.fit(mc.FETs)    
    plt.hist(mc.FETs,bins=20,normed=True)
    x = np.linspace(mc.FETs.min(),mc.FETs.max(), 100)
    plt.plot(x, sts.expon.pdf(x,loc=loc,scale=scale),'r-', label=r"Fitted exponential ($\lambda \exp(-\lambda x)$)  :  $\lambda$=%5.3f" % (1/scale))
    plt.title("Distribution of the First Time Encounter")
    plt.xlabel('time (sec)')
    plt.ylabel('distribution')
    plt.legend()
    plt.show()
    
    halfCI = 1.96*np.std(mc.FETs)/np.sqrt(numRealisations)
    print('Mean FTE : '+str(np.mean(mc.FETs))+' ± '+str(halfCI))
    
    #
    events = mc.events
    plot_bar_from_counter(events)
    plt.show()
    
    return mc
    
#    proba = events.get('Repair')/sum(events.values())

def proba_vs_genomicDistance(nb_monomers,gmax,numRealisations):
#    gmax = nb_monomers - 3
    gmin = 2
    
    probas = []
    demiCI = []
    
    plt.figure()
    plt.xlabel('genomic distance (in number of monomers)')
    plt.ylabel(r'$\mathbb{P}$(Repair)',rotation=0)
    
    for eps in [0.2,0.1,0.05,0.001]:
        
        probas = []
        demiCI = []
        for g in np.arange(gmin,gmax,dtype=int):
            p0 = RCLPolymer(nb_monomers, dimension, b, numConectors)
            mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, g, eps, waitingSteps)
            mc.run()
            probas.append(mc.repairProba[0])
            demiCI.append(mc.repairProba[1])
        
        probas = np.array(probas)
        demiCI = np.array(demiCI)
        
        plt.errorbar(x=np.arange(gmin,gmax), y=probas, yerr=demiCI,
                     fmt='-o', label=r'$\varepsilon = $ '+str(eps), capsize = 4)
#        plt.plot(np.arange(gmin,gmax),probas-demiCI,'r--',lw=0.4)
#        plt.plot(np.arange(gmin,gmax),probas+demiCI,'r--',lw=0.4)

    plt.legend()        
    plt.show()
    
    
    return (probas,demiCI)


def proba_vs_encounterDistance(numRealisations):
    probas = []
    demiCI = []
    
    testDistances = np.arange(0.001,b,0.002)
 
    for epsilon in testDistances:
        p0 = RCLPolymer(nb_monomers, dimension, b, numConectors)
        mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, genomicDistance, epsilon,waitingSteps)
        mc.run()
        probas.append(mc.repairProba[0])
        demiCI.append(mc.repairProba[1])
    
    probas = np.array(probas)
    
    plt.figure()
    plt.errorbar(x=testDistances, y=probas, yerr=demiCI,
                 fmt='-o', capsize = 4)
    plt.show()


def proba_vs_conectorsNumber(numRealisations,nb_monomers):
    probas = []
    demiCI = []
    plt.figure()
    
    for genomicDistance in [2,5,10,50]:
        probas = []
        demiCI = []
        
        for Nc in range(5,nb_monomers):
            p0 = RCLPolymer(nb_monomers, dimension, b, Nc)
            mc = EncounterSimulation(dt, diffusionCte, p0, dt_relax, numRealisations, maxIterationsPerExperiment, 2, genomicDistance, encounterDistance, waitingSteps)
            mc.run()
            probas.append(mc.repairProba[0])
            demiCI.append(mc.repairProba[1])
        
        probas = np.array(probas)
        demiCI = np.array(demiCI)
        
        plt.errorbar(x=np.arange(5,nb_monomers), y=probas, yerr=demiCI,
                     fmt='-o', label=r'$g = $ '+str(genomicDistance), capsize = 4)
        #plt.plot(np.arange(2,2*nb_monomers),probas-demiCI,'r--',lw=0.4)
        #plt.plot(np.arange(2,2*nb_monomers),probas+demiCI,'r--',lw=0.4)
        
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