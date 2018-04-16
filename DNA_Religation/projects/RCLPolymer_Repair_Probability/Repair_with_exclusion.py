# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 12:01:07 2018

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
import scipy.stats as sts
import pickle

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
                        numRealisations   = 500, 
                        dt                = 0.01,
                        dt_relax          = 0.01,
                        numSteps          = 500,
                        excludedVolumeCutOff = 0.2,
                        waitingSteps = 200,
                        numMaxSteps = 20000,
                        encounterDistance = 0.1,
                        genomicDistance = 10,
                        Nb = 2
                        )

############################################################################
############################################################################


############################################################################
# FUNCTIONS ################################################################
############################################################################

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


def plotFET(mc):
    FETs = mc.results['FETs']
    FETs = FETs.flatten()
    plt.figure()
    loc, scale = sts.expon.fit(FETs)    
    plt.hist(FETs,bins='auto',normed=True)
    x = np.linspace(FETs.min(),FETs.max(), 100)
    plt.plot(x, sts.expon.pdf(x,loc=loc,scale=scale),'r-', label="Fitted exponential ($\lambda \exp(-\lambda x)$)  :  $\lambda$=%5.3f" % (1/scale))
    plt.title("Distribution of the First Time Encounter")
    plt.xlabel('time (sec)')
    plt.ylabel('distribution')
    plt.legend()
    plt.show()


def FET_Simulation(polymerParams,simulationParams):
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    p0 = RCLPolymer(**polymerParams)
    results = {}
    mc = Experiment(p0, results, simulationParams,"Encounter_withExcludedVolume")
        
    plotFET(mc)
    
    halfCI = 1.96*np.std(mc.results['FETs'])/np.sqrt(simulationParams['numRealisations'])
    print('Mean FTE : '+str(np.mean(mc.results['FETs']))+' ± '+str(halfCI))
    
    events = mc.results['eventsCounter']
    plot_bar_from_counter(events)
    plt.show()

    filepath = 'results/FET_Distribution_Experiment_'+date+'.pkl'
    
    with open(filepath, 'wb') as output:
        pickle.dump(mc, output, pickle.HIGHEST_PROTOCOL)
    
    return mc


def saveasPickle(filepath, object2save):
    with open(filepath, 'wb') as output:
        pickle.dump(object2save, output, pickle.HIGHEST_PROTOCOL)
        
def openPickle(file):
    with open(file, 'rb') as input:
        return pickle.load(input)



def proba_vs_genomicDistance(polymerParams,simulationParams,gmax,gStep,test_epsilons,errorbars=False):
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
        
        simulationParams['encounterDistance'] = eps
        
        probas = []
        demiCI = []
        for g in xg:    
            print("Simulation for g =",g,"and epsilon =",eps)
            p0 = RCLPolymer(**polymerParams)
            results = {}
            simulationParams['genomicDistance'] = g
            mc = Experiment(p0, results, simulationParams,"Encounter_withExcludedVolume")
            
            probas.append(mc.results['repair_probability_CI'][0])
            demiCI.append(mc.results['repair_probability_CI'][1])
        
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

    np.save('results/proba_vs_genomicDistance_ExcludedVolume__'+date+'.npy',output)
    
    plt.legend()        
    plt.show()
    

def proba_vs_ExclusioncutoffRadius(polymerParams,simulationParams,test_cutoffs,errorbars=False):
    
    Allresults = []
    
    probas = []
    demiCI = []
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    plt.figure()
    plt.xlabel('Cutoff radius (micron)')
    plt.ylabel(r'$\mathbb{P}$(Correct repair)')
    
    
    for threshold in test_cutoffs:
        
        print('Simulation for a cutoff of',threshold)
        
        p0 = RCLPolymer(**polymerParams)
        simulationParams['excludedVolumeCutOff'] = threshold
        results = {}
        mc = Experiment(p0, results, simulationParams, "Encounter_withExcludedVolume")
    
        Allresults.append(mc.results)
        
        probas.append(mc.results['repair_probability_CI'][0])
        demiCI.append(mc.results['repair_probability_CI'][1])
    
    if errorbars:
        plt.errorbar(x=test_cutoffs, y = probas, yerr=demiCI,
                     fmt='-o', capsize = 4)
    else:
        plt.plot(test_cutoffs, probas,'-o')
    
    filepath = 'results/Proba_vs_cutoffs__'+date+'.pkl'
    
    saveasPickle(filepath, Allresults)
    
        
#    plt.legend()        
    plt.show()
    

############################################################################
############################################################################


############################################################################
# TESTS ####################################################################
############################################################################

if __name__ == '__main__':
    
    ### SIMPLE SIMULATIONS WITH EXCLUDED VOLUME
    
#    mc = FET_Simulation(polymerParams,simulationParams)
    
#    proba_vs_genomicDistance(polymerParams,simulationParams,15,2,[0.2,0.1],errorbars=True)
    
    proba_vs_ExclusioncutoffRadius(polymerParams,simulationParams,np.arange(0,0.2*5,0.01),errorbars=True)
    
#    # Without excluded volume
#    p0 = RCLPolymer(**polymerParams)
#    print(len(p0.forces), 'added forces.')
#    results = {}
#    mc_normal = Experiment(p0, results, simulationParams)
#    print('MRG without volume exclusion :', np.sqrt(mc_normal.results['MSRG']))    
    
############################################################################
############################################################################