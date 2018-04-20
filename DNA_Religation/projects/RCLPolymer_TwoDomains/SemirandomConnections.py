# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 17:44:03 2018

@author: ignacio
"""


############################################################################
# PACKAGES AND MODULES #####################################################
############################################################################

## To import the Polymer and Experiments modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from modules.Polymers import RCLPolymer
from modules.Experiment import Experiment
#from modules.Forces import ExcludedVolume, LocalExcludedVolume

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import strftime
import scipy.stats as sts
import pickle
import csv

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

NcMatrix = np.ones((2,2),dtype=int)*0
NcMatrix[0,0] = 25
NcMatrix[1,1] = 0
                     
polymerParams = dict(numMonomers = 100, # TODO (.., ..., ...)
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 15,
                     keepCL      = False
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 500, 
                        dt                = 0.01,
                        dt_relax          = 0.01,
#                        numSteps          = 500,
                        excludedVolumeCutOff = 0.1,
                        waitingSteps = 200,
                        numMaxSteps = 1000,
                        encounterDistance = 0.1,
                        genomicDistance = 10,
                        Nb = 2
#                        selectedSubDomain = 0
                        )

test_distances = np.arange(1,30,3,dtype = int)

gmax = 20
gStep = 1
test_epsilons = [0.1]

x_Nc = np.arange(3,30,2)

errorbars = True

############################################################################
############################################################################



# Define two domains


def proba_v_interDomainDistance(polymerParams, simulationParams, TAD_size, test_distances, errorbars):
    """
    Probability vs interDomainDistance
    """
    
    plt.figure()
    rcParams.update({'axes.labelsize': 'xx-large'})
    plt.xlabel('genomic distance between TADs')
    plt.ylabel(r'$\mathbb{P}$(Repair)')
    
    xg = test_distances
    
    date = strftime("%Y_%m_%d_%H_%M")
    
    output = np.empty((3,len(xg)))

    probas = []
    demiCI = []
    
    for g in xg:    
        print("Simulation for 2 TADs of size %s, at a distance % s" % (TAD_size,g))
        tad1_right = 50 - np.ceil(g/2).astype(int)
        tad2_left = 50 + np.floor(g/2).astype(int)
        polymerParams['TADs'] = [(tad1_right-TAD_size,tad1_right),(tad2_left,tad2_left+TAD_size)]
        p0 = RCLPolymer(**polymerParams)
        results = {}
        mc = Experiment(p0, results, simulationParams,"TAD_Repair")
        
        probas.append(mc.results['repair_probability_CI'][0])
        demiCI.append(mc.results['repair_probability_CI'][1])
    
    probas = np.array(probas)
    demiCI = np.array(demiCI)
    
    output[0] = xg
    output[1] = probas
    output[2] = demiCI
    
    if errorbars:
        plt.errorbar(x=xg, y=probas, yerr=demiCI,
                     fmt='-o', capsize = 4)
    else:
        plt.plot(xg,probas,'-o')

    np.save('results/proba_vs_interTADdistance__'+date+'.npy',output)
    
    plt.legend()        
    plt.show()
    
    return mc


def proba_v_connectivityFraction(polymerParams, simulationParams, test_xi, test_TADsize, errorbars):
    """
    Probability vs interDomainDistance
    """
    
    plt.figure()
    rcParams.update({'axes.labelsize': 'xx-large'})
    plt.xlabel(r'$\xi_{TAD}$')
    plt.ylabel(r'$\mathbb{P}$(Repair)')
        
    date = strftime("%Y_%m_%d_%H_%M")
    
    output = np.empty((len(test_TADsize),3,6))

       
    N = polymerParams['numMonomers']

    for i, size in enumerate(test_TADsize):
        print('***************************************')
        print('SIMULATION FOR 2 TADs of size',size)
        Nl = (size-2)*(size-1)/2
        
        probas = []
        demiCI = []
    
        xi_min = 1/Nl
        xi_max = (size/1.5)/Nl
        
        x_xi = np.linspace(xi_min,xi_max,6)
        
        for xi in x_xi:    
            print("Simulation for 2 TADs cross-linked at %s percent" % (xi*100))
            Nc = np.ceil(xi*Nl).astype(int)
            print(" Nc = %s  random cross-links" % Nc)
            
            polymerParams['TADs'] = [(N//2-size,N//2),(N//2+1,N//2+1+size)]
            polymerParams['TADs_Nc'] = [Nc, Nc]
            
            p0 = RCLPolymer(**polymerParams)
            results = {}
            mc = Experiment(p0, results, simulationParams,"TAD_Repair")
            
            print(" Repair proba : ", mc.results['repair_probability_CI'][0])
            probas.append(mc.results['repair_probability_CI'][0])
            demiCI.append(mc.results['repair_probability_CI'][1])
    
        probas = np.array(probas)
        demiCI = np.array(demiCI)
        
        output[i][0] = x_xi
        output[i][1] = probas
        output[i][2] = demiCI
    
        if errorbars:
            plt.errorbar(x=x_xi, y=probas, yerr=demiCI,
                     fmt='-o', label = 'TAD size = '+str(size),capsize = 4)
        else:
            plt.plot(x_xi, probas,'-o',label = 'TAD size = '+str(size))


    np.save('results/proba_vs_TADsizeANDconnectivity__'+date+'.npy',output)
    
    plt.legend()        
    plt.show()
    
    return mc    


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
    rcParams.update({'axes.labelsize': 'xx-large'})
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
    mc = Experiment(p0, results, simulationParams,"EncounterSimulation")
     
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



def proba_vs_genomicDistance(experimentName,polymerParams,simulationParams,gmax,gStep,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba_vs_genomicDistance_' + experimentName + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
        
        gmin = 2
          
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        plt.xlabel('genomic distance (in number of monomers)')
        plt.ylabel(r'$\mathbb{P}$(Repair)')
        
        xg = np.arange(gmin,gmax,gStep,dtype=int)
        
        for i, keepCL in enumerate([True,False]):
            
            if keepCL:
                labelkeep = 'Keeping CL in damage zone'
            else:
                labelkeep = 'Removing CL in damage zone'
            
            polymerParams['keepCL'] = keepCL
            
            probas = np.zeros(len(xg))
            demiCI = np.zeros(len(xg))
            
            for j, g in enumerate(xg):    
                print("Simulation for g =",g,"and ",labelkeep)
                p0 = RCLPolymer(**polymerParams)
                simulationParams['genomicDistance'] = g
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,experimentName)
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
                
                if i == 0 and g == gmin:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                writer.writerow({**{'experimentSetID' : str(i)}, **mc.results})
            
            
            probas = np.array(probas)
            demiCI = np.array(demiCI)
            
            if errorbars:
                plt.errorbar(x=xg, y=probas, yerr=demiCI,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(xg,probas,'-o',label=labelkeep)
        
        plt.legend()        
        plt.show()


def proba_vs_genomicDistance_andVE(polymerParams,simulationParams,gmax,gStep,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_VE_proba_vs_genomicDistance_' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
        
        gmin = 2
        
        
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        plt.xlabel('genomic distance (in number of monomers)')
        plt.ylabel(r'$\mathbb{P}$(Repair)')
        
        xg = np.arange(gmin,gmax,gStep,dtype=int)
        
        for i, VolumeExclusion in enumerate([True,False]):
            
            if VolumeExclusion:
                labelkeep = 'Excluding volume with a cutoff of ' + str(simulationParams['excludedVolumeCutOff']) + ' μm'
                experimentName = "Encounter_withExcludedVolume"
            else:
                labelkeep = 'Without excluded volume'
                experimentName = "EncounterSimulation"        
                simulationParams['excludedVolumeCutOff'] = 0
            probas = np.zeros(len(xg))
            demiCI = np.zeros(len(xg))
            for j, g in enumerate(xg):    
                print("Simulation for g =",g,"and ",labelkeep)
                p0 = RCLPolymer(**polymerParams)
                simulationParams['genomicDistance'] = g
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,experimentName)
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
                
                if i == 0 and g == gmin:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
            if errorbars:
                plt.errorbar(x=xg, y=probas, yerr=demiCI,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(xg,probas,'-o',label=labelkeep)
        
        plt.legend()        
        plt.show()


def proba_vs_Nc_andVE(polymerParams,simulationParams,x_Nc,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_VE_proba_vs_Nc' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        rcParams.update({'legend.fontsize': 'large'})
        plt.xlabel('genomic distance (in number of monomers)')
        plt.ylabel(r'$\mathbb{P}$(Repair)')
        
        first_time = True
        for i, VolumeExclusion in enumerate([True,False]):
            
            if VolumeExclusion:
                labelkeep = 'Excluding volume with a cutoff of ' + str(simulationParams['excludedVolumeCutOff']) + ' μm'
                experimentName = "Encounter_withExcludedVolume"
            else:
                labelkeep = 'Without excluded volume'
                experimentName = "EncounterSimulation"        
                simulationParams['excludedVolumeCutOff'] = 0
            probas = np.zeros(len(x_Nc))
            demiCI = np.zeros(len(x_Nc))
            for j, Nc in enumerate(x_Nc):    
                print("Simulation for Nc =",Nc,"and",labelkeep)
                polymerParams['Nc'] = Nc
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,experimentName)
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
                
                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
            if errorbars:
                plt.errorbar(x=x_Nc, y=probas, yerr=demiCI,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(x_Nc,probas,'-o',label=labelkeep)
        
        plt.legend()        
        plt.show()

def opencsv(file,wantedResults):
    wantedResults = {res : [] for res in wantedResults}
    with open(file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for res in wantedResults:
                try:
                    wantedResults[res].append(float(row[res]))
                except:
                    wantedResults[res].append(row[res])
    return wantedResults
            

if __name__ == "__main__":
    
    file = 'results/2018_04_20_12_55_VE_proba_vs_genomicDistance_.csv'
    wanted = {'genomicDistance',
              'excludedVolumeCutOff',
              'repair_probability',
              'repair_CIhalflength',
              'FETs'}
    res = opencsv(file,wanted)

    a = [row.split(' ')[1:] for row in res['FETs'][0].split('\n')]
    b = [s for r in a for s in r]

    
    
#    mc = proba_v_interDomainDistance(polymerParams, simulationParams, TAD_size, test_distances, errorbars)
#    mc = proba_v_connectivityFraction(polymerParams, simulationParams, test_xi, test_TADsize, errorbars)
#    mc = FET_Simulation(polymerParams,simulationParams)

#    proba_vs_genomicDistance(polymerParams,simulationParams,gmax,gStep,errorbars)
#    proba_vs_genomicDistance_andVE(polymerParams,simulationParams,gmax,gStep,errorbars)
#    proba_vs_Nc_andVE(polymerParams,simulationParams,x_Nc,errorbars)
#    print(__name__)