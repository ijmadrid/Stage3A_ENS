# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:30:51 2018

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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import strftime
#import pickle
import csv

############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

                     
polymerParams = dict(numMonomers = 100, # np.array([100,100]), # TODO (.., ..., ...)
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 20, #NcMatrix,
                     keepCL      = False
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 100, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
                        numSteps          = 12000,
                        excludedVolumeCutOff = 0.005,
#                        excludedVolumeSpringConstant = 0.30,
                        waitingSteps = 0,
#                        numMaxSteps = 500,
                        encounterDistance = 0.045,
                        genomicDistance = 20,
                        Nb = 2,
                        Nc_inDamageFoci = 2
#                        selectedSubDomain = 0
                        )


#x_Nc = np.arange(3,11)
#x_sigma = np.linspace(0.03,0.25,num=5)
#x_Nc = np.array([3,5,7,9,11,13,15]) #np.arange(3,20,3)
#x_Nd = np.array([0,1,2,3])
#gmax = 12
#gStep = 1

#x_kappa = (3*0.008/(0.2**2))*np.linspace(0, 4, num = 30)

x_g = [2,4,20]
x_Nc = np.arange(20,41,5)
errorbars = True

############################################################################
############################################################################


def adaptiveEpsilon(xi, N, b):
    y = 1 + (N*xi)/(2*(1 -xi))
    return 2*np.sqrt(6 * b**2 / (N * 2*(1-xi) * np.sqrt(y**2-1)))

def manipulate(function, arguments):
    return 

def proba_vs_keepDFCL(polymerParams,simulationParams,x_Nc,errorbars=False):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + 'proba_VE_vs_Nc_NcinDF_adptEPSILONandSIGMA' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize' : 'xx-large'})
        rcParams.update({'legend.fontsize': 'xx-large'})
        rcParams.update({'xtick.labelsize': 'large'}) 
        plt.xlabel('Number of connectors in DF')
        plt.ylabel(r'$\mathbb{P}$(Repair)') #('Mean first encounter time (sec)') #
        
        first_time = True
        
        N = polymerParams['numMonomers']

        for i, keepCL in enumerate([True, False]):
            
            if keepCL:
                labelkeep = "Keeping CLs in DF"
                
            else:
                labelkeep = 'Removing CLs in DF'       
            polymerParams['keepCL'] = keepCL
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            probas = np.zeros(len(x_Nc))
            demiCI = np.zeros(len(x_Nc))
            for j, Nc in enumerate(x_Nc):    
                print("Simulation for %s CL in DF" % Nc)
                simulationParams['Nc_inDamageFoci'] = Nc

                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(polymerParams['Nc']+Nc)/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
#                ### ADAPTIVE CUTOFF RADIUS
#                simulationParams['excludedVolume'] = 2 * adaptiveEpsilon(xi, N, polymerParams['b'])
#                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s and σ = %s" % (simulationParams['encounterDistance'],simulationParams['excludedVolumeCutOff']))
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,"EncounterSimulation")
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
#                mfet[j] = mc.results['meanFET']
#                efet[j] = mc.results['halfCI_FET']
                
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
        
        
    
def proba_vs_NcNd_andKeepCL(polymerParams,simulationParams,x_Nc,x_Nd,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + 'proba_VE_vs_Nc_NcinDF_adptEPSILONandSIGMA' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize' : 'xx-large'})
        rcParams.update({'legend.fontsize': 'xx-large'})
        rcParams.update({'xtick.labelsize': 'large'}) 
        plt.xlabel(r'Number of total connectors $(N_c + N_d)$')
        plt.ylabel(r'$\mathbb{P}$(Repair)') #('Mean first encounter time (sec)') #
        
        first_time = True
        
        N = polymerParams['numMonomers']
#        Nc0 = x_Nc[0]
#        xi0 = 2*Nc0/((N-1)*(N-2))
#        eps0 = simulationParams['encounterDistance']
#        sigma0 = simulationParams['excludedVolumeCutOff']
          
        polymerParams['keepCL'] = False
        
        for i, Nd in enumerate(x_Nd):
            
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            probas = np.zeros(len(x_Nc))
            demiCI = np.zeros(len(x_Nc))
            for j, Nc in enumerate(x_Nc):    
                print("Simulation for total Nc = %s where Nd = %s" % (Nc,Nd))
                polymerParams['Nc'] = Nc - Nd
                simulationParams['Nc_inDamageFoci'] = Nd
                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(Nc+Nd)/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
#                ### ADAPTIVE CUTOFF RADIUS
#                simulationParams['excludedVolume'] = 2 * adaptiveEpsilon(xi, N, polymerParams['b'])
#                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s " % (simulationParams['encounterDistance']))
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,'EncounterSimulation')
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
#                mfet[j] = mc.results['meanFET']
#                efet[j] = mc.results['halfCI_FET']
                
                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
            if errorbars:
                plt.errorbar(x=x_Nc, y=probas, yerr=demiCI,
                             fmt='-o', label=r"$N_d = %s $" % Nd, capsize = 4)
            else:
                plt.plot(x_Nc,probas,'-o',label=r"$N_d = %s $" % Nd)
        
        plt.legend()        
        plt.show()

def proba_vs_Nc_andKeepCL(polymerParams,simulationParams,x_Nc,errorbars=False):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + 'proba_VE_vs_Nc_NcinDF_adptEPSILONandSIGMA' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize' : 'xx-large'})
        rcParams.update({'legend.fontsize': 'xx-large'})
        rcParams.update({'xtick.labelsize': 'xx-large'}) 
        rcParams.update({'ytick.labelsize': 'xx-large'}) 
        plt.xlabel(r'$N_c$')
        plt.ylabel(r'$\mathbb{P}$(Repair)') #('Mean first encounter time (sec)') #
        
        first_time = True
        
        N = polymerParams['numMonomers']

        for i, keepCL in enumerate([True, False]):
            
            if keepCL:
                labelkeep = "Keeping CLs in DF"
                
            else:
                labelkeep = 'Removing CLs in DF'       
            polymerParams['keepCL'] = keepCL
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            probas = np.zeros(len(x_Nc))
            demiCI = np.zeros(len(x_Nc))
            for j, Nc in enumerate(x_Nc):    
                print("Simulation for %s CL in DF" % Nc)
                polymerParams['Nc'] = Nc

                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(Nc)/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
#                ### ADAPTIVE CUTOFF RADIUS
#                simulationParams['excludedVolume'] = 2 * adaptiveEpsilon(xi, N, polymerParams['b'])
#                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s" % (simulationParams['encounterDistance']))
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,"EncounterSimulation")
    #            oneTAD_Repair
                probas[j] = mc.results['repair_probability']
                demiCI[j] = mc.results['repair_CIhalflength']
#                mfet[j] = mc.results['meanFET']
#                efet[j] = mc.results['halfCI_FET']
                
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
        
        
        
def proba_vs_genomicDistance(polymerParams,simulationParams,gmax,gStep,errorbars=False):

    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba_vs_genomicDistance_Ncs_' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
        
        gmin = 2
        
        
        plt.figure()
        rcParams.update({'axes.labelsize': 'xx-large'})
        plt.xlabel('genomic distance (in number of monomers)')
        plt.ylabel(r'$\mathbb{P}$(Repair)')
        
        xg = np.arange(gmin,gmax,gStep,dtype=int)
        
        first_time = True
        
        for i, keepCL in enumerate([True, False]):
            
            if keepCL:
                labelkeep = "Keeping CLs in DF"
                
            else:
                labelkeep = 'Removing CLs in DF'       
            polymerParams['keepCL'] = keepCL
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            probas = np.zeros(len(xg))
            demiCI = np.zeros(len(xg))  

            for j, g in enumerate(xg):    
                print("Simulation for g = %s " % g)
                simulationParams['genomicDistance'] = g

                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,'EncounterSimulation')
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
                plt.errorbar(x=xg, y=probas, yerr=demiCI,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(xg,probas,'-o',label=labelkeep)
        
        plt.legend()        
        plt.show()


############################################################################
############################################################################
###################     FET STUDY ##########################################
############################################################################
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

#    frequencies = counter.values()
#    names = counter.keys()
#
#    x_coordinates = np.arange(len(counter))
#    ax.bar(x_coordinates, frequencies, align='center')
    
    ax.bar(0, counter['Repair'], align='center', color = 'green')
    
    cols = ('yellow', 'orange', 'tomato', 'red')
    b = 0
    for i, mismatch in enumerate(('a1-b1', 'a1-b2', 'a2-b1', 'a2-b2')):
        ax.bar(1,counter['Misrepair_'+mismatch], color = cols[i], bottom = b, label=mismatch)
        b+=counter['Misrepair_'+mismatch]
#    ax.bar(1, counter['Misrepair_a1-b1'], color='yellow')
#    ax.bar(1, counter['Misrepair_a1-b2'], color='orange', bottom=counter['Misrepair_a1-b1'])
#    ax.bar(1, counter['Misrepair_a2-b1'], color='red', bottom=counter['Misrepair_a1-b2'])
#    ax.bar(1, counter['Misrepair_a2-b2'], color='tomato', bottom=counter['Misrepair_a2-b1'])
#    pna = ax.bar(0, counter['Repair'], align='center', color = 'g')

    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(('Repair','Misrepair'))
    plt.legend()#(pm1[0], pm2[0], pm3[0], pm4[0]), ('A1-B1', 'A1-B2', 'A2-B1', 'A2-B2'))
    
    return ax


def plotFET(mc):
    rcParams.update({'axes.labelsize': 'xx-large'})
    rcParams.update({'legend.fontsize': 'xx-large'})
    rcParams.update({'xtick.labelsize': 'xx-large'}) 
    rcParams.update({'ytick.labelsize': 'xx-large'})
    FETs = mc.results['FETs']
    FETs = FETs.flatten()
    FETs = FETs[~np.isnan(FETs)]
    plt.figure()
    fitAmplitude = mc.results['expFit_Amplitude']
    fitLambda    = mc.results['expFit_Rate']
    def exponential(x, amplitude, lambd):
        return amplitude * np.exp(-lambd*x)
    plt.hist(FETs,bins='auto',normed=True)
    x = np.linspace(FETs.min(),FETs.max(), 100)
    plt.plot(x, exponential(x,fitAmplitude,fitLambda),'r-', lw = 4, label="Fit ($ %5.3f \exp(- %5.3f x)$)" % (fitAmplitude, fitLambda))
    plt.title("FET Distribution")
    plt.xlabel('FET (sec)')
    plt.ylabel('Distribution')
    plt.legend()
    plt.show()


def FET_Simulation(polymerParams,simulationParams):
    
#    date = strftime("%Y_%m_%d_%H_%M")
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '__FETSimulationResults' + '.csv'
    
    p0 = RCLPolymer(**polymerParams)
    results = {**polymerParams, **simulationParams}
    mc = Experiment(p0, results, simulationParams,"EncounterSimulation")

    ### SAVE RESULTS
    with open('results/'+filename, 'w') as csvfile:
        fieldnames = list(mc.results)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({**mc.results})
    
    ### PLOT RESULTS
    plotFET(mc)
    
    halfCI = 1.96*np.nanstd(mc.results['FETs'])/np.sqrt(simulationParams['numRealisations'])
    print('Mean FTE : '+str(np.nanmean(mc.results['FETs']))+' ± '+str(halfCI))
    
    events = mc.results['eventsCounter']
    plot_bar_from_counter(events)
    plt.show()

#    filepath = 'results/FET_Distribution_Experiment_'+date+'.pkl'
    
#    with open(filepath, 'wb') as output:
#        pickle.dump(mc, output, pickle.HIGHEST_PROTOCOL)
    
    return mc


def mFET_vs_NcinDF(polymerParams,simulationParams,x_Nc,errorbars=False):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + 'mFET_vs_NcinDF__adptEPSILONandSIGMA' + '.csv'
    
    with open('results/'+filename, 'w') as csvfile:
                
        plt.figure()
        rcParams.update({'axes.labelsize' : 'xx-large'})
        rcParams.update({'legend.fontsize': 'xx-large'})
        rcParams.update({'xtick.labelsize': 'xx-large'}) 
        rcParams.update({'ytick.labelsize': 'xx-large'}) 
        plt.xlabel('Number of connectors in DF')
        plt.ylabel('Mean FET (sec)') #('Mean first encounter time (sec)') #
        
        first_time = True
        
        N = polymerParams['numMonomers']

        for i, keepCL in enumerate([True, False]):
            
            if keepCL:
                labelkeep = "Keeping CLs in DF"
                
            else:
                labelkeep = 'Removing CLs in DF'       
            polymerParams['keepCL'] = keepCL
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
            mfet = np.zeros(len(x_Nc))
            mfet_dx = np.zeros(len(x_Nc))
            for j, Nc in enumerate(x_Nc):    
                print("Simulation for %s CL in DF" % Nc)
                simulationParams['Nc_inDamageFoci'] = Nc

                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(polymerParams['Nc']+Nc)/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE CUTOFF RADIUS
                simulationParams['excludedVolumeCutOff'] = 2 * adaptiveEpsilon(xi, N, polymerParams['b'])
                
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                print("ε = %s and σ = %s" % (simulationParams['encounterDistance'],simulationParams['excludedVolumeCutOff']))
                
                p0 = RCLPolymer(**polymerParams)
                
                
                
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,"EncounterSimulation")
    #            oneTAD_Repair
                mfet[j] = mc.results['meanFET']
                mfet_dx[j] = mc.results['halfCI_FET']
#                mfet[j] = mc.results['meanFET']
#                efet[j] = mc.results['halfCI_FET']
                
                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
            if errorbars:
                plt.errorbar(x=x_Nc, y=mfet, yerr=mfet_dx,
                             fmt='-o', label=labelkeep, capsize = 4)
            else:
                plt.plot(x_Nc,mfet,'-o',label=labelkeep)
        
        plt.legend()        
        plt.show()

################################################################################
################################################################################
########################### PROBA VS SIGMA #####################################
################################################################################
################################################################################
################################################################################  

def proba_v_sigma(polymerParams,simulationParams,x_sigma,errorbars=False):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba-v_sigma' + '.csv'

#    plt.figure()
#    rcParams.update({'axes.labelsize' : 'xx-large'})
#    rcParams.update({'legend.fontsize': 'xx-large'})
#    rcParams.update({'xtick.labelsize': 'xx-large'}) 
#    rcParams.update({'ytick.labelsize': 'xx-large'}) 
#    plt.xlabel('Number of connectors in DF')
#    plt.ylabel('Mean FET (sec)') #('Mean first encounter time (sec)') #

    first_time = True
        
    with open('results/'+filename, 'w') as csvfile:

        for i, keepCL in enumerate([True]):
            
#            if keepCL:
#                labelkeep = "Keeping CLs in DF"
#                
#            else:
#                labelkeep = 'Removing CLs in DF'       
            polymerParams['keepCL'] = keepCL
#            mfet = np.zeros(len(x_Nc))
#            efet = np.zeros(len(x_Nc))
#            p = np.zeros(len(x_sigma))
#            dp = np.zeros(len(x_sigma))
            for j, sigma in enumerate(x_sigma):    
                print("Simulation for σ = %s " % sigma)
                simulationParams['excludedVolumeCutOff'] = sigma
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                if sigma == 0:
                    expName = "EncounterSimulation"
                else:
                    expName = "Encounter_withRepairSphere"
                mc = Experiment(p0, results, simulationParams,expName)

#                p[j] = mc.results['repair_probability']
#                dp[j] = mc.results['repair_CIhalflength']
          
                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
            
#            if errorbars:
#                plt.errorbar(x=x_sigma, y=p, yerr=dp,
#                             fmt='-o', label=labelkeep, capsize = 4)
#            else:
#                plt.plot(x_sigma,p,'-o',label=labelkeep)
##        
#        plt.legend()        
#        plt.show()    
#        

def proba_v_VEkappa(polymerParams,simulationParams,x_sigma,x_kappa,errorbars=False):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba-v_kappa' + '.csv'

    first_time = True
        
    with open('results/'+filename, 'w') as csvfile:

        for i, sigma in enumerate(x_sigma):
  
            print("Simulation for σ = %s " % sigma)
            simulationParams['excludedVolumeCutOff'] = sigma

            for j, kappa in enumerate(x_kappa):    
                print("Simulation for κ = %s " % kappa)
                simulationParams['excludedVolumeSpringConstant'] = kappa
                
                p0 = RCLPolymer(**polymerParams)
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,"Encounter_withRepairSphere")

                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})
  


def proba_v_gNc(polymerParams,simulationParams,x_g,x_Nc,errorbars=False):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_proba-v_gNc_withVE' + '.csv'

    first_time = True
    N = polymerParams['numMonomers']
    
    with open('results/'+filename, 'w') as csvfile:

        for i, g in enumerate(x_g):
  
            print("Simulation for g = %s " % g)
            simulationParams['genomicDistance'] = g

            for j, nc in enumerate(x_Nc):    
                print("Simulation for Nc = %s " % nc)
                polymerParams['Nc'] = nc
                
                ### ADAPTIVE ENCOUNTER DISTANCE
                xi = 2*(nc + simulationParams['Nc_inDamageFoci'] )/((N-1)*(N-2))
#                scaleFactor = np.sqrt( (1-xi0)*np.sqrt(xi0) / ((1-xi)*np.s1qrt(xi)) )
                simulationParams['encounterDistance'] = adaptiveEpsilon(xi, N, polymerParams['b'])
                
#                simulationParams['excludedVolumeCutOff'] = 3*simulationParams['encounterDistance']
                ### ADAPTIVE dt
                simulationParams['dt'] = np.round((0.2*simulationParams['encounterDistance'])**2/(2*simulationParams['diffusionConstant']),decimals=4)-0.0001                
                
                simulationParams['numMaxSteps'] = int(60//simulationParams['dt'])
                
                print("ε = %s" % (simulationParams['encounterDistance']))
                
                p0 = RCLPolymer(**polymerParams)
                
                simulationParams['waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt_relax']).astype(int)
                
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,"EncounterSimulation")

                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})


                
def watchOneSimulation(polymerParams, simulationParams):
    p0 = RCLPolymer(**polymerParams)
    results = {}
    mc = Experiment(p0, results, simulationParams, "watchEncounter")
    print(mc.results['FET'])
    return mc

def trackAlpha_v_gNc(polymerParams,simulationParams,x_g,x_Nc):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_trackMSD' + '.csv'

    first_time = True
    N = polymerParams['numMonomers']
    
    with open('results/'+filename, 'w') as csvfile:

        for i, g in enumerate(x_g):
  
            print("Simulation for g = %s " % g)
            simulationParams['genomicDistance'] = g

            for j, nc in enumerate(x_Nc):    
                print("Simulation for Nc = %s " % nc)
                polymerParams['Nc'] = nc
                                
                p0 = RCLPolymer(**polymerParams)
                
                simulationParams['waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt_relax']).astype(int)
                
                results = {**polymerParams, **simulationParams}
                mc = Experiment(p0, results, simulationParams,"trackMSD")

                if first_time:
                    fieldnames = ['experimentSetID']+list(mc.results)
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    first_time = False
                writer.writerow({**{'experimentSetID' : str(i)+'_'+str(j)}, **mc.results})


###############################################################################
## Experiments with TADs ######################################################
###############################################################################



################################################################################
################################################################################
########################### MAIN ###############################################
################################################################################
################################################################################
################################################################################    

if __name__ == "__main__":
        
#    proba_vs_keepDFCL(polymerParams,simulationParams,x_Nc,errorbars)   ###########
#    proba_vs_genomicDistance(polymerParams,simulationParams,gmax,gStep,errorbars)
#    proba_vs_Nc_andKeepCL(polymerParams,simulationParams,x_Nc,errorbars)
#    proba_v_sigma(polymerParams,simulationParams,x_sigma,errorbars)
#    proba_v_VEkappa(polymerParams,simulationParams,x_sigma,x_kappa,errorbars)
#    proba_v_gNc(polymerParams,simulationParams,x_g,x_Nc,errorbars)
#    FET_Simulation(polymerParams,simulationParams)
#    mFET_vs_NcinDF(polymerParams,simulationParams,x_Nc,errorbars)
    mc = watchOneSimulation(polymerParams, simulationParams)
    ani = mc.plot_trajectoire(show=True)
    
    plt.figure()
    plt.plot(mc.results['realtime'],mc.results['a1MSD'])
    plt.plot(mc.results['realtime'],mc.results['a2MSD'])
    plt.plot(mc.results['realtime'],mc.results['b1MSD'])
    plt.plot(mc.results['realtime'],mc.results['b2MSD'])
    plt.plot(mc.results['realtime'],mc.results['polymerMSD'])
    plt.show()
    
#    print("hey")
#    trackAlpha_v_gNc(polymerParams,simulationParams,x_g,x_Nc)