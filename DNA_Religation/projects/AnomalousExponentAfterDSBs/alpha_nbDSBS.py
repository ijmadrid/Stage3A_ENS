# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:18:44 2018

@author: ignacio
"""

#############################################################################

## To import the Polymer and Experiments modules
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from modules.Polymers import RCLPolymer
from modules.Experiment import Experiment
from modules.Forces import RepairSphere

import numpy as np
from time import strftime
import csv

import matplotlib.pyplot as plt

def tMSDestimator(trajectory,lag):
    """
    MSD AC estimator
    Return an estimation trajectory-wise of MSD(lag)
    """
    Np = len(trajectory)
    dr = trajectory[0:Np-lag] - trajectory[lag:Np]
    return np.mean(np.linalg.norm(dr, axis = 1)**2)
    

def trackAlpha(xp):
    """
    Track the anomalous exponent of each monomer from a RCL polymer with two 
    defined DSB at positions A1 and B1
    """
    
    # Empty arrays to save results
    MSDfit_A = np.zeros(xp.polymer.numMonomers)
    MSDfit_alpha = np.zeros(xp.polymer.numMonomers)
#    ACMSDfit_A = np.zeros(xp.polymer.numMonomers)
#    ACMSDfit_alpha = np.zeros(xp.polymer.numMonomers)
    monomerSDs = np.zeros((xp.numRealisations, xp.polymer.numMonomers, xp.numSteps))
#    monomerACs = monomerSDs.copy()
    
    removedCLs = 0
    
    for i in range(xp.numRealisations):
        
        #Relax and cut
        removedCLs += xp.makeDefinedBreak(xp.A1, xp.B1)
        
        #Prepare array to save sds of each monomer
        sd = np.zeros((xp.numSteps, xp.polymer.numMonomers))    

        if xp.excludedVolumeCutOff > 0:
            # ADD EXCLUDED VOLUME
            kappa = 3*xp.diffusionConstant/(xp.polymer.b**2)
            
            repulsionForce = lambda polymer : - kappa * RepairSphere(polymer, xp.polymer.freeMonomers, xp.excludedVolumeCutOff)   
            xp.polymer.addnewForce(repulsionForce)

            # Wait some more time (relaxation)
            xp.polymer.step(xp.VE_waitingSteps,xp.dt,xp.diffusionConstant)
                
            
        #Save initial position
        r0 = xp.polymer.get_r().copy()
        
#        trajectory = np.zeros((xp.numSteps, xp.polymer.numMonomers, 3))
            
        for t in range(xp.numSteps):
            # Square displacement of each monomer at time t
            sd[t] = np.linalg.norm( xp.polymer.get_r() - r0 , axis = 1)**2
            # Save position
#            trajectory[t] = xp.polymer.get_r().copy()
            # Position actualization
            xp.polymer.step(1,xp.dt,xp.diffusionConstant)            

        print("\rMotion simulation has finished.",end='\r')
        # SD of each monomer
        monomerSDs[i] = sd.transpose() # N x numSteps
        
        # SD via AC
#        trajectory = trajectory.swapaxes(0,1)
#        for m, trajectory_m in enumerate(trajectory):
##            print(m)
#            for t in range(xp.numSteps):
##                print("AC for t = "+str(t),end='\r')
#                monomerACs[i][m][t] = tMSDestimator(trajectory_m, t)
        
        print("\rSimulation "+str(i)+" of "+str(xp.numRealisations)+".............................",end="\r")
        
        # Use new polymer for the next one
        xp.polymer = xp.polymer.new()
            
    xp.addResults('removedCLs',removedCLs/xp.numRealisations)
    
    # Average over realisations
    monomersMSD = np.mean(monomerSDs, axis = 0)    # size: N x numSteps
#    monomersAC = np.mean(monomerACs, axis = 0)

    xp.addResults('monomersMSD',monomersMSD)
    
#    xp.addResults('monomersAC',monomersAC)

    # Fit and extraction of anomalous exponent
    print("\rFitting....................................................\r")
    realtime = np.arange(xp.numSteps) * xp.dt

    # A and alpha of each monomer
    for m in range(xp.polymer.numMonomers):
        MSDfit_A[m], MSDfit_alpha[m] = A, alpha = xp.getMSDfit(realtime,monomersMSD[m])
#        ACMSDfit_A[m], ACMSDfit_alpha[m] = A_ac, alpha_ac = xp.getMSDfit(realtime,monomersAC[m])
        xp.addResults('A.'+str(m), A)
        xp.addResults('alpha.'+str(m), alpha)
#        xp.addResults('ac_A.'+str(m), A_ac)
#        xp.addResults('ac_alpha.'+str(m), alpha_ac)
    
    print("\rFitting ready!.............................................\r")
    
    xp.addResults('A_monomers', MSDfit_A)
    xp.addResults('alpha_monomers', MSDfit_alpha)
#    
#    xp.addResults('ac_A_monomers', ACMSDfit_A)
#    xp.addResults('ac_alpha_monomers', ACMSDfit_alpha)



def trackAlpha_vsNc(x_Nc):
    date = strftime("%Y_%m_%d_%H_%M")
    filename = date + '_trackAlphaAfterDSB_vsNc'

    with open('results/'+filename+'.csv', 'w') as csvfile:
        
        first_time = True
        
        for i, Nc in enumerate(x_Nc):
            
            polymerParams['Nc'] = Nc
            p0 = RCLPolymer(**polymerParams)
            
            # Waiting time will be relaxation time (Nc dependant)
            simulationParams['waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt_relax']).astype(int)
            simulationParams['VE_waitingSteps'] = np.ceil(p0.relaxTime(simulationParams['diffusionConstant'])/simulationParams['dt']).astype(int)


            results = {**polymerParams, **simulationParams}
            mc = Experiment(p0, results, simulationParams, trackAlpha)
        
            if first_time:
                fieldnames = ['experimentSetID']+list(mc.results)
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                first_time = False
                os.makedirs('results/'+filename+'_figures')
            
            writer.writerow({**{'experimentSetID' : str(i)}, **mc.results})

                # Save MSD curve
            plt.figure()
            plt.plot(mc.results['monomersMSD'].transpose())
            legend = ['_o']*polymerParams['numMonomers']
            legend[simulationParams['A1']] = 'a1'
            legend[simulationParams['A1']+1] = 'a2'
            legend[simulationParams['B1']] = 'b1'
            legend[simulationParams['B1']+1] = 'b1'
            legend = tuple(legend)
            plt.legend(legend)
            plt.savefig('results/'+filename+'_figures/nc'+str(Nc)+'.png')
            plt.close()


############################################################################
############################################################################


############################################################################
# PARAMETERS ###############################################################
############################################################################

                     
polymerParams = dict(numMonomers = 100,
                     dim         = 3,
                     b           = 0.2,
                     Nc          = 5,
                     keepCL      = False
                     )

simulationParams = dict(# Physicial parameters
                        diffusionConstant = 0.008,
                        # Numerical parameters
                        numRealisations   = 500, 
                        dt                = 0.005,
                        dt_relax          = 0.01,
                        excludedVolumeCutOff = 0.1,
                        numSteps = 10000,
#                        waitingSteps = "relax",
#                        encounterDistance = 0.10,
#                        genomicDistance = 10,
                        Nb = 2,
                        Nc_inDamageFoci = 0,
                        A1 = 46,
                        B1 = 53
                        )

###########################################################################
###########################################################################

if __name__ == "__main__":
    
#    date = strftime("%Y_%m_%d_%H_%M")
#    filename = date + '_trackAlphaAfterDSB' + '.csv'
#
#    with open('results/'+filename, 'w') as csvfile:
#        results = {**polymerParams, **simulationParams}
#        p0 = RCLPolymer(**polymerParams)
#        mc = Experiment(p0, results, simulationParams, trackAlpha)
#        
#        fieldnames = ['experimentSetID']+list(mc.results)
#        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#        writer.writeheader()
#        writer.writerow({**{'experimentSetID' : 1}, **mc.results})

    import warnings
    warnings.filterwarnings("ignore")
    
    x_Nc = np.arange(2,60,2)
#    trackAlpha_vsNc(x_Nc)

    # DEFINE THE BREAKS
    polymerParams['numMonomers'] = 100
    simulationParams['A1'] = 30
    simulationParams['B1'] = 68
    
    ### Impose CLs in the DF ...
    simulationParams['Nc_inDamageFoci'] = 3
    
#     ... and keep them
    polymerParams['keepCL'] = False
    trackAlpha_vsNc(x_Nc)
    
    # ... and remove them
#    polymerParams['keepCL'] = False
#    trackAlpha_vsNc(x_Nc)

#    simulationParams['A1'] = 10
#    simulationParams['B1'] = 25
#
#    # ... and keep them
#    polymerParams['keepCL'] = True
#    trackAlpha_vsNc(x_Nc)
#    
#    # ... and remove them
#    polymerParams['keepCL'] = False
#    trackAlpha_vsNc(x_Nc)
    
    
#    ### With selective VE and keeping CLs
#    simulationParams['excludedVolumeCutOff'] = 0.1
#    polymerParams['keepCL'] = True
#    trackAlpha_vsNc(x_Nc)
#    
#    # selective VE and remove CLs
#    polymerParams['keepCL'] = False
#    trackAlpha_vsNc(x_Nc)    
#    
#    ### With selective VE and keeping CLs
#    simulationParams['excludedVolumeCutOff'] = 0.2
#    polymerParams['keepCL'] = True
#    trackAlpha_vsNc(x_Nc)
#    
#    # selective VE and remove CLs
#    polymerParams['keepCL'] = False
#    trackAlpha_vsNc(x_Nc)