# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 18:20:51 2018

@author: ignacio
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from .Forces import ExcludedVolume, LocalExcludedVolume, RepairSphere
from time import time

class Experiment():
    """
    Monte-Carlo set of experiments of the polymer dynamics
    Initate with:
        polymer         :   Polymer object which contains the relevant polymer 
                            parameters (number of monomers, standard deviation 
                            of distance between monomers, number of random 
                            cross-links).
        resultsDataFrame:   Dictionary which will contains all the relevant
                            measures
        params          :   Dictionary with the experiment params:
                                diffusionConstant
                                etc.
    """
    
    def __init__(self, polymer, resultsDataFrame, params, customExperiment = "SimpleDynamicSimulation"):
        # Polymer to simulate
        self.polymer = polymer
        self.trajectoire = []
        # Results data file
        self.results = resultsDataFrame
        self.params = params
        
        self.monomers2follow = []
        
        # Unpackage params values
        for paramKeys, paramValues in params.items(): 
            exec('self.'+paramKeys+'=paramValues')

#        assert np.sqrt(2*self.diffusionConstant*self.dt) <= 0.2*self.encounterDistance, "dt should verify: sqrt(2*D*dt) < 0.2*epsilon"
        
        if customExperiment == "SimpleDynamicSimulation":
            print("Simulation of a Simple Dynamic (no breaks)")
            self.runSimpleDynamic()
        elif customExperiment == "EncounterSimulation":
            print("Simulation of Encounter after two DSBs")
            self.runEncounterSimulation()
        elif customExperiment == "twoDSB":
            print("Simultation of two DSBs until relaxation time")
            self.runTwoRandomBreakSimulation()
        elif customExperiment == "break":
            print('Break and wait')
            self.randomBreaks_SingleStep()
        elif customExperiment == "Encounter_withExcludedVolume":
            print("Simulation of Encounter after two DSBs adding exclusion forces")
            self.BreakAndExclude()
        elif customExperiment == "Encounter_withRepairSphere":
            print("Simulation of Encounter after two DSBs adding selective exclusion forces")
            self.SelectiveExclude()
        elif customExperiment == "TAD_Repair":
            print("Simulation of TAD Repair with one DSB in each TAD")
            self.TAD_repair()
        elif customExperiment == "oneTAD_Repair":
            print("Simulation of the repair of a TAD with two DSB")
            self.oneTADrepair()
        elif customExperiment == 'persistentDSB':
            print("Simulation of two persistent DSB motion")
            self.persistentDSB()
        elif customExperiment == 'watchEncounter':
            self.watchEncounter()
            
        elif customExperiment == 'trackMSD':
            self.MSD_track()
            
        elif customExperiment == 'trackDSB':
            self.trackDSB_distanceDistribution()
        
        elif customExperiment == 'measureMSRG':
            self.measureStatisticalProperties()
        
        elif callable(customExperiment):
            # customExperiment should be a function that has a 
            # Experiment object as unique argument
            print("Custom experiment")
            customExperiment(self)
        
        else:
            print("Custom Experiment should be a function or the name of an implemented experiment!!!")
    
    def get_params(self):
        return self.params
    
    def addResults(self, name, value):
        self.results[name] = value

    def plot_trajectoire(self,show=False):
        fig = plt.figure()
        ax = Axes3D(fig)
#        ax.auto_scale_xyz([-10, 10], [-10, 10], [-10, 10])
        
        
        
        X = self.trajectoire[0][:,0]
        Y = self.trajectoire[0][:,1]
        Z = self.trajectoire[0][:,2]
        line, = ax.plot(X, Y, Z)
        dots = ax.scatter(X, Y, Z, c=self.polymer.colors, marker='o', s = 30)
    
#        cls = []    
#        clpairs = self.polymer.offDiagPairs()
#        for clpair in clpairs:
#            mx = X[clpair[0]]
#            nx = X[clpair[1]]
#            my = Y[clpair[0]]
#            ny = Y[clpair[1]]
#            mz = Z[clpair[0]]
#            nz = Z[clpair[1]]
#            cl, = ax.plot([mx,nx],[my,ny],[mz,nz],color = 'g')
#            cls.append(cl)
        
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0

        mid_x = (X.max()+X.min()) * 0.5
        mid_y = (Y.max()+Y.min()) * 0.5
        mid_z = (Z.max()+Z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)


        def animate(i):
            ri = self.trajectoire[i]
            line.set_xdata(ri[:,0])  # update the data
            line.set_ydata(ri[:,1])  # update the data
            line.set_3d_properties(ri[:,2])
            dots._offsets3d = (ri[:,0], ri[:,1], ri[:,2])
#            for p, cl in enumerate(cls):
#                cl.set_xdata(ri[:,0][clpairs[p][0]],ri[:,0][clpairs[p][1]])
#                cl.set_ydata(ri[:,1][clpairs[p][0]],ri[:,1][clpairs[p][1]]) 
#                cl.set_3d_properties(ri[:,2][clpairs[p][0]],ri[:,2][clpairs[p][1]]) 
            return line,
                
        def init():
            line.set_data([],[])
            line.set_3d_properties([])
#            for cl in cls:
#                cl.set_data([],[])
#                cl.set_3d_properties([])
            return line,
        
        ani = animation.FuncAnimation(fig, animate, frames=self.numSteps, init_func=init,
                                              interval=self.dt, blit=False, repeat = True)
    
        plt.show()
        if not show:
            plt.close(fig)
        
        return ani

    
    def addMonomer2follow(self,monomer):
        self.monomers2follow.append(monomer)
    
    def reduceResults(self,experiment2):
        """
        Combine the results of two experiments (self and experiment2)
        """
        assert experiment2.get_params() == self.params, "Experiments don't have same parameters. Results should not be reduced into a single result object."
        results = {}
        results['iterationsNumber'] = self.numRealisations + experiment2.numRealisations
        results['lastTrajectory'] = experiment2.trajectory
        
        return Experiment()
        
        
    ####################################################
    ############# IMPLEMENTED EXPERIMENTS ##############
    ####################################################    
        
    def runSimpleDynamic(self):        
        
        self.addResults("iterationsNumber",self.numRealisations)
        msrg = np.zeros(self.numRealisations)   
        
        for i in range(self.numRealisations):            
            # Burn in until relaxation time
            relaxSteps = np.ceil(self.polymer.relaxTime(self.diffusionConstant)/self.dt_relax).astype(int)
            self.polymer.step(relaxSteps, self.dt_relax, self.diffusionConstant)        
            # Once relaxes calcule some statistical properties
            msrg[i] = self.polymer.get_msrg()
      
            # Main simulation           
            # Save initial position
            trajectory = np.zeros((self.numSteps+1,self.polymer.numMonomers,3))
            trajectory[0] = self.polymer.get_r()
            
            for t in range(self.numSteps):
                self.polymer.step(1, self.dt, self.diffusionConstant)
                trajectory[t+1] = self.polymer.get_r()
            
            # At the end of each iteration make a new polymer
            # for the next one
            self.polymer = self.polymer.new()
        
        # Save results
        self.trajectoire = trajectory
        self.addResults("lastTrajectory", trajectory)
        self.addResults("MSRG", np.mean(msrg))
   

    def makeDefinedBreak(self, A1, B1):
        breakLoci = (A1,B1)

        # Add CLs to the Damage-to-be foci
        if(self.Nc_inDamageFoci > 0):
            self.polymer.imposeDSBconnections(self.Nc_inDamageFoci,breakLoci)
        
        # Verify is polymer is splittable for the prepeared DSBs
        trialNb = 0
        
        while(not(self.polymer.isSplittable(breakLoci))):
            # if not, make new connections (#TODO: try an heuristic maybe?)
            trialNb +=  1
#            print("Try",trialNb,"for making a valid polymer.")
            self.polymer.reset()
            if(self.Nc_inDamageFoci > 0):
                self.polymer.imposeDSBconnections(self.Nc_inDamageFoci,breakLoci)
        
#        self.results["TotalRejectedPolymers"] += trialNb
        # Once the polymer is splittable:
        # Burn in until relaxation time
        relaxSteps = np.ceil(self.polymer.relaxTime(self.diffusionConstant)/self.dt_relax).astype(int)
#        self.addResults("relaxSteps", relaxSteps)
        self.polymer.step(relaxSteps,self.dt_relax,self.diffusionConstant)
        
        # Induce DSBs
        self.polymer.cutNow(breakLoci,definitive=True)
        # Remove the CL concerning the cleavages if it is the case
        if self.polymer.keepCL == False:
            removedNum = self.polymer.removeCL()
        else:
            removedNum = 0
            
        # Wait some more time
        self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)
        
        return removedNum
        
        
        
    def randomBreaks_SingleStep(self):
        
        # Prepare the random DSBs
        breakLoci = self.polymer.randomCuts(self.genomicDistance,self.Nb)
        
        # Add CLs to the Damage-to-be foci
        if(self.Nc_inDamageFoci > 0):
            self.polymer.imposeDSBconnections(self.Nc_inDamageFoci,breakLoci)
        
        # Verify is polymer is splittable for the prepeared DSBs
        trialNb = 0
        
        while(not(self.polymer.isSplittable(breakLoci))):
            # if not, make new connections (#TODO: try an heuristic maybe?)
            trialNb +=  1
#            print("Try",trialNb,"for making a valid polymer.")
            self.polymer.reset()
            if(self.Nc_inDamageFoci > 0):
                self.polymer.imposeDSBconnections(self.Nc_inDamageFoci,breakLoci)
        
#        self.results["TotalRejectedPolymers"] += trialNb
        # Once the polymer is splittable:
        # Burn in until relaxation time
        relaxSteps = np.ceil(self.polymer.relaxTime(self.diffusionConstant)/self.dt_relax).astype(int)
#        self.addResults("relaxSteps", relaxSteps)
        self.polymer.step(relaxSteps,self.dt_relax,self.diffusionConstant)
        
        # Induce DSBs
        self.polymer.cutNow(breakLoci,definitive=True)
        # Remove the CL concerning the cleavages if it is the case
        if self.polymer.keepCL == False:
            removedNum = self.polymer.removeCL()
        else:
            removedNum = 0
            
        # Wait some more time
        self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)
        
        return removedNum
    
    
    
    def runTwoRandomBreakSimulation(self):
        
        self.addResults("iterationsNumber",self.numRealisations)
        msrg = np.zeros(self.numRealisations)   
        
        for i in range(self.numRealisations): 
            self.randomBreaks_SingleStep()
            msrg[i] = self.polymer.get_msrg()
            self.polymer = self.polymer.new()
            
        self.addResults("MSRG", np.mean(msrg))
        halfCI = 1.96*np.std(msrg)/np.sqrt(len(msrg))
        self.addResults("MSRG_95CI", halfCI)
        
    
    
    def runEncounterSimulation(self):
        self.addResults("TotalRejectedPolymers", 0)
        self.addResults("iterationsNumber",self.numRealisations)
        FETs = np.ones(self.numRealisations)*np.nan
        events = np.repeat('NA',self.numRealisations).astype('<U15')
        removedNums = np.ones(self.numRealisations)*np.nan
        post_msrgs = np.ones(self.numRealisations)*np.nan
        
        for i in range(self.numRealisations):
                        
            # Simulates the break and some waiting time
            removedNum = self.randomBreaks_SingleStep()
            
            ##################################################                 
            # Simulation until encounter              
            t = 0
            didEncounter = self.polymer.anyEncountered(self.encounterDistance)
            while(not(didEncounter[0]) and t < self.numMaxSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                didEncounter = self.polymer.anyEncountered(self.encounterDistance)
                t += 1
                
            if( t <self.numMaxSteps):
                FETs[i] = t 
                events[i] = didEncounter[1]
                removedNums[i] = removedNum
                msrg_post_encounter = self.polymer.get_msrg() 
                post_msrgs[i] = msrg_post_encounter
        
            self.polymer = self.polymer.new()
            ##################################################  

        self.saveEncounterResults(FETs, events, removedNums, post_msrgs)
#        self.addResults('MSRG_atEncounter', post_msrgs)
#        self.addResults('Ensemble_MSRG', np.nanmean(post_msrgs))
    
    
    
    def BreakAndExclude(self):
        """
        After inducing the DSBs add exclusion forces (excluded volume forces)
        from the cleaved monomers.
        """
        self.addResults("TotalRejectedPolymers", 0)
        self.addResults("iterationsNumber",self.numRealisations)
        FETs = np.ones(self.numRealisations)*np.nan
        events = np.repeat('NA',self.numRealisations).astype('<U15')
        removedNums = np.ones(self.numRealisations)*np.nan
        post_msrgs = np.ones(self.numRealisations)*np.nan

        if self.numRealisations >= 10: k = self.numRealisations//10        
        else: k = 1
        for i in range(self.numRealisations):
            
            if not(i%k):
                print("|",'='*(i//k),'-'*(10-i//k),"| Simulation", i+1, "of", self.numRealisations)

            # Simulates the break and some waiting time:
            removedNum = self.randomBreaks_SingleStep()

            # ADD EXCLUDED VOLUME       
            kappa = 3*self.diffusionConstant/(self.polymer.b**2)
            repulsionForce = lambda polymer : - kappa * LocalExcludedVolume(polymer, self.polymer.freeMonomers, self.excludedVolumeCutOff)   
            self.polymer.addnewForce(repulsionForce)

            # Wait some more time
            self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)

            ##################################################                 
            # Simulation until encounter              
            t = 0
            didEncounter = self.polymer.anyEncountered(self.encounterDistance)
            while(not(didEncounter[0]) and t < self.numMaxSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                didEncounter = self.polymer.anyEncountered(self.encounterDistance)
                t += 1
                
            if( t <self.numMaxSteps):
                FETs[i] = t 
                events[i] = didEncounter[1]
                removedNums[i] = removedNum
                msrg_post_encounter = self.polymer.get_msrg() 
                post_msrgs[i] = msrg_post_encounter
        
            self.polymer = self.polymer.new()
            ##################################################  

        self.saveEncounterResults(FETs, events, removedNums, post_msrgs)
#        self.addResults('MSRG_atEncounter', post_msrgs)
#        self.addResults('Ensemble_MSRG', np.nanmean(post_msrgs))


    def SelectiveExclude(self):
        """
        Add a repair spehere, repulsive for the non-cleaved monomers
        """
        self.addResults("TotalRejectedPolymers", 0)
                
        self.addResults("iterationsNumber",self.numRealisations)
        FETs = np.ones(self.numRealisations)*np.nan
        events = np.repeat('NA',self.numRealisations).astype('<U15')
        removedNums = np.ones(self.numRealisations)*np.nan
        post_msrgs = np.ones(self.numRealisations)*np.nan

        if self.numRealisations >= 10: k = self.numRealisations//10        
        else: k = 1
        
        for i in range(self.numRealisations):
            
            if not(i%k):
                print("|",'='*(i//k),'-'*(10-i//k),"| Simulation", i+1, "of", self.numRealisations)

            # Simulates the break and some waiting time:
            removedNum = self.randomBreaks_SingleStep()

            # ADD EXCLUDED VOLUME
            if(self.excludedVolumeCutOff > 0):
                try:
                    kappa = self.excludedVolumeSpringConstant
                except:
                    kappa = 3*self.diffusionConstant/(self.polymer.b**2)
                
                repulsionForce = lambda polymer : - kappa * RepairSphere(polymer, self.polymer.freeMonomers, self.excludedVolumeCutOff)   
                self.polymer.addnewForce(repulsionForce)
    
                # Wait some more time
                self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)

            ##################################################                 
            # Simulation until encounter              
            t = 0
            didEncounter = self.polymer.anyEncountered(self.encounterDistance)
            while(not(didEncounter[0]) and t < self.numMaxSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                didEncounter = self.polymer.anyEncountered(self.encounterDistance)
                t += 1
                
            if( t <self.numMaxSteps):
                FETs[i] = t 
                events[i] = didEncounter[1]
                removedNums[i] = removedNum
                msrg_post_encounter = self.polymer.get_msrg() 
                post_msrgs[i] = msrg_post_encounter
        
            self.polymer = self.polymer.new()
            ##################################################  

        self.saveEncounterResults(FETs, events, removedNums, post_msrgs)
        

    def TAD_repair(self):
        """
        One independant break in each TAD
        """
        
        self.addResults("iterationsNumber",self.numRealisations)
        FETs = np.ones(self.numRealisations)*np.nan
        events = np.repeat('NA',self.numRealisations).astype('<U15')
        removedNums = np.ones(self.numRealisations)*np.nan
        
        for i in range(self.numRealisations):
            
            
            ##################################################            
            # Prepare the random DSBs
            breakLoci = []
            for TAD in self.polymer.TADs:
                breakLoci.append(self.polymer.inregionCut(TAD))
                
            # Verify is polymer is splittable for the prepeared DSBs
            while(not(self.polymer.isSplittable(breakLoci))):
                # if not, make new connections (#TODO: try an heuristic maybe?)
                self.polymer.reset()
            
            # Once the polymer is splittable:
            # Burn in until relaxation time
            relaxSteps = np.ceil(self.polymer.relaxTime(self.diffusionConstant)/self.dt_relax).astype(int)
            self.polymer.step(relaxSteps,self.dt_relax,self.diffusionConstant)
            
            # Induce DSBs
            self.polymer.cutNow(breakLoci,definitive=True)
            # Remove the CL concerning the cleavages if it is the case
            if self.polymer.keepCL == False:
                removedNum = self.polymer.removeCL()
        
            # Wait some more time
            self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)
            ##################################################  
        
            ##################################################                 
            # Simulation until encounter              
            t = 0
            didEncounter = self.polymer.anyEncountered(self.encounterDistance)
            while(not(didEncounter[0]) and t < self.numMaxSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                didEncounter = self.polymer.anyEncountered(self.encounterDistance)
                t += 1
                
            if( t <self.numMaxSteps):
                FETs[i] = t 
                events[i] = didEncounter[1]
                removedNums[i] = removedNum
        
            self.polymer = self.polymer.new()
            ##################################################  
        
        self.saveEncounterResults(FETs, events, removedNums)
        
    
    def oneTADrepair(self):
        """
        Induce Nb DSBs separated by a given genomic distance on the first TAD
        """

        self.addResults("iterationsNumber",self.numRealisations)
        FETs = np.ones(self.numRealisations)*np.nan
        events = np.repeat('NA',self.numRealisations).astype('<U15')
        removedNums = np.ones(self.numRealisations)*np.nan
        
        boundaries = np.append(0, self.polymer.subDomainnumMonomers.cumsum())
        k = self.selectedSubDomain
        l = boundaries[k]
        r = boundaries[k+1]

        for i in range(self.numRealisations):
            
            starttime = time()
            
            ##################################################            
            # Prepare the random DSBs
            breakLoci = self.polymer.inregionCut(l,r,self.genomicDistance,self.Nb)
            
            # Verify is polymer is splittable for the prepeared DSBs
            while(not(self.polymer.isSplittable(breakLoci))):
                # if not, make new connections (#TODO: try an heuristic maybe?)
                self.polymer.reset()
            
            # Once the polymer is splittable:
            # Burn in until relaxation time
            relaxSteps = np.ceil(self.polymer.relaxTime(self.diffusionConstant)/self.dt_relax).astype(int)
            self.polymer.step(relaxSteps,self.dt_relax,self.diffusionConstant)
            
            # Induce DSBs
            self.polymer.cutNow(breakLoci,definitive=True)
            # Remove the CL concerning the cleavages if it is the case
            if self.polymer.keepCL == False:
                removedNum = self.polymer.removeCL()
            else:
                removedNum = 0
        
            # Wait some more time
            self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)
            ##################################################  
        
            ##################################################                 
            # Simulation until encounter              
            t = 0
            didEncounter = self.polymer.anyEncountered(self.encounterDistance)
            while(not(didEncounter[0]) and t < self.numMaxSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                didEncounter = self.polymer.anyEncountered(self.encounterDistance)
                t += 1
                
            if( t <self.numMaxSteps):
                FETs[i] = t  
                events[i] = didEncounter[1]
                removedNums[i] = removedNum
        
            self.polymer = self.polymer.new()
            ##################################################  
            
            print('\r' + 'Simulation ' + str(i) + ') DSB at '+str(breakLoci)+' | Execution time : ' + str(time() - starttime), end='\r')
            
        self.saveEncounterResults(FETs, events, removedNums)
    
    

    def persistentDSB(self):
        """
        Simulate the motion of persistent DSBs, one in each TAD
        """
        
        assert type(self.polymer.subDomainnumMonomers) == np.ndarray, "The polymer must contain sub-domains"
        
        self.addResults("iterationsNumber",self.numRealisations)
        
        
        msrg = np.zeros(self.numRealisations)   
        
        from scipy.special import comb
        encounternumber = comb(2*len(self.polymer.subDomainnumMonomers),2).astype(int)
        interbreakDistance = np.zeros((self.numRealisations,self.numSteps+1,encounternumber))


        boundaries = np.append(0, self.polymer.subDomainnumMonomers.cumsum())
        leftEnds = boundaries[:-1]
        rightEnds = boundaries[1:]
        
        for i in range(self.numRealisations):
            
            ##################################################            
            # Prepare the random DSBs
            breakLoci = np.zeros(len(self.polymer.subDomainnumMonomers)).astype(int)
            for i in range(len(self.polymer.subDomainnumMonomers)):
                breakLoci[i] = self.polymer.inregionCut(leftEnds[i],rightEnds[i],1,1)
                
            # Verify is polymer is splittable for the prepeared DSBs
            while(not(self.polymer.isSplittable(breakLoci))):
                # if not, make new connections (#TODO: try an heuristic maybe?)
                self.polymer.reset()
            
            # Once the polymer is splittable:
            # Burn in until relaxation time
            relaxSteps = np.ceil(self.polymer.relaxTime(self.diffusionConstant)/self.dt_relax).astype(int)
            self.polymer.step(relaxSteps,self.dt_relax,self.diffusionConstant)
            
            # Induce DSBs
            self.polymer.cutNow(breakLoci,definitive=True)
            # Remove the CL concerning the cleavages if it is the case
            if self.polymer.keepCL == False:
                self.polymer.removeCL()

            # ADD EXCLUDED VOLUME
            if self.excludedVolumeCutOff > 0:
                kappa = 3*self.diffusionConstant/(self.polymer.b**2)
                repulsionForce = lambda polymer : - kappa * LocalExcludedVolume(polymer, self.polymer.freeMonomers, self.excludedVolumeCutOff)   
                self.polymer.addnewForce(repulsionForce)
            
            # Wait some more time
#            self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)
            ##################################################  
            
            # Once relaxes calcule some statistical properties
            msrg[i] = self.polymer.get_msrg()
      
            # Main simulation           
            interbreakDistance[i][0] = self.polymer.interBreakDistance()
            
            for t in range(self.numSteps):
                self.polymer.step(1, self.dt, self.diffusionConstant)
                interbreakDistance[i][t+1] = self.polymer.interBreakDistance()
                
            self.polymer = self.polymer.new()
            ##################################################  
        
        self.addResults('MSRG_atEncounter', msrg)
        self.addResults('Ensemble_MSRG', np.mean(msrg))
        self.addResults('Ensemble_MSRG_dx', 1.96*np.std(msrg)/np.sqrt(len(msrg)))
        self.addResults('MeanInterBreakDistance',np.mean(interbreakDistance,axis=0))


    def computeSD(self,t):
        """
        Square Deviation of each monomer at time t
        """
        r0 = self.trajectoire[0]
        rt = self.trajectoire[t]
        return np.linalg.norm( rt - r0 , axis = 1)**2
        
    def get_avg_msd(self):
        return np.mean(np.mean(self.sds,axis=0),axis=1)
 
    def get_msd_per_monomer(self):
        return np.mean(self.sds,axis=0)
    
    def MSD_track(self):
        """
        Track of statistical properties of the breaks troughout time
        Goal: Extract the anomalous exponent of each monomer after the break
        """
        
        self.addResults("iterationsNumber",self.numRealisations)
        polymerSDs = np.zeros((self.numRealisations,self.numSteps))
        breaksSDs = np.zeros((self.numRealisations,4,self.numSteps))
        #events = np.repeat('NA',self.numRealisations).astype('<U15')
        #removedNums = np.ones(self.numRealisations)*np.nan
        #post_msrgs = np.ones(self.numRealisations)*np.nan

        if self.numRealisations >= 10: k = self.numRealisations//10        
        else: k = 1
        
        for i in range(self.numRealisations):
            if not(i%k):
                print("|",'='*(i//k),'-'*(10-i//k),"| Simulation", i+1, "of", self.numRealisations)

            # Simulates the break and some waiting time:
            removedNum = self.randomBreaks_SingleStep()

            if self.excludedVolumeCutOff > 0:
                # ADD EXCLUDED VOLUME
                try:
                    kappa = self.excludedVolumeSpringConstant
                except:
                    kappa = 3*self.diffusionConstant/(self.polymer.b**2)
                
                repulsionForce = lambda polymer : - kappa * RepairSphere(polymer, self.polymer.freeMonomers, self.excludedVolumeCutOff)   
                self.polymer.addnewForce(repulsionForce)
    
                # Wait some more time
                self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)

            ##################################################                 
            # Simulation until NumSteps
            
            # Will contain the Square Displacement of each monomer
            sd = np.zeros((self.numSteps, self.polymer.numMonomers))
            
            r0 = self.polymer.get_r().copy()
            for t in range(self.numSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                # Square displacement of each monomer at time t
                sd[t] = np.linalg.norm( self.polymer.get_r() - r0 , axis = 1)**2


            # Polymer SD
            polymerSDs[i] = np.mean(sd,axis=1)

            # Breaks SD
            #breaksSDs[i] = np.zeros((4,self.numSteps))
            for j, brokenMonomer in enumerate(self.polymer.freeMonomers):
                breaksSDs[i][j] = sd[:,brokenMonomer]
            
                       

            self.polymer = self.polymer.new()
            ##################################################  

#        self.saveEncounterResults(FETs, events, removedNums, post_msrgs)
        
        ### Computation of MSD averaging over all realisations
        polymerMSD = np.mean(polymerSDs, axis = 0) # size: numSteps
        breaksMSD = np.mean(breaksSDs, axis = 0) # size: 4 x numSteps
        
        ### Extraction of anomalous exponent
        
        realtime = np.arange(self.numSteps) * self.dt
        
        polymerMSDfit_amplitude, polymerMSDfit_alpha = self.getMSDfit(realtime, polymerMSD)
        breakMSDfit_amplitude = np.zeros(4)
        breakMSDfit_alpha = np.zeros(4)
        for b in range(4):
            breakMSDfit_amplitude[b], breakMSDfit_alpha[b] = self.getMSDfit(realtime, breaksMSD[b])
        
        self.addResults("polymerMSD",polymerMSD)
        self.addResults("A1_MSD",breaksMSD[0])
        self.addResults("A2_MSD",breaksMSD[1])
        self.addResults("B1_MSD",breaksMSD[2])
        self.addResults("B2_MSD",breaksMSD[3])
        
        self.addResults("polymerAlpha",polymerMSDfit_alpha)
        self.addResults("A1_Alpha",breakMSDfit_alpha[0])
        self.addResults("A2_Alpha",breakMSDfit_alpha[1])
        self.addResults("B1_Alpha",breakMSDfit_alpha[2])
        self.addResults("B2_Alpha",breakMSDfit_alpha[3])
        
        self.addResults("polymer_MSDamplitude",polymerMSDfit_amplitude)
        self.addResults("A1_MSDamplitude",breakMSDfit_amplitude[0])
        self.addResults("A2_MSDamplitude",breakMSDfit_amplitude[1])
        self.addResults("B1_MSDamplitude",breakMSDfit_amplitude[2])
        self.addResults("B2_MSDamplitude",breakMSDfit_amplitude[3])
        
        

    def exponentialFit(self,x):
        from scipy.stats import expon
        
        try:
            loc, scale = expon.fit(x)
        
            amplitude = np.exp(loc/scale)/scale
            rate      = 1/scale
        
        except:
            amplitude = np.nan
            rate      = np.nan
            
        popt = (amplitude, rate)
        
#        from scipy.optimize import curve_fit
#        
#        def exponential(x, amplitude, lambd):
#            return amplitude * np.exp(-lambd*x)
#        
#        bin_heights, bin_borders = np.histogram(x, normed='True')
#        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
#        try:
#            popt, _ = curve_fit(exponential, bin_centers, bin_heights, p0=[x[0],1/np.mean(x)])
#        except:
#            popt = (np.nan, np.nan)
            
        return popt
    
    
    def powerFit(self, x):
        from scipy.optimize import curve_fit
        
        def power(x, A, mu):
            return A * x ** ( -1 - mu)
        
        bin_heights, bin_borders = np.histogram(x, normed='True')
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        try:
            popt, _ = curve_fit(power, bin_centers, bin_heights, p0=[x[0],1/np.mean(x)])
        except:
            popt = (np.nan, np.nan)
        return popt
    
    
    def getMSDfit(self,time,MSDarray):
        from scipy.optimize import curve_fit
        
        def power(t, A, alpha):
            return A*t**alpha
        
        try:
            popt, _ = curve_fit(power, time, MSDarray)
        except:
            popt = (np.nan, np.nan)
        
        return popt
        #TODO
        
#        x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
#        plt.plot(x_interval_for_fit, exponential(x_interval_for_fit, *popt), label='Fit %s exp(- %s t)' % (*popt))
#        plt.legend()



    
    
    def trackDSB_distanceDistribution(self):
        """
        Get histogram of times where broken monomers and near and/or far from each other
        """

        self.addResults("iterationsNumber",self.numRealisations)
        polymerSDs = np.zeros((self.numRealisations,self.numSteps))
        centerSDs = np.zeros((self.numRealisations,self.numSteps))
        breaksSDs = np.zeros((self.numRealisations,2*self.Nb,self.numSteps))
#
#        from scipy.special import comb
#        encounternumber = comb(2*self.Nb,2).astype(int)
        interbreakDistance = np.zeros((self.numRealisations,self.numSteps,6))
        
#        if self.numRealisations >= 10: k = self.numRealisations//10        
#        else: k = 1
        
        start = time()
        for i in range(self.numRealisations):
#            if not(i%k):
#                print("|",'='*(i//k),'-'*(10-i//k),"| Simulation", i+1, "of", self.numRealisations)


            # Simulates the break and some waiting time:
            
            try:
                self.makeDefinedBreak(self.A1, self.B1)
                self.genomicDistance = self.B1 - self.A1 - 1
            except:
                self.randomBreaks_SingleStep()

            if self.excludedVolumeCutOff > 0:
                # ADD EXCLUDED VOLUME
                try:
                    kappa = self.excludedVolumeSpringConstant
                except:
                    kappa = 3*self.diffusionConstant/(self.polymer.b**2)
                
                repulsionForce = lambda polymer : - kappa * RepairSphere(polymer, self.polymer.freeMonomers, self.excludedVolumeCutOff)   
                self.polymer.addnewForce(repulsionForce)
    
                # Wait some more time
                self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)

            ##################################################                 
            # Simulation until NumSteps
            
            # Will contain the Square Displacement of each monomer
            sd = np.zeros((self.numSteps, self.polymer.numMonomers))
            # Will contain the SD of the center of mass
            sdCenter = np.zeros(self.numSteps)
            
            # Will contain the distance between breaks
            
            
            r0 = self.polymer.get_r().copy()
            rc0 = self.polymer.getCenterOfMass().copy()
            for t in range(self.numSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                # Center of mass
                sdCenter[t] = np.linalg.norm( self.polymer.getCenterOfMass() - rc0 )**2
                # Square displacement of each monomer at time t
                sd[t] = np.linalg.norm( self.polymer.get_r() - r0 , axis = 1)**2
                # Inter break distance
                interbreakDistance[i][t] = self.polymer.interBreakDistance()

            # Polymer SD
            polymerSDs[i] = np.mean(sd,axis=1)

            # Center of Mass SD
            centerSDs[i] = sdCenter
            
            # Breaks SD
            breaksSDs[i] = np.zeros((4,self.numSteps))
            for j, brokenMonomer in enumerate(self.polymer.freeMonomers):
                breaksSDs[i][j] = sd[:,brokenMonomer]
            
                       

            self.polymer = self.polymer.new()
            
            transcurredTime = time() - start
            totalTime = transcurredTime * self.numRealisations/(i+1)
            print('\r' + 'Simulation ' + str(i) + ' of ' + str(self.numRealisations) + ' | ' +
                  str(round(totalTime - transcurredTime, 2)) + ' sec still.', end='\r')
            
            ##################################################  

        
        ### Save d(m,n)'s
        ### VARIANCE ALONG TIME SERIES
        ### AVERAGE VARIANCE ALONG REALIZATIONS

        ### Get break names
        from itertools import combinations
        monomerNames = [chr(97 + i//2) + str(1 + i%2) for i in range(2*self.Nb)]
        breakNames = [i for i in combinations(monomerNames, r = 2)]
        combinationNames = ['%s-%s' % breakNames[i] for i in range(len(breakNames))]
        
        # Variance and mean over time
        var_t  = np.var(interbreakDistance, axis = 1)   # numRealisations x (2Nb 2)
        mean_t = np.mean(interbreakDistance, axis = 1)   # numRealisations x (2Nb 2)
        
        # Mean variance and variance of means (calculated above) over realisations
        meanvarDistance = np.mean(var_t, axis = 0)  # (2Nb 2) x 1
        varmeanDistance = np.var(mean_t, axis = 0)  # (2Nb 2) x 1
        
        # Total Variance
        totalvarDistance = meanvarDistance + varmeanDistance
        
        for i in range(6):
            self.addResults(combinationNames[i]+"_interbreakDistance_totalVariance", totalvarDistance[i])
            self.addResults(combinationNames[i]+"_interbreakDistance_meanOfVariances", meanvarDistance[i])
#        self.addResults('MeanInterBreakDistance',np.mean(interbreakDistance,axis=0))
        

        ### Get interval lengths
        aboveLengths = {mn:[] for mn in combinationNames}
        belowLengths = {mn:[] for mn in combinationNames}
        for k, key in enumerate(combinationNames):
            for i in range(self.numRealisations):
                up, lo = self.getThresholdedIntervalLengths(interbreakDistance[i,:,k])
                aboveLengths[key].extend(up)
                belowLengths[key].extend(lo)
        
        

        ### Fit Exponential Distribution to the Time Intervals Distribution
        above_fit = {}
        below_fit = {}
        for mn in aboveLengths.keys():
            aboveLengths[mn] = np.array(aboveLengths[mn])*self.dt
            belowLengths[mn] = np.array(belowLengths[mn])*self.dt
            above_fit[mn] = a = self.exponentialFit(aboveLengths[mn])
            below_fit[mn] = b = self.exponentialFit(belowLengths[mn])
            
            # Histogram
            bin_heights, bin_borders = np.histogram(aboveLengths[mn], normed='True')
            bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
            self.addResults(mn+"_HistogramCenters", bin_centers)
            self.addResults(mn+"_HistogramHeights", bin_heights)
            
        
            self.addResults(mn+"_ComebackRate", a[1])
            self.addResults(mn+"_ComebackAmplitude", a[0])
            self.addResults(mn+"_meanComebackTime", aboveLengths[mn].mean())
            self.addResults(mn+"_stdComebackTime", aboveLengths[mn].std())
            self.addResults(mn+"_meanComebackTime_dx", 1.96*aboveLengths[mn].std()/np.sqrt(len(aboveLengths[mn])))
            
            self.addResults(mn+"_TakeoffRate", b[1])
            self.addResults(mn+"_TakeoffAmplitude", b[0])
            self.addResults(mn+"_meanTakeoffTime", belowLengths[mn].mean())
            self.addResults(mn+"_stdTakeoffTime", belowLengths[mn].std())
            self.addResults(mn+"_meanTakeoffTime_dx", 1.96*belowLengths[mn].std()/np.sqrt(len(belowLengths[mn])))
            

#        self.addResults("aboveTimes", aboveLengths)
#        self.addResults("belowTimes", belowLengths)
        
        meanTimes = np.array([np.mean(aboveLengths[i]) for i in aboveLengths.keys()])
        self.addResults("estimatedProba_byMeanTimes", (1/meanTimes[0] + 1/meanTimes[-1])/np.sum(1/meanTimes))
        
        rates = np.array([above_fit[i][1] for i in aboveLengths.keys()])
        self.addResults("estimatedProba_byRates", (rates[0] + rates[-1])/np.nansum(rates))

        
#
#        self.addResults("aboveFit", above_fit)
#        self.addResults("belowFit", below_fit)        
        
#        ### Computation of MSD averaging over all realisations
        polymerMSD = np.mean(polymerSDs, axis = 0) # size: numSteps
        breaksMSD = np.mean(breaksSDs, axis = 0) # size: 2*Nb x numSteps
        centerMSD = np.mean(centerSDs, axis = 0)
#        
        ### Extraction of anomalous exponent
        
        realtime = np.arange(self.numSteps) * self.dt
        
        polymerMSDfit_amplitude, polymerMSDfit_alpha = self.getMSDfit(realtime, polymerMSD)
        centerMSDfit_amplitude, centerMSDfit_alpha = self.getMSDfit(realtime, centerMSD)
        breakMSDfit_amplitude = np.zeros(2*self.Nb)
        breakMSDfit_alpha = np.zeros(2*self.Nb)
        for b in range(2*self.Nb):
            breakMSDfit_amplitude[b], breakMSDfit_alpha[b] = self.getMSDfit(realtime, breaksMSD[b])
        
        self.addResults("polymerMSD",polymerMSD)
        self.addResults("centerMSD",centerMSD)
        for im, m in enumerate(monomerNames):
            self.addResults(m+"_MSD", breaksMSD[im])
        
        self.addResults("polymerAlpha",polymerMSDfit_alpha)
        self.addResults("centerAlpha",centerMSDfit_alpha)
        for im, m in enumerate(monomerNames):
            self.addResults(m+"_Alpha", breakMSDfit_alpha[im])
        
        self.addResults("polymer_MSDamplitude",polymerMSDfit_amplitude)
        self.addResults("center_MSDamplitude",centerMSDfit_amplitude)
        for im, m in enumerate(monomerNames):
            self.addResults(m+"_Amplitude", breakMSDfit_amplitude[im])


    def getThresholdedIntervalLengths(self, distance):
        """
        Return the lengths of time intervals in which 
        |Rm - Rn| > threshold
        distance is 1D
        """
        thresholded = (distance >= self.distanceThreshold)
        edges = np.where(np.diff(thresholded))[0]
        intervalLengths = np.diff(np.append(0,edges+1))
        aboveIntervalLengths = intervalLengths[1+thresholded[0]::2]
        belowIntervalLengths = intervalLengths[2-thresholded[0]::2]
        return (aboveIntervalLengths, belowIntervalLengths)


    def measureStatisticalProperties(self):
        """
        Measure statistical properties at given times
        """
        self.addResults("iterationsNumber",self.numRealisations)
        msrg = np.zeros((self.numRealisations,len(self.times2sample)))
        breaks2center = np.zeros((self.numRealisations,len(self.times2sample),2*self.Nb))
                
        start = time()
        for i in range(self.numRealisations):
#            if not(i%k):
#                print("|",'='*(i//k),'-'*(10-i//k),"| Simulation", i+1, "of", self.numRealisations)


            # Simulates the break and some waiting time:
            
            try:
                self.makeDefinedBreak(self.A1, self.B1)
                self.genomicDistance = self.B1 - self.A1 - 1
            except:
                self.randomBreaks_SingleStep()

            if self.excludedVolumeCutOff > 0:
                # ADD EXCLUDED VOLUME
                try:
                    kappa = self.excludedVolumeSpringConstant
                except:
                    kappa = 3*self.diffusionConstant/(self.polymer.b**2)
                
                repulsionForce = lambda polymer : - kappa * RepairSphere(polymer, self.polymer.freeMonomers, self.excludedVolumeCutOff)   
                self.polymer.addnewForce(repulsionForce)
    
                # Wait some more time
                self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)

            ##################################################                 
            # Simulation until NumSteps


            T0 = 0
            for it, T in enumerate(self.times2sample):

                self.polymer.step(T-T0,self.dt,self.diffusionConstant)

                # Distance breaks to center
                breaksPosition = self.polymer.get_r()[self.polymer.freeMonomers]
                breaks2center[i][it] = np.linalg.norm(breaksPosition - self.polymer.getCenterOfMass(), axis = 1)
                    
                msrg[i][it] = self.polymer.get_msrg()

                T0 = T
                       

            self.polymer = self.polymer.new()
            
            transcurredTime = time() - start
            totalTime = transcurredTime * self.numRealisations/(i+1)
            print('\r' + 'Simulation ' + str(i) + ' of ' + str(self.numRealisations) + ' | ' +
                  str(round(totalTime - transcurredTime, 2)) + ' sec still.', end='\r')
        
        mean_msrg = np.mean(msrg, axis = 0)
        mean_breaks2center = np.mean(breaks2center, axis = 0)
        
        msrg_dx = 1.96*np.std(msrg, axis = 0)/np.sqrt(self.numRealisations)
        mean_breaks2center_dx = 1.96*np.std(breaks2center, axis = 0)/np.sqrt(self.numRealisations)
        
        self.addResults("realtime",np.arange(self.times2sample[-1])*self.dt)
        for it, T in enumerate(self.times2sample):
            real = str(round(T*self.dt,1))
            self.addResults("MSRG_at_"+real,mean_msrg[it])
            self.addResults("MSRGdx_at_"+real,msrg_dx[it])
            for b in range(4):
                self.addResults("Distance_"+str(b)+"_2center_at_"+real,mean_breaks2center[it][b])
                self.addResults("Distance_"+str(b)+"_2center_at_"+real+"_dx",mean_breaks2center_dx[it][b])
                    
        

    def watchEncounter(self):

        # Simulates the break and some waiting time:
        removedNum = self.randomBreaks_SingleStep()

        # ADD EXCLUDED VOLUME
        if self.excludedVolumeCutOff > 0:
            try:
                kappa = self.excludedVolumeSpringConstant
            except:
                kappa = 3*self.diffusionConstant/(self.polymer.b**2)    
            repulsionForce = lambda polymer : - kappa * RepairSphere(polymer, self.polymer.freeMonomers, self.excludedVolumeCutOff)   
            self.polymer.addnewForce(repulsionForce)

        self.polymer.colors = ['y']*self.polymer.numMonomers
        for m in self.polymer.freeMonomers:
            self.polymer.colors[m] = 'r'
        self.polymer.colors[:self.polymer.freeMonomers[0]] = ['g']*(self.polymer.freeMonomers[0])
        self.polymer.colors[(self.polymer.freeMonomers[1]+1):self.polymer.freeMonomers[2]] = ['b']*(self.polymer.freeMonomers[2]-self.polymer.freeMonomers[1]-1)
#        self.polymer.colors[(self.polymer.freeMonomers[2]+1):self.polymer.freeMonomers[3]] = ['b']*(self.polymer.numMonomers-self.polymer.freeMonomers[3]-1)
        # Wait some more time
        self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)

        ##################################################                 
        # Simulation until encounter              
        t = 0
        trajectory = np.zeros((self.numSteps+1,self.polymer.numMonomers,3))
        trajectory[0] = self.polymer.get_r()
        
#        didEncounter = self.polymer.anyEncountered(self.encounterDistance)
        while(t < self.numSteps):
            self.polymer.step(1,self.dt,self.diffusionConstant)
#            didEncounter = self.polymer.anyEncountered(self.encounterDistance)
            trajectory[t+1] = self.polymer.get_r()
            t += 1
        
        if(t < self.numSteps):
            trajectory[t+1:-1] = self.polymer.get_r()
        
        self.trajectoire = trajectory
        
        msd = np.zeros(t)
        msrg = np.zeros(t)
        freeEnds_msd = np.zeros((len(self.polymer.freeMonomers), t))
        for it in range(t):
            msd[it] = np.mean(np.linalg.norm( trajectory[it] - trajectory[0] , axis = 1)**2)
            msrg[it] = np.mean( np.linalg.norm(trajectory[it] - np.mean(trajectory[it], axis = 0), axis = 1)**2)
            for im, m in enumerate(self.polymer.freeMonomers):
                freeEnds_msd[im][it] = (np.linalg.norm( trajectory[it][m] - trajectory[0][m])**2)

        
        
        self.addResults("realtime",np.arange(t)*self.dt)
        self.addResults("polymerMSD",msd)
        self.addResults("MSRG",msrg)
        for im, m in enumerate(self.polymer.freeMonomers):
            self.addResults(self.polymer.freeMonomersNames[m]+"MSD", freeEnds_msd[im])
        self.addResults("trajectory",self.trajectoire)
        self.addResults("removedCLs",removedNum)
        self.addResults("FET",t)
        
        ##################################################  

    def saveEncounterResults(self, FETs, events, removedNums, post_msrgs = None):
        # Prepare results 
        FETs = FETs*self.dt
        self.addResults("events",events)
        eventsC = Counter(events)
        total_valid_experiments = sum(eventsC.values())-eventsC['NA']
        
        if total_valid_experiments == 0:
            repairProbas = 0.
            repairHalfCI = 0.
            mFET = np.nan
            repairMFET = np.nan
            misrepairMFET = np.nan
            halfCI_mFET = 0.
            meanremovedCLs = 0.
            fitAmplitude = fitRate = 0.
            print('No valid experiments!')
        else:
            proba = eventsC['Repair']/total_valid_experiments
            repairProbas = proba
            repairHalfCI = 1.96*np.sqrt((proba - proba**2)/total_valid_experiments)
            mFET = np.nanmean(FETs)
            repairMFET = np.nanmean(FETs[events == 'Repair'])
            misrepairMFET = np.nanmean(FETs[(events != 'Repair') * (events != 'NA')])
            halfCI_mFET = 1.96*np.nanstd(FETs)/np.sqrt(total_valid_experiments)
            meanremovedCLs = np.nanmean(removedNums)
            fitAmplitude, fitRate = self.exponentialFit(FETs[~np.isnan(FETs)])

        
        # Save results
        self.addResults("FETs",FETs)
        self.addResults("meanFET",mFET)
        self.addResults("repairMFET",repairMFET)
        self.addResults("misrepairMFET",misrepairMFET)
        self.addResults("halfCI_FET",halfCI_mFET)
        self.addResults("expFit_Amplitude",fitAmplitude)
        self.addResults("expFit_Rate",fitRate)
        self.addResults("eventsCounter",eventsC)
        self.addResults("repair_probability",repairProbas)
        self.addResults("repair_CIhalflength",repairHalfCI)
        self.addResults("mean_removedCLs",meanremovedCLs)
        if post_msrgs is not None:
            self.addResults('MSRG_atEncounter', post_msrgs)
            self.addResults('Ensemble_MSRG', np.nanmean(post_msrgs))
            self.addResults('Ensemble_MSRG_dx', 1.96*np.nanstd(post_msrgs)/np.sqrt(total_valid_experiments))
