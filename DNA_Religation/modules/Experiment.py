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
from .Forces import ExcludedVolume, LocalExcludedVolume
from .Polymers import RCLPolymer
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
        
        assert np.sqrt(2*self.diffusionConstant*self.dt) <= 0.2*self.encounterDistance, "dt should verify: sqrt(2*D*dt) < 0.2*epsilon"
        
        if customExperiment == "SimpleDynamicSimulation":
            print("Simulation of a Simple Dynamic (no breaks)")
            self.runSimpleDynamic()
        elif customExperiment == "EncounterSimulation":
            print("Simulation of Encounter after two DSBs")
            self.runEncounterSimulation()
        elif customExperiment == "twoDSB":
            print("Simultation of two DSBs until relaxation time")
            self.runTwoRandomBreakSimulation()
        elif customExperiment == "Encounter_withExcludedVolume":
            print("Simulation of Encounter after two DSBs adding exclusion forces")
            self.BreakAndExclude()
        elif customExperiment == "TAD_Repair":
            print("Simulation of TAD Repair with one DSB in each TAD")
            self.TAD_repair()
        elif customExperiment == "oneTAD_Repair":
            print("Simulation of the repair of a TAD with two DSB")
            self.oneTADrepair()
        elif customExperiment == 'persistentDSB':
            print("Simulation of two persistent DSB motion")
            self.persistentDSB()
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

    def plot_trajectoire(self):
        fig = plt.figure()
        ax = Axes3D(fig)
#        ax.auto_scale_xyz([-10, 10], [-10, 10], [-10, 10])
        
        
        
        X = self.trajectoire[0][:,0]
        Y = self.trajectoire[0][:,1]
        Z = self.trajectoire[0][:,2]
        line, = ax.plot(X, Y, Z)
        dots = ax.scatter(X, Y, Z, c=self.polymer.colors, marker='o')

#        annotations = []
#        for m in self.monomers2follow:
#            rm = self.trajectoire[0][m,:]
#            annotations.append(ax.text(rm[0],rm[1],rm[2],  '%s' % (str(m)), size=5, zorder=1, color='k') )
        
#        fT = self.trajectoire[0][self.monomers2follow,:]
#        fX = fT[:,0]
#        fY = fT[:,1]
#        fZ = fT[:,2]
#        dots2follow = ax.scatter(fX, fY, fZ, c='r', marker='o')
        
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
            
#            ri2follow = ri[self.monomers2follow,:]
#            dots2follow._offsets3d = (ri2follow[:,0], ri2follow[:,1], ri2follow[:,2])
            
#            annotations[0].set_x(ri[m,0])
#            annotations[0].set_y(ri[m,1])
#            annotations[0].set_3d_properties(ri[m,3])
                
            return line, 
                
        def init():
            line.set_data([],[])
            line.set_3d_properties([])
            return line,
        
        ani = animation.FuncAnimation(fig, animate, frames=self.numSteps, init_func=init,
                                              interval=self.dt, blit=False, repeat = True)
    
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
   

    
    def randomBreaks_SingleStep(self):
        
        # Prepare the random DSBs
        breakLoci = self.polymer.randomCuts(self.genomicDistance,self.Nb)
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

        self.saveEncounterResults(FETs, events, removedNums)
        self.addResults('MSRG_atEncounter', post_msrgs)
        self.addResults('Ensemble_MSRG', np.nanmean(post_msrgs))
    
    
    
    def BreakAndExclude(self):
        """
        After inducing the DSBs add exclusion forces (excluded volume forces)
        from the cleaved monomers.
        """
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

        self.saveEncounterResults(FETs, events, removedNums)
        self.addResults('MSRG_atEncounter', post_msrgs)
        self.addResults('Ensemble_MSRG', np.nanmean(post_msrgs))


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
        self.addResults('MeanInterBreakDistance',np.mean(interbreakDistance,axis=0))


    def saveEncounterResults(self, FETs, events, removedNums):
        # Prepare results 
        FETs = FETs*self.dt
        self.addResults("events",events)
        events = Counter(events)
        total_valid_experiments = sum(events.values())-events['NA']
        
        if total_valid_experiments == 0:
            repairProbas = 0.
            repairHalfCI = 0.
            mFET = np.nan
            halfCI_mFET = 0.
            meanremovedCLs = 0.
            print('No valid experiments!')
        else:
            proba = events['Repair']/total_valid_experiments
            repairProbas = proba
            repairHalfCI = 1.96*np.sqrt((proba - proba**2)/total_valid_experiments)
            mFET = np.nanmean(FETs)
            halfCI_mFET = 1.96*np.nanstd(FETs)/np.sqrt(total_valid_experiments)
            meanremovedCLs = np.nanmean(removedNums)
        
        # Save results
        self.addResults("FETs",FETs)
        self.addResults("meanFET",mFET)
        self.addResults("halfCI_FET",halfCI_mFET)
        self.addResults("eventsCounter",events)
        self.addResults("repair_probability",repairProbas)
        self.addResults("repair_CIhalflength",repairHalfCI)
        self.addResults("mean_removedCLs",meanremovedCLs)

        