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
import multiprocessing as mp

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
        # Unpackage params values
        for paramKeys, paramValues in params.items(): 
            exec('self.'+paramKeys+'=paramValues')
            
            
        if customExperiment == "SimpleDynamicSimulation":
            print("Simple Dynamic Simulation")
            self.runSimpleDynamic()
        elif customExperiment == "EncounterSimulation":
            print("After two DSBs Encounter Simulation")
            self.runEncounterSimulation()
        elif customExperiment == "twoDSB":
            print("Two DSBs until relaxation time (no encounter!)")
            self.runTwoRandomBreakSimulation()
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
        
        x = self.trajectoire[0][:,0]
        line, = ax.plot(x, self.trajectoire[0][:,1], self.trajectoire[0][:,2])
                
        def animate(i):
            ri = self.trajectoire[i]
            line.set_xdata(ri[:,0])  # update the data
            line.set_ydata(ri[:,1])  # update the data
            line.set_3d_properties(ri[:,2])
            return line,
                
        def init():
            line.set_data([],[])
            line.set_3d_properties([])
            return line,
        
        ani = animation.FuncAnimation(fig, animate, frames=self.numSteps, init_func=init,
                                              interval=self.dt, blit=False, repeat = True)
    
        return ani

    
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
   

    
    def randomBreak_SingleSimulation(self):
        
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
            self.polymer.removeCL()
        
        # Wait some more time
        self.polymer.step(self.waitingSteps,self.dt_relax,self.diffusionConstant)
    
    
    
    
    def runTwoRandomBreakSimulation(self):
        
        self.addResults("iterationsNumber",self.numRealisations)
        msrg = np.zeros(self.numRealisations)   
        
        for i in range(self.numRealisations): 
            self.randomBreak_SingleSimulation()
            msrg[i] = self.polymer.get_msrg()
            if(self.polymer.get_msrg() > 10):
                print("SHITY POLYMER")
                self.polymer.plot()
                print("Positions:")
                print(self.polymer.positions)
                print("Params:")
                print(self.polymer.get_params())
                print("Laplacian matrix:")
                print(self.polymer.LaplacianMatrix)
                print("Relaxation time:")
                print(self.polymer.relaxTime(self.diffusionConstant))  
            self.polymer = self.polymer.new()
            
        self.addResults("MSRG", np.mean(msrg))
        halfCI = 1.96*np.std(msrg)/np.sqrt(len(msrg))
        self.addResults("MSRG_95CI", halfCI)
        
    
    
    def runEncounterSimulation(self):
            
        self.addResults("iterationsNumber",self.numRealisations)
        FETs = []
        events = []
        repairProba = []
        
        for i in range(self.numRealisations):
            
            # Simulates the break and some waiting time
            self.randomBreak_SingleSimulation()
            
            # Simulation until encounter              
            t = 0
            while(not(self.polymer.anyEncountered(self.encounterDistance)[0]) and t < self.numMaxSteps):
                self.polymer.step(1,self.dt,self.diffusionConstant)
                t += 1
                
            if( t <self.numMaxSteps):
                FETs.append(t)
                events.append(self.polymer.anyEncountered(self.encounterDistance)[1])
        
            self.polymer = self.polymer.new()

        
        # Prepare results 
        FETs = np.array(FETs)*self.dt
        events = Counter(events)
        total_valid_experiments = sum(events.values())
        
        if total_valid_experiments == 0:
            repairProba = [0.5,0.5]
            print('No valid experiments!')
        else:
            proba = events['Repair']/total_valid_experiments
            repairProba = [proba,
                           1.96*np.sqrt((proba - proba**2)/total_valid_experiments)]
        
        # Save results
        self.addResults("FETs",FETs)
        self.addResults("eventsCounter",events)
        self.addResults("repair_probability_CI",repairProba)
        