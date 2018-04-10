# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:48:32 2018

@author: ignacio
"""

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from .Polymers import *

class Simulation():
    
    def __init__(self,numSteps,dt,diffusionConst,polymer,dt_relax,numRealisations,**kwargs):
        """
        initialize with numSteps, diffusionCoeff
        dt, and initialize polymer class
        """
        self.numSteps = numSteps
        self.dt = dt
        self.polymer = polymer
        self.trajectoire = []
        self.dt_relax = dt_relax
        self.D = diffusionConst
        self.numRealisations = numRealisations

        self.msrg = np.zeros(numRealisations) 
#        self.sds = np.zeros((numRealisations,numSteps,polymer.numMonomers)) 
#        self.b2e = np.zeros(numRealisations)
        self.timeline = np.arange(0,numSteps*dt,dt)
        self.wasRun = False

    def relaxSteps(self):
        z = 2*self.polymer.Nc/((self.polymer.numMonomers-1)*(self.polymer.numMonomers-2))
        return int(np.ceil(self.polymer.b**2/(self.polymer.dim*self.D*(self.polymer.numMonomers*z + 4.0*(1-z)*np.sin(np.pi/(2.0*self.polymer.numMonomers))**2))/self.dt_relax))
    
    def add_step(self,r1):
        self.trajectoire.append(r1)

    def run(self):
        for i in range(self.numRealisations):
            # Burn in until relaxation time
            self.polymer.step(self.relaxSteps(),self.dt_relax,self.D)

            self.msrg[i] = self.polymer.get_msrg()
#            self.b2e[i] = self.b2()
            
            # Simulation
            self.add_step(self.polymer.get_r().copy())
            for t in range(self.numSteps):
                self.polymer.step(1,self.dt,self.D)
                self.add_step(self.polymer.get_r().copy())
        
                # Results collect
#                self.sds[i,t] = self.computeSD(t)
                    
            self.polymer = self.polymer.new()
            self.trajectoire = []
        
        self.wasRun = True
    

    def computeMSRG(self):
        """
        Mean Square Radius of Gyration at final time
        """
        r = self.polymer.get_r() # last position
        return np.mean(np.linalg.norm(r - np.mean(r, axis=0), axis = 1)**2)
    
    def computeSD(self,t):
        """
        Square Deviation of each monomer at time t
        """
        r0 = self.trajectoire[0]
        rt = self.trajectoire[t]
        return np.linalg.norm( rt - r0 , axis = 1)**2
        
    def b2(self):
        r = self.trajectoire[0]
        return np.mean(np.linalg.norm(r[1:-1]-r[0:-2], axis=1)**2)

    def get_avg_msd(self):
        return np.mean(np.mean(self.sds,axis=0),axis=1)
 
    def get_msd_per_monomer(self):
        return np.mean(self.sds,axis=0)
    
    def get_msrg(self):
        return np.mean(self.msrg)
    
    def get_b2e(self):
        return np.mean(self.b2e)


    def MSD_fit_mean(self):
        """
        Fits A and alpha to MSD(t) = f(t) = A*t**alpha
        """
        f = lambda t, A, alpha : A*t**alpha
        popt, pcov = curve_fit(f, self.timeline, self.get_avg_msd())
        return popt
    
    def MSD_fit_monomer(self,i):
        """
        Fits A and alpha to MSD(t) = f(t) = A*t**alpha
        """
        f = lambda t, A, alpha : A*t**alpha
        popt, pcov = curve_fit(f, self.timeline, self.sampler.get_msd_per_monomer()[:,i])
        return popt


#    def get_distVar_vector(self):
#        return np.mean(self.distVars,axis=0)
#    def distVar(self,m,n):
#        """
#        Variance of the distances between monomer m and monomer n
#        """
#        dists = self.distances[:,m,n]
#        return np.mean(dists**2)


    def plot_polymer_at(self,t):
        fig = plt.figure()
        x = self.trajectoire[t].get_r()[:,0]
        y = self.trajectoire[t].get_r()[:,1]
        if(self.polymer.dim == 3):
            ax = fig.gca(projection='3d')
            z = self.trajectoire[t].get_r()[:,2]
            ax.plot(x,y,z)
        else:
            plt.plot(x,y)
    
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
        

class EncounterSimulation(Simulation):
    
    def __init__(self,dt,diffusionConst,polymer,dt_relax,Nrealisations,maxIterationsPerExperiment,numBreaks,genomicDistance,encounterDistance,waitingSteps,simul_until=True):
        Simulation.__init__(self,maxIterationsPerExperiment,dt,diffusionConst,polymer,dt_relax,Nrealisations)
        self.genomicDistance   = genomicDistance
        self.encounterDistance = encounterDistance
        self.FETs = []
        self.events = []
        self.repairProba = []
        self.msrg = []
        self.Nb = numBreaks
        self.waitingSteps = waitingSteps
        self.simulate_until_encounter = simul_until
    
    def run(self):
        
        for i in range(self.numRealisations):
            
            # Prepare the random DSBs
            breakLoci = self.polymer.randomCuts(self.genomicDistance,self.Nb)
                        
            # Verify is polymer is splittable for the prepeared DSBs
            while(not(self.polymer.isSplittable(breakLoci))):
                # if not, make new connections (#TODO: try an heuristic maybe?)
                self.polymer.reset()
            
            # Once the polymer is splittable:
            # Burn in until relaxation time
            self.polymer.step(self.relaxSteps(),self.dt_relax,self.D)
            
            # Induce DSBs
            self.polymer.cutNow(breakLoci,definitive=True)
            # Remove the CL concerning the cleavages if it is the case
            if self.polymer.keepCL == False:
                self.polymer.removeCL()
            
            # Wait some more time
            self.polymer.step(self.waitingSteps,self.dt_relax,self.D)
            
            self.msrg.append(self.polymer.get_msrg())
            
            if self.simulate_until_encounter:
            # Simulation until encounter
                self.add_step(self.polymer.get_r().copy())
                t = 0
                while(not(self.polymer.anyEncountered(self.encounterDistance)[0]) and t < self.numSteps):
                    self.polymer.step(1,self.dt,self.D)
                    self.add_step(self.polymer.get_r().copy())
            
                    # Results collect
    #                self.sds[i,t] = self.computeSD(t)
                    t += 1
                
                if(t<self.numSteps):
                    self.FETs.append(t)
                    self.events.append(self.polymer.anyEncountered(self.encounterDistance)[1])
    
                    
    #            self.b2e[i] = self.b2()
        
            self.polymer = self.polymer.new()
            self.trajectoire = []
        
        self.FETs = np.array(self.FETs)*self.dt
        self.wasRun = True
        
        self.events = Counter(self.events)
        total_valid_experiments = sum(self.events.values())
        
        if total_valid_experiments == 0:
            self.repairProba = [0.5,0.5]
            print('No valid experiments!')
        else:
            proba = self.events['Repair']/total_valid_experiments
            self.repairProba = [proba,
                                1.96*np.sqrt((proba - proba**2)/total_valid_experiments)]
        