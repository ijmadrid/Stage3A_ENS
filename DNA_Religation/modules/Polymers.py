# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:36:41 2018

@author: ignacio
"""

from .Graph import Graph
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
from scipy.spatial.distance import pdist, squareform
from itertools import combinations, product
import copy


class Polymer(Graph):
    
    def __init__(self, numMonomers, dim, b, **kwargs):
        """
        initialize with numMonomers and b, the std of distance between adjacent monomers
        return initialized polymer class, i.e. initial position
        """
        self.numMonomers = numMonomers
        self.dim = dim
        position = kwargs.get('position',None)
        if(position is not None):
            self.positions = position
        else:
            self.positions = b*np.cumsum(np.random.randn(numMonomers,dim),0)
        self.b = b
        self.LaplacianMatrix = self.Linear_Rouse_Matrix()
        Graph.__init__(self,self.numMonomers,adjmatrix=self.LaplacianMatrix.copy())
        self.Nc = 0
        self.freeMonomers = []
        self.possibleEncounters = []
        self.freeMonomersNames = {}
        
        self.forces = []
        
        self.colors = ['r']*numMonomers
        
    def get_r(self):
        return self.positions

    def plot(self):
        fig = plt.figure()
        x = self.get_r()[:,0]
        y = self.get_r()[:,1]
        if(self.dim == 3):
            ax = fig.gca(projection='3d')
            z = self.get_r()[:,2]
            
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0

            mid_x = (x.max()+x.min()) * 0.5
            mid_y = (y.max()+y.min()) * 0.5
            mid_z = (z.max()+z.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)

            ax.plot(x,y,z)
            ax.scatter(x, y, z, c=self.colors, marker='o')
            if self.Nc > 0:
                clpairs = self.offDiagPairs()
                for clpair in clpairs:
                    mx = x[clpair[0]]
                    nx = x[clpair[1]]
                    my = y[clpair[0]]
                    ny = y[clpair[1]]
                    mz = z[clpair[0]]
                    nz = z[clpair[1]]
                    ax.plot([mx,nx],[my,ny],[mz,nz],color = 'r')
        else:
            plt.plot(x,y)

    def distance(self,m,n):
        """
        Return the euclidean distance between the monomers m and n
        """
        if type(m) == np.ndarray:
            return np.linalg.norm( self.positions[m] - self.positions[n] , axis = 1)
        
        return np.linalg.norm(self.positions[m]-self.positions[n])
    
    def distanceMatrix(self):
        return squareform(pdist(self.positions,'euclidean'))
    
    def interBreakDistance(self):
        return self.distance(self.possibleEncounters[:,0],self.possibleEncounters[:,1])        
    
    def Linear_Rouse_Matrix(self):
        m = self.numMonomers
        M = np.diag(-1*np.ones(m-1),-1) + np.diag(2*np.ones(m),0) + np.diag(-1*np.ones(m-1),1)
        M[0,0] = 1
        M[m-1,m-1] = 1
        return M 
    
    def showMatrix(self):
#        plt.figure()
#        ax = fig.add_subplot(111,projection='3d')    
#        X = range(self.numMonomers)
#        Y = range(self.numMonomers)
#        X, Y = np.meshgrid(X, Y)
#        Z = self.LaplacianMatrix.flatten('F')      
#        #Define colormap and get values for barcolors from it
#        cmap = plt.cm.RdYlBu_r
#        norm = mpl.colors.Normalize(vmin=Z.min(), vmax=Z.max())
#        barcolors = plt.cm.ScalarMappable(norm, cmap)
#        ax.bar3d(X.flatten('F'),Y.flatten('F'),np.zeros_like(X.flatten('F')),0.5,0.5,Z,color=barcolors.to_rgba(Z),alpha=0.4)
        laplacian = plt.matshow(self.LaplacianMatrix, cmap='Spectral')
        plt.colorbar(laplacian)
        plt.show()

    
    def relaxTime(self, D):
        z = 2*self.Nc/((self.numMonomers-1)*(self.numMonomers-2))
        return self.b**2/(2*D*(self.numMonomers*z)) # + 4.0*(1-z)*np.sin(np.pi/(2.0*self.numMonomers))**2))
        #changed to fit PRE paper
        
    def step(self,numsteps,dt,D):
        
#        springConstant = (self.dim*D/(self.b**2))
        
        if len(self.forces) > 0:
            # Other potential gradients have been added to the dynamics
            for j in range(numsteps):
#                print([f(self)])
                self.positions += - (np.sum([f(self) for f in self.forces],axis=0) + (self.dim*D/(self.b**2))*(np.dot(self.LaplacianMatrix,self.get_r())))*dt + np.random.randn(self.numMonomers,self.dim)*np.sqrt(2.0*D*dt)
        
        else:
            # Rouse Harmonic potential only
            for j in range(numsteps):
                self.positions += -(self.dim*D/(self.b**2))*np.dot(self.LaplacianMatrix,self.get_r())*dt + np.random.randn(self.numMonomers,self.dim)*np.sqrt(2.0*D*dt)
        
    def addnewForce(self, potentialGradient):
        self.forces.append(potentialGradient)
    
    def getCenterOfMass(self):
        return np.mean(self.positions, axis = 0)
    
    def haveEncountered(self,mono1,mono2,eps):
        """
        Return True if mono1 and mono2 have encountered, i.e., |R_mono1 - R_mono2| < eps
        """        
        return self.distance(mono1,mono2) < eps
    
    def anyEncountered(self,eps):
        """
        Return the tuple (Encounter, Event) where :
        Encounter is True if any pair of free monomers that have encountered
        If Encounter is True, Event is a String:
            Event is Repair if the encountered monomers are A1 and A2 or B1 and B2
            Event is Fail if not
        """
        
        endsDistance = self.distance(self.possibleEncounters[:,0],self.possibleEncounters[:,1])
        
        if np.min(endsDistance) < eps:
            #pair = self.possibleEncounters[np.argmin(endsDistance)]
            selectedPair = np.random.choice(np.arange(0,len(endsDistance))[endsDistance < eps]) # coin toss
            pair = self.possibleEncounters[selectedPair]
            if(pair[1]-pair[0]==1 or pair[0]-pair[1]==1):
                return (True,'Repair')
            else:
                tag = self.freeMonomersNames[pair[0]]+"-"+self.freeMonomersNames[pair[1]]
#                tag = ""
                return (True,'Misrepair_'+tag)     
            
        return (False,None)
        #return np.prod(self.haveEncountered(possibleEncounters[:,0],possibleEncounters[:,1],eps))
        
    def get_msrg(self):
        return np.mean(np.linalg.norm(self.positions - np.mean(self.positions, axis=0), axis = 1)**2)
    
    
    def get_params(self):
        return dict(b                   = self.b,
                    monomersNumber      = self.numMonomers,
                    crosslinksNumber    = self.Nc)

    
    def simulateExperiment(self,experimentSet,results=dict()):
        """
        Simulate the experiment defined by the Experiment object experimentSet
        """
        results["PolymerParams"]    = self.get_params()
        results["ExperimentParams"] = experimentSet.get_params()
        experimentSet.run(self, results)
    
    
    def potentialEnergy(self):
        """
        Returns the polymer total elastic potencial energy
        """
        return 0.5*(pdist(self.positions)**2).sum()

class RCLPolymer(Polymer):
    
    def __init__(self,numMonomers,dim,b,Nc, keepCL=False, **kwargs):
        
        
        assert (np.size(Nc) == np.size(numMonomers)**2), "numMonomers %s and Nc %s must have according shape!" % (np.shape(numMonomers), np.shape(Nc))
                   
        self.keepCL = keepCL
 
        if type(numMonomers) == np.ndarray :
            Polymer.__init__(self,numMonomers.sum(),dim,b,**kwargs)
            subdomainBoundaries = numMonomers.cumsum()
            start = 0
            for i, end in enumerate(subdomainBoundaries):
                self.domainConnect(start, end, Nc[i,i])
                for j in range(i+1, len(numMonomers)):
                    start2 = subdomainBoundaries[j-1]
                    end2 = subdomainBoundaries[j]
                    self.interDomainConnect(start, end, start2, end2, Nc[i,j])
                start = end
            self.Nc = Nc.sum()
            
            self.subDomainNc = Nc
            self.subDomainnumMonomers = numMonomers
       
        else:
            Polymer.__init__(self,numMonomers,dim,b,**kwargs)
            self.Nc = Nc
            self.randomConnect()
            
            self.subDomainNc = Nc
            self.subDomainnumMonomers = numMonomers

        

    def new(self):
        return RCLPolymer(self.subDomainnumMonomers, self.dim, self.b, self.subDomainNc, keepCL = self.keepCL)
    
    def reset(self):
        """
        Reset the connections, keeping the polymer in its current position
        """
        self.__init__(self.subDomainnumMonomers, self.dim, self.b, self.subDomainNc, self.keepCL, position = self.positions)

    def copy(self):
        return copy.deepcopy(self)
    
    def connect(self, selected):
        """
        Connect the pairs indicated in selected
        """
        # Mise a jour de la matrice laplacienne
        self.LaplacianMatrix[selected[0],selected[1]] = -1
        self.LaplacianMatrix[selected[1],selected[0]] = -1
        for i in selected[0] : self.LaplacianMatrix[i,i] += 1 
        for i in selected[1] : self.LaplacianMatrix[i,i] += 1 
        # Mise a jour de la liste d'adjacence
        self.addEdge(selected[0],selected[1])
        self.addEdge(selected[1],selected[0])
        
    def randomConnect(self):
        """
        Induce self.Nc random connections in the whole polymer
        """
        if self.Nc == 0:
            return
        else:
            possible_pairs = np.vstack(np.triu_indices(self.numMonomers,k=2)).T
            Nl = len(possible_pairs)
            selected = possible_pairs[np.random.choice(Nl,size=self.Nc,replace=False)].T
            self.connect(selected)
    
    def domainConnect(self,left,right,Nc):
        """
        Induce Nc random connections between the monomers in [left,right[
        """
        if Nc == 0:
            return
        else:
            domainLength = right - left
            # Tirage des arets aleatoires
            possible_pairs = np.vstack(np.triu_indices(domainLength,k=2)).T
            Nl = len(possible_pairs) #Nl = int((domainLength-2)*(domainLength-1)/2)
            selected = left + possible_pairs[np.random.choice(Nl,size=Nc,replace=False)].T
            self.connect(selected)
            # Color the domain
            self.colors[left:right] = ['g']*(right-left)
        
    
    def interDomainConnect(self,l1,r1,l2,r2,Nc):
        """
        Induce Nc random connecions between [l1, r1[ and [l2, r2[
        """
        if Nc == 0:
            return
        else:
            x = np.arange(l1,r1)
            y = np.arange(l2,r2)
            possible_pairs = [[a,b] for a,b in product(x,y)]
            if r1 == l2:
                possible_pairs.pop((r1-l1-1)*(r2-l2))
            Nl = len(possible_pairs)
            selected = np.array(possible_pairs)[np.random.choice(Nl,size=Nc,replace=False)].T
            self.connect(selected)
    
    def extraDomainConnect(self,Nc):
        """
        Induce connectors between monomers which are not from the same domain
        """
        if Nc == 0:
            return
        else:
            possible_pairs = np.vstack(np.triu_indices(self.numMonomers,k=2)).T.tolist()
            for l,r in self.TADs:
                length = r-l
                for i in np.arange(l,r):    
                    for j in range(length-i):
                        possible_pairs.remove([i,i+2+j])
            
            possible_pairs = np.array(possible_pairs)
            
            selected = possible_pairs[np.random.choice(len(possible_pairs),size=Nc,replace=False)].T
            
            self.connect(selected)


    def imposeDSBconnections(self,Nc,breakLoci):
        """
        Imposes connecitivy in the Damage Foci (monomers in the DSBs)
        NB : breakLoci = [A1, B1]
        """
        #TODO
        #Adapt to any Nb
        if Nc == 0:
            return
        else:
            breakLoci = np.array(breakLoci)
            breakLoci = np.vstack((breakLoci,breakLoci+1)).T.flatten() # [A1, B1] -> [A1, A2, B1, B2]
            possible_pairs = [[a,b] for a,b in product(breakLoci, np.arange(0,self.numMonomers))]

            k0 = 0
            k1 = 4
            
            if breakLoci[0] == 0:
                for _ in range(2):
                    possible_pairs.pop(0)
                k0 += 1
            
            if breakLoci[-1] == 99:
                for _ in range(2):
                    possible_pairs.pop(-1)
                k1 -= 1
                
            for i in range(k0,k1):
                for _ in range(3):
                    possible_pairs.pop(self.numMonomers*i + breakLoci[i] - 1 - 3*i + k0)
                
            for _ in range(2):
                possible_pairs.pop(breakLoci[2]-3 + k0)
            for _ in range(2):
                possible_pairs.pop(self.numMonomers + breakLoci[2] - 2*3 + k0 - 2)
            selected = np.array(possible_pairs)[np.random.choice(len(possible_pairs),size=Nc,replace=False)].T
            self.connect(selected)                
            self.Nc += Nc
        
    
    def randomCuts(self,g,Nb):
        """
        Define Nb random DSB in the polymer each at distance g
        0--...--o--o--o--A1 x A2--o--...--o--B1 x B2--o--o--o--...--N
                                <---- g ---->
                           
        g  : genomic distance
        
        Return the upstream monomer (the left one in the sense 0 --> N)
        of each DSB
        """
        # A1 ~ Unif[0,N-1-(Nc-1)(g-1)[
        A1 = np.random.randint(1, self.numMonomers-1-(Nb-1)*(g+1)-1)
        return A1 + np.arange(Nb)*(1+g)
    
    
    def inregionCut(self,l,r,g,Nb):
        """
        Define a random cut in the given region
        """
        A1 = np.random.randint(l+1-self.keepCL, r-1-(Nb-1)*(g+1)-(1-self.keepCL))
        return A1 + np.arange(Nb)*(1+g)
        
        
    def cutNow(self,leftMonomers,definitive=False):
        """
        Make the cuts defined by self.randomCuts(g,Nb)
        If defintive = True, save the DSB loci ans generates
        the possible pairs for the encounters (A1-B1, A1-B2, etc.)
        """
        # A1 ~ Unif[0,N-1-(Nc-1)(g-1)[
        for A1 in leftMonomers:
            A2 = A1 + 1
            # Mise a jour de la matrice laplacienne
            self.LaplacianMatrix[A1,A2] = 0
            self.LaplacianMatrix[A2,A1] = 0
            self.LaplacianMatrix[A1,A1] -= 1 
            self.LaplacianMatrix[A2,A2] -= 1 
            # Mise a jour de la liste d'adjacence
            self.cutEdge(A1,A2)
            # Add new free ends to freeMonomers list
            self.freeMonomers.extend([A1,A2])
        
        if definitive:
            self.generatePossibleEncounters()
#            
            for i in range(len(self.freeMonomers)):
                self.freeMonomersNames[self.freeMonomers[i]] = chr(97 + i//2) + str(1 + i%2)


    def isSplittable(self,cutLoci):
        """
        Verify is the RCL Polymer can be splittable by the random cuts
        for which left monomers are in cutLoci
        """
        mirror = self.copy()
        # Make new random cuts
        mirror.cutNow(cutLoci)
        # Remove the CL concerning the cleavages
        if self.keepCL == False:
            mirror.removeCL()        
        return mirror.isConnected()
        
        
    def removeCL(self):
        """
        Remove all cross links between for all the monomers which have
        been cleaved
        """
        removedNum = 0
        for m in self.freeMonomers: #[::2]:
            # Remove all CL from and to the concerned monomers
            self.LaplacianMatrix[m][:m-1] = 0 
            self.LaplacianMatrix[m][m+2:] = 0
            self.LaplacianMatrix[:,m][:m-1] = 0 
            self.LaplacianMatrix[:,m][m+2:] = 0
            self.LaplacianMatrix[m,m] = 1 
            removedNum += self.cutAllEdgesWith(m)
        # Update the diagonal of the Laplacian
        np.fill_diagonal(self.LaplacianMatrix, 0) # dummy, not sure if it is worth to think a better way
        np.fill_diagonal(self.LaplacianMatrix, -1*self.LaplacianMatrix.sum(axis = 1))
        
        return removedNum

                
    def offDiagPairs(self):
        """
        Return the new connected (off-3diagonal) pairs 
        """
        return np.transpose(np.nonzero(np.triu(self.LaplacianMatrix,k=2)))

        
    
    def generatePossibleEncounters(self):
        """
        Once the polymer is cut, this method returns the possible combinations
        of encounters
        """
        possibleEncounters = [i for i in combinations(self.freeMonomers,r=2)]
        self.possibleEncounters = np.array(possibleEncounters)
        

        

class BetaPolymer(Polymer):
    
    def __init__(self,numMonomers,dim,b,beta):
        Polymer.__init__(self,numMonomers,dim,b)
        self.beta = beta
        self.modes = self.position_to_Fourier(self.positions)
        self.springConstant = 1.
        self.alphas = self.alphaMatrix()
        
        A = np.zeros((self.numMonomers,self.numMonomers))
    
        for p in range(self.numMonomers):
            for n in range(self.numMonomers):
                A[p,n] = self.alpha(p,n+1)
        
        self.alphas = A
    
    def make(self):
        #TODO
        raise NotImplementedError

    def alpha(self,p,n): 
        if p == 0:
            return np.sqrt(1/self.numMonomers)
        else:
            return np.sqrt(2/self.numMonomers)*np.cos((n-0.5)*p*np.pi/self.numMonomers)
     
    def eigenvalues(self):
        return  self.springConstant*np.array([4*np.sin(p*np.pi/2*self.numMonomers)**self.bet for p in range(self.numMonomers)])
        
    def spectral_step(self,numsteps,dt,D):
        Dp = np.ones(self.numMonomers)*D
        Dp[0] /= self.numMonomers
        eigenvalues = self.eigenvalues()
        for j in range(numsteps):
            self.modes += -Dp*eigenvalues*self.modes*dt + \
            np.random.randn(self.numMonomers,self.dim)*np.sqrt(2.0*Dp*dt)
    
    def position_to_Fourier(self):
        """
        Return the normal coordiantes (Fourier transform) for a given R
        NB : The shape of R is numMonomers x 3 (dim)
        """
        #TODO Try to do it with FFT 
        U = self.alphas @ self.positions
        
        return U
    
    def getPositions_from_Modes(self):
        
        return np.linalg.inv(self.alphas) @ self.modes


#
#class DomainRCLPolymer(Polymer):
#    """
#    aa
#    """
#    