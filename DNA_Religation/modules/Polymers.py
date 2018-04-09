# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:36:41 2018

@author: ignacio
"""

from .Graph import Graph
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
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
        
#    def new(self):
#        return Polymer(self.numMonomers, self.dim, self.b)

    def get_r(self):
        return self.positions

    def plot(self):
        fig = plt.figure()
        x = self.get_r()[:,0]
        y = self.get_r()[:,1]
        if(self.dim == 3):
            ax = fig.gca(projection='3d')
            z = self.get_r()[:,2]
            ax.plot(x,y,z)
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
        
    def Linear_Rouse_Matrix(self):
        m = self.numMonomers
        M = np.diag(-1*np.ones(m-1),-1) + np.diag(2*np.ones(m),0) + np.diag(-1*np.ones(m-1),1)
        M[0,0] = 1
        M[m-1,m-1] = 1
        return M 
    
    def showMatrix(self):
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')    
        X = range(self.numMonomers)
        Y = range(self.numMonomers)
        X, Y = np.meshgrid(X, Y)
        Z = self.LaplacianMatrix.flatten('F')
        
        #Define colormap and get values for barcolors from it
        cmap = plt.cm.RdYlBu_r
        norm = mpl.colors.Normalize(vmin=Z.min(), vmax=Z.max())
        barcolors = plt.cm.ScalarMappable(norm, cmap)


        ax.bar3d(X.flatten('F'),Y.flatten('F'),np.zeros_like(X.flatten('F')),0.5,0.5,Z,color=barcolors.to_rgba(Z),alpha=0.4)
        plt.show()


    def makeRCL(self,Nc):
        return RCLPolymer(self.numMonomers,self.dim,self.b,Nc)
    
    def step(self,numsteps,dt,D):
        for j in range(numsteps):
            self.positions += -(self.dim*D/(self.b**2))*np.dot(self.LaplacianMatrix,self.get_r())*dt + np.random.randn(self.numMonomers,self.dim)*np.sqrt(2.0*D*dt)
    

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
                return (True,'Fail')            
            
        return (False,None)
        #return np.prod(self.haveEncountered(possibleEncounters[:,0],possibleEncounters[:,1],eps))

class RCLPolymer(Polymer):
    
    def __init__(self,numMonomers,dim,b,Nc,keepCL=False,**kwargs):
        Polymer.__init__(self,numMonomers,dim,b,**kwargs)
        self.Nc = Nc
        self.connect()
        self.keepCL = keepCL
        
        #TODO
        # Faire la matrice B et la construire sans toucher la matrice de Rouse
        # Puis faire Laplacian = B + Rouse

    def new(self):
        return RCLPolymer(self.numMonomers, self.dim, self.b, self.Nc)
    
    def reset(self):
        """
        Reset the connections, keeping the polymer in its current position
        """
        self.__init__(self.numMonomers, self.dim, self.b, self.Nc, position = self.positions)

    def copy(self):
        return copy.deepcopy(self)
    
    def connect(self):
        # Tirage des arets aleatoires
        possible_pairs = np.vstack(np.triu_indices(self.numMonomers,k=2)).T
        Nl = int((self.numMonomers-2)*(self.numMonomers-1)/2)
        selected = possible_pairs[np.random.choice(Nl,size=self.Nc,replace=False)].T
        # Mise a jour de la matrice laplacienne
        self.LaplacianMatrix[selected[0],selected[1]] = -1
        self.LaplacianMatrix[selected[1],selected[0]] = -1
        for i in selected[0] : self.LaplacianMatrix[i,i] += 1 
        for i in selected[1] : self.LaplacianMatrix[i,i] += 1 
        # Mise a jour de la liste d'adjacence
        self.addEdge(selected[0],selected[1])
        self.addEdge(selected[1],selected[0])
    
    
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
        A1 = np.random.randint(0,self.numMonomers-1-(Nb-1)*(g+1))
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
        for m in self.freeMonomers[::2]:
            # Remove all CL from and to the concerned monomers
            self.LaplacianMatrix[m][:m-1] = 0 
            self.LaplacianMatrix[m][m+2:] = 0
            self.LaplacianMatrix[:,m][:m-1] = 0 
            self.LaplacianMatrix[:,m][m+2:] = 0
            self.LaplacianMatrix[m,m] = 1 
            self.cutAllEdgesWith(m)

                
    def offDiagPairs(self):
        """
        Return the new connected (off-3diagonal) pairs 
        """
        return np.transpose(np.nonzero(np.triu(self.LaplacianMatrix,k=2)))

    
    def validDSB(self,g,Nb):
        """
        If the DSB splits the polymer, new DSB are tried until the cleaved 
        polymer is fully connected
        """
        # Make new random cuts
        self.randomCut(g,Nb)
        # Remove the CL concerning the cleavages
        self.removeCL()
        
        if(self.isConnected()):
            self.generatePossibleEncounters()
            return
        
        else:
#            token = 0
            while(not(self.isConnected())):
                # Reconnect the RCL polymer
                self.reset()
                # Make new random cuts
                self.randomCut(g,Nb)
                # Remove the CL concerning the cleavages
                self.removeCL()
#                token += 1
#            print(str(token)+' essais to make a valid DSB')
            self.generatePossibleEncounters()
            return
        
    
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