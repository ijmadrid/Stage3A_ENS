# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:34:41 2018

@author: ignacio
"""

from collections import defaultdict
import numpy as np

class Graph():
    
    def __init__(self,vertices,**kwargs):
        self.V= vertices #No. of vertices
        self.adj = defaultdict(list) # adjacency list
        adjmatrix = kwargs.get('adjmatrix',None)
        if(adjmatrix is not None):
            np.fill_diagonal(adjmatrix,0)
            self.toList(adjmatrix)
  
    # function to add an edge to graph
    def addEdge(self,u,v):
        if type(u) == np.ndarray:
            for i, ui in enumerate(u): 
                    self.adj[ui].append(v[i])
        else:
            self.adj[u].append(v)
    
    # function to convert adjacency matrix to adjacency list
    def toList(self,adjmatrix):
        for i in range(self.V):
            for j in range(self.V):
                if adjmatrix[i,j]!=0:
                    self.addEdge(i,j)
    
    # Function that returns reverse (or transpose) of this graph
    def getTranspose(self):
        g = Graph(self.V)
        # Recur for all the vertices adjacent to this vertex
        for i in self.adj:
            for j in self.adj[i]:
                g.addEdge(j,i)
        return g
    
    def DFSUtil(self,v,visited):
        # Mark the current node as visited and print it
        visited[v]= True
        #Recur for all the vertices adjacent to this vertex
        for i in self.adj[v]:
            if visited[i]==False:
                self.DFSUtil(i,visited)
    
    
    def isConnected(self):       
        # Step 1: Mark all the vertices as not visited (For first DFS)\
        visited =[False]*(self.V)     
        # Step 2: Do DFS traversal starting from first vertex.  
        self.DFSUtil(0,visited)   
        # If DFS traversal doesn’t visit all vertices, then return false.  
        for i in range(self.V):
            if (not visited[i]):
                return False 
#        # Step 3: Create the transposed graph
#        g = self.getTranspose() 
#        # Step 4:  Mark all the vertices as not visited (For second DFS)\
#        visited =[False]*(self.V) 
#        # Step 5: Do DFS for reversed graph starting from the same first vertex.  
#        g.DFSUtil(0, visited);   
#        # If 2nd DFS traversal doesn’t visit all vertices, then return false.  
#        for i in range(self.V):
#            if (not visited[i]):
#                return False    
        return True
    
    
    def cutEdge(self,i,j):
        self.adj[i].remove(j)
        self.adj[j].remove(i)
    
    def cutAllEdgesWith(self,i):
        assert i != 0  and i != self.V , 'Trying to remove cross-links from the extremes'
        for j in self.adj[i][1:]:
            self.cutEdge(i,j)
    

