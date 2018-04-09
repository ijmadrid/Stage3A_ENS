# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 14:32:29 2018

@author: ignacio
"""

import numpy as np

N = 5

R = np.cumsum(np.random.randn(N,3),0)


## Handmade
def alpha(p,n): 
    if p == 0:
        return np.sqrt(1/N)
    else:
        return np.sqrt(2/N)*np.cos((n-0.5)*p*np.pi/N)




A = np.zeros(N,N)

for p in range(N):
    for n in range(N):
        A[p,n] = alpha(p,n)

U = A @ R

#R_reconstructed = np.linalg.inv(A) @ U

## FFT

Ufft = np.fft.rfft(R,norm='ortho')
