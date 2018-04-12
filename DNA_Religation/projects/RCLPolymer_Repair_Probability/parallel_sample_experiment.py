# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 11:50:21 2018

@author: ignacio
"""

import sys
import numpy as np
from pathos.multiprocessing import Pool
from time import time

def integrate(its):
    # I totally cheated and tweaked the number of chunks
    # to get the fastest result
    chunks = 10000
#    chunk_size = its / chunks

    np.random.seed()  # Each process needs a different seed!

    sum = 0.0
    for i in range(chunks):  # For each chunk...
        # ...do a vectorised Monte Carlo calculation
        u = np.random.uniform(size=int(its/chunks))
        sum += np.sum(np.exp(-u * u))  # Do the Monte Carlo

	# We did 'its' total iterations in this process, so
    # normalise the result and return it
    return sum / float(its)

if __name__ == '__main__':
    start = time()
    num_procs = int(sys.argv[1])

    iters = 1000000000
    its = int(iters / num_procs)  # Each process gets a share of the iterations

    pool = Pool(processes=num_procs)
    
    # Each process calls 'integrate' with 'its' iterations
    result = pool.map(integrate, [its] * num_procs)
    
    # pool.map returns a list of length 'num_procs', with
    # element 'i' being the return value of 'integrate' from
    # process 'i'

    # Renormalise by the number of processors
    print (sum(result) / float(num_procs))
    print ('Time : '+str(time()-start))