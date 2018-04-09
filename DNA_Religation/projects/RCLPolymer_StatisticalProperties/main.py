# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 10:40:11 2018

@author: ignacio
"""

import numpy as np
import matplotlib.pyplot as plt

from Simulation import Simulation
from Polymers import RCLPolymer


nb_monomers         = 30

numRealisations     = 50
numConectors        = 0 #number of added connectors

numSteps            = 200
dt                  = 0.01
dt_relax            = 0.05

dimension           = 3
diffusionCte        = 0.008
b                   = 0.2


p0 = RCLPolymer(nb_monomers, dimension, b, numConectors)

mc = Simulation(numSteps, dt, diffusionCte, p0, dt_relax, numRealisations)
mc.run()

msd = mc.get_msd_per_monomer()
mmsd = mc.get_avg_msd()
msrg = mc.get_msrg()
b2e = mc.get_b2e()

# Function to fit to the MSD
f = lambda t, A, alpha : A*t**alpha

# Fitted parameters
fit_params = mc.MSD_fit_mean()

print('Estimated MSRG : '+str(msrg))
print('Real MSRG      : '+str(nb_monomers*b*b/6))
print('Estimated b2   : '+str(b2e))
print('Real b2        : '+str(b*b))

plt.figure()
plt.xlabel('time (sec)')
plt.ylabel('MSD (mu^2)')
for mm in range(nb_monomers):
    plt.plot(np.arange(0,numSteps*dt,dt),msd[:,mm],label='monomer '+str(mm))
#    popt = stats.MSD_fit_monomer(mm)
#    plt.plot(mc.times, f(mc.times, *popt),'b--',
#         label='fit: A=%5.3f, alpha=%5.3f' % tuple(popt))

#train_times = np.repeat([mc.times],nb_monomers,axis=0)
#train_msd = msd.T
#popt, pcov = curve_fit(f, train_times.flatten(), msd.T.flatten())
#plt.plot(mc.times, f(mc.times, *popt), 'r--',
#         label='fit: A=%5.3f, alpha=%5.3f' % tuple(popt))

plt.plot(mc.timeline, f(mc.timeline, *fit_params), 'r--',
         label='fit: A=%5.3f, alpha=%5.3f' % tuple(fit_params))

plt.legend()
plt.show()

plt.figure()
plt.xlabel('time (sec)')
plt.ylabel('MSD (mu^2)')
plt.plot(np.arange(0,numSteps*dt,dt),mmsd,label='avg of all monomers')


plt.plot(mc.timeline, f(mc.timeline, *fit_params), 'r--',
         label='fit: A=%5.3f, alpha=%5.3f' % tuple(fit_params))
plt.legend()
plt.show()

#plt.figure()
#plt.xlabel('Monomer number')
#plt.ylabel('EP')
#monoref = 0
#monoms = np.arange(1,nb_monomers)
#plt.plot(monoms,ep1[1:])
#plt.show()
