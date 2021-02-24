# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:54:31 2020

@author: marem
"""


import matplotlib.pyplot as plt
import numpy as np
from math import exp, cos, sin, log
from math import pi, sqrt, exp, cos, sin
from cmath import exp

NN = 10000
dt = 0.01
TT = np.arange(0,dt*NN,dt)
time1 = np.zeros(NN)
time2 = np.zeros(NN)
time3 = np.zeros(NN)
time4 = np.zeros(NN)

A = np.matrix('0 1 0 0;-.4 -9.04 0 0;0 0 0 1;-1 -2 -2 -3')
B1 = np.matrix('0; 1; 0;0')
B2 = np.matrix('0; 0; 0; 0')
x = np.matrix('0;0;0;0')

f = np.ones(NN)
for n in range(0,1000):
    f[n] = 1
for n in range(1000,NN):
    f[n] = 0
    
g = np.zeros(NN)
for n in range(1000,NN):
    g[n] = .5



nsteps = NN
for i in range(nsteps):
    time1[i] = x[0]
    time2[i] = x[1]
    time3[i] = x[2] 
    time4[i] = x[3]
    
    x = x + dt*A*x + dt*B1*f[i] +dt*B2*g[i]    
    
plt.subplot(111)
plt.title('Problem 6.2.2 simulation')
plt.plot(TT,time1,'k',label='x1')
plt.plot(TT,time2,'g',label='x2')
plt.plot(TT,time3,'b',label='x3')
plt.plot(TT,time4,'r',label='x4')
#plt.plot(TT,f)
#plt.axis([0,15,-.15,.1])
plt.legend()
plt.grid()
plt.savefig('Example-6-2-1.png',dpi=300)