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

w1 = 2*pi*.3
w2 = 2*pi*.15

A = np.matrix('0 1 0;-.03 -.05 1;0 0 -.198')
B1 = np.matrix('0; 0; 1')
B2 = np.matrix('0; 0; 1')
x = np.matrix('0;0;0')

f = np.zeros(NN)
for n in range(0,2000):
    f[n] = sin(w1*TT[n])
for n in range(2000,NN):
    f[n] = 0

g = np.zeros(NN)
for n in range(2000,4000):
    g[n] = sin(w2*TT[n])
for n in range(4000,NN):
    g[n] = 0

x[0] = 0
x[1] = 0
x[2] = 0

nsteps = NN
for i in range(nsteps):
    time1[i] = x[0]
    time2[i] = x[1]
    time3[i] = x[2] 
    
    x = x + dt*A*x + dt*B1*f[i] +dt*B2*g[i]
    
    
plt.subplot()
plt.title('Problem 6.2.4 simulation')
plt.plot(TT,time1,'k',label='x1')
plt.plot(TT,time2,'g',label='x2')
plt.plot(TT,time3,'b',label='x3')
#plt.plot(TT,g)
#plt.axis([0,15,-.15,.1])
plt.legend()
plt.grid()
plt.savefig('Example-6-2-4.png',dpi=300)
