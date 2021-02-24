# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
Hw 4
Problem 6.3.1
"""


import matplotlib.pyplot as plt
import numpy as np
from math import exp, cos, sin, log
from math import pi, sqrt, exp, cos, sin
from cmath import exp

NN = 10000
dt = 0.05
TT = np.arange(0,dt*NN,dt)
time1 = np.zeros(NN)
time2 = np.zeros(NN)
time3 = np.zeros(NN)
time4 = np.zeros(NN)
wtime1 = np.zeros(NN)
wtime2 = np.zeros(NN)
wtime3 = np.zeros(NN)
wtime4 = np.zeros(NN)

A = np.matrix('0 1 0 0;-9.04 -.4 0 0;0 0 0 1;-1 -2 -2 -3')
B1 = np.matrix('0; 1; 0;0')
B2 = np.matrix('0; 0; 0; 1')
x = np.matrix('0;0;0;0')
w = np.matrix('0; 0; 0; 0')


f = np.zeros(NN)
for n in range(0,200):
    f[n] = 1

    
g = np.zeros(NN)
for n in range(2000,NN):
    g[n] = .5


nsteps = NN
for i in range(nsteps):
    time1[i] = x[0]
    time2[i] = x[1]
    time3[i] = x[2] 
    time4[i] = x[3]
    
    x = x + dt*A*x + dt*B1*f[i] + dt*B2*g[i]    
    
plt.subplot(111)
plt.title('Problem 6.3.1 - state space formulation simulation (dt=0.1)')
plt.plot(TT,time1,'k',label='x1')
plt.plot(TT,time2,'g',label='x2')
plt.plot(TT,time3,'b',label='x3')
plt.plot(TT,time4,'r',label='x4')
#plt.xlim(0,50)
plt.legend()
plt.grid()
plt.savefig('Example-6-3-1.png',dpi=300)


I4 = np.eye(4)
A2 = A*A
A3 = A*A2
A4 = A2*A2
A5 = A*A4
A6 = A4*A2
F = (I4 + A*dt + 0.5*A2*(dt**2) + (1/6)*A3*(dt**3) + (1/24)*A4*(dt**4) + (1/120)*A5*(dt**5) + (1/720)*A6*(dt**6))
Ainv = np.linalg.inv(A)
G1 = (F-I4)*Ainv*B1
G2 = (F-I4)*Ainv*B2

nsteps = NN
for i in range(nsteps):
    wtime1[i] = w[0]
    wtime2[i] = w[1]
    wtime3[i] = w[2]
    wtime4[i] = w[3]
    
    w = F*w + G1*f[i] + G2*g[i] 
    
plt.figure()    
plt.subplot()
plt.title('Problem 6.3.1 - state transition formulation simulation (dt= 0.01)')
plt.plot(TT,wtime1,'k',label='x1')
plt.plot(TT,wtime2,'g',label='x2')
plt.plot(TT,wtime3,'b',label='x3')
plt.plot(TT,wtime4,'r',label='x4')
plt.plot(TT,f)
#plt.xlim(0,50)
plt.legend()
plt.grid()
plt.savefig('Example-6-3-1.png',dpi=300)

