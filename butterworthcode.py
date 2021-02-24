# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 18:33:33 2020

@author: marem
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sig
from math import pi, exp, cos, sin, log, sqrt

""" Two poles """
num = [0, 0, 1]
den = [1, 1.41, 1]

system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)

plt.subplot(211)
plt.semilogx(w,10**(0.05*Hmag),'k')  # Plot amplitude, not dB.
plt.title('Butter_filter')
plt.axis([0.1, 10, 0, 1])
plt.yticks([0, 0.1, .5, .707, 1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([.1,.8,1,2,10])
plt.savefig('Butter.png',dpi=300)
plt.show()

dt = 0.001
NN = 50000
TT = np.arange(0,NN*dt,dt)
y = np.zeros(NN)
f = np.zeros(NN)

A,B,C,D = sig.tf2ss(num,den)
x = np.zeros(np.shape(B))

""" The sinusodal input """
omega = 1
for n in range(NN):
    f[n] = sin(omega*n*dt)

plt.subplot(211)
plt.plot(TT,f,'k')
plt.yticks([-1, 0, 1 ])
plt.axis([0, NN*dt,-1,1])
plt.grid()
plt.ylabel('f')

for m in range(NN):
    x = x + dt*A.dot(x) + dt*B*f[m]
    y[m] = C.dot(x) + D*f[m]

plt.subplot(212)    
plt.plot(TT,y,'k')
plt.axis([0, NN*dt,-1,1])
plt.yticks([-1, -.707, 0, .707, 1 ])
plt.text(10,.5,'omega = {}'.format(round(omega,1)),fontsize=12)
plt.grid()
plt.xlabel('T (sec)')
plt.ylabel('y')
plt.savefig('y_vs_f.png',dps=300)
plt.show()
    
    
	

