# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
HW 8
Problem 8.1.2
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sig
from math import pi, exp, cos, sin, log, sqrt

num = np.zeros(5)
den = np.zeros(5)

num[2] = 90000
den[0] = 1
den[1] = 423
den[2] = (2e6)+90000
den[3] = 423e6
den[4] = (10**6)**2

print(num,'\n',den)
system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)


############# Plotting Bode #############################
Ampl = 10**(0.05*Hmag)
#plt.subplot(211)
plt.semilogx(w,Ampl,'k')  # Plot amplitude, not dB.
plt.title('Butter_filter')
plt.axis([600, 1400, 0, 1])
#plt.yticks([0,.707, 1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([850,1000,1150])
plt.show()

########### Finding frequency at given H ###################
H_index = 0
H = .707
for i in range(99):
    if(Ampl[i] > H and Ampl[i+1] < H):
        H_index = i

print('|H| = ',Ampl[H_index],'dB','\n\u03c9 = ',w[H_index],'rad/s')


############### The sinusodal input ##############################
dt = 0.0001
NN = 5000 

TT = np.arange(0,NN*dt,dt)
y = np.zeros(NN)
f = np.zeros(NN)

A,B,C,D = sig.tf2ss(num,den)
x = np.zeros(np.shape(B))


omega = w[H_index]

for n in range(NN):
    f[n] = sin(omega*n*dt)

plt.subplot(211)
plt.plot(TT,f,'k')
#plt.yticks([-1, 0, 1 ])
#plt.axis([0, NN*dt,-1,1])
plt.grid()
plt.ylabel('f')

for m in range(NN):
    x = x + dt*A.dot(x) + dt*B*f[m]
    y[m] = C.dot(x) + D*f[m]

plt.subplot(212)    
plt.plot(TT,y,'k')
#plt.axis([0, NN*dt,-1,1])
#plt.yticks([-1, -.707, 0, .707, 1 ])
#plt.text(10,.5,'omega = 1.52',fontsize=12) 
plt.grid()
plt.xlabel('T (sec)')
plt.ylabel('y')
plt.show()
    
    
