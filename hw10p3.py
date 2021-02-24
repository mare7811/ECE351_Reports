# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
HW 10
Problem 8.3.6
"""
import matplotlib.pyplot as plt
from scipy import signal as sig
import math
import numpy as np


num = np.zeros(5)
den = np.zeros(5)

num[2] = 1.43*(500**2)
num[4] = 0
den[0] = 1
den[1] = 1.43*500
den[2] = 2*(5000**2) + 1.52*(500**2)
den[3] = 1.43*500*(5000**2)
den[4] = (5000**2)**2

print('num = ',num,'\n','den = ',den)
system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)



############# Plotting Bode #############################
Ampl = 10**(0.05*Hmag)
#plt.subplot(211)
plt.semilogx(w,Ampl,'k')  # Plot amplitude, not dB.
plt.title('Problem 8.3.6')
plt.axis([10e2, 10e4, 0, 1])
#plt.yticks([0,Hp,Hs, 1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([5000])
plt.show()



############### The sinusodal input ##############################
dt = 0.01
NN = 500 
DT = dt*1000 
TT = np.arange(0,NN*dt,dt)
y = np.zeros(NN)
f = np.zeros(NN)

A,B,C,D = sig.tf2ss(num,den)
x = np.zeros(np.shape(B))


#omega = ws

#for n in range(NN):
#    f[n] = math.sin(omega*n*dt)

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
    
    
