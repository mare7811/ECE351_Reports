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

#H(s) calculator for a Butterworth filter
n = 2
Dtheta = 180/n
if(n%2 == 1):
    angle1 = Dtheta*np.pi/180
    angle2 = (2*Dtheta)*np.pi/180 
else:
    angle1 = (Dtheta/2)*np.pi/180
    angle2 = ((Dtheta/2)+Dtheta)*np.pi/180
H12 = [1,2*cos(angle1),1]
#H34 = [1,2*cos(angle2),1] #must be excluded if n < 4
H_odd = [1,1]

def poly_mult(A,B):
    prod = np.zeros(np.size(A)+np.size(B)-1)   
    
    for i in range(0,np.size(A)):
        for j in range(0,np.size(B)):
            prod[i+j] += A[i] * B[j]
    return prod

    
################ Setting num and den #########################
#den = poly_mult(H12,H34)    #must be excluded if n < 4
den = H12
if((n%2) == 1):
    den = poly_mult(den,H_odd) 

num = np.zeros(np.size(den))
num[np.size(den)-1] = 1

num[0] = 1
num[2] = 0

#Scaling
den[1] = 50e3*den[1]
den[2] = (50e3)**2*den[2]

print(num,'\n',den)
system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)


############# Plotting Bode #############################
Ampl = 10**(0.05*Hmag)
#plt.subplot(211)
plt.semilogx(w,Ampl,'k')  # Plot amplitude, not dB.
plt.title('Butter_filter')
#plt.axis([10e1, 10e4, 0, 1])
plt.yticks([0,.707, 1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([50000])
plt.show()

########### Finding frequency at given H ###################
H_index = 0
H = .707
for i in range(99):
    if(Ampl[i] > H and Ampl[i+1] < H):
        H_index = i

print('|H| = ',Ampl[H_index],'dB','\n\u03c9 = ',w[H_index],'rad/s')


############### The sinusodal input ##############################
dt = 0.01
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
    
    
