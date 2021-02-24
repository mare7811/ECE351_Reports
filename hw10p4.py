# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
HW 10
Problem 8.4.2
"""
import matplotlib.pyplot as plt
from scipy import signal as sig
import math
import numpy as np


Hp = .95
Hs = .15
ws = 1
wp = .5
eps = np.sqrt(Hs**2/(1-Hs**2))

n1 = math.acosh(np.sqrt((1/eps**2)*(1/(1-Hp**2))))
n2 = 1/(math.acosh(1/wp))
n = n1*n2
n = math.ceil(n)
print('n = ',n)

alpha = (1/eps)+np.sqrt(1+1/eps**2)
a = .5*(alpha**(1/n)-alpha**(-1/n))
b = .5*(alpha**(1/n)+alpha**(-1/n))


Dtheta = 180/n
angle1 = (Dtheta/2)*np.pi/180
angle2 = ((Dtheta/2)+Dtheta)*np.pi/180
angle3 = ((Dtheta/2)+(2*Dtheta))*np.pi/180

s1 = a*np.cos(angle1) + 1j*b*np.sin(angle1)
q1 = 1/s1
q1c = np.conjugate(q1)
s2 = a*np.cos(angle2) + 1j*b*np.sin(angle2)
q2 = 1/s2
q2c = np.conjugate(q2)
s3 = a*np.cos(angle3) + 1j*b*np.sin(angle3)
q3 = 1/s3
q3c = np.conjugate(q3)

d1 = [1,np.real(q1+q1c),np.real(q1*q1c)]
d2 = [1,np.real(q2+q2c),np.real(q2*q2c)]
d3 = [1,np.real(q3+q3c),np.real(q3*q3c)]
d12 = np.convolve(d1,d2)
den = np.convolve(d12,d3)

num = np.zeros(n+1) 
omega1 = 1/np.cos(np.pi/8)
omega2 = 1/np.cos(3*np.pi/8)
omega3 = 1/np.cos(5*np.pi/8)

K = den[6]*(omega1**2)*(omega2**2)*(omega3**2)

n1 = [1,0,omega1**2]
n2 = [1,0,omega2**2]
n3 = [1,0,omega3**2]

n12 = np.convolve(n1,n2)
num = K*np.convolve(n12,n3)



print('num = ',num,'\n','den = ',den)
system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)



############# Plotting Bode #############################
Ampl = 10**(0.05*Hmag)
#plt.subplot(211)
plt.semilogx(w,Ampl,'k')  # Plot amplitude, not dB.
plt.title('Problem 8.4.2')
#plt.axis([10e3, 10e6, 0, 1])
#plt.yticks([0,Hp,Hs, 1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([.8,1])
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


omega = ws

for n in range(NN):
    f[n] = math.sin(omega*n*dt)

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
plt.text(10,.5,'omega = 1.52',fontsize=12) 
plt.grid()
plt.xlabel('T (sec)')
plt.ylabel('y')
plt.show()
    
    
