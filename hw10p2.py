# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
HW 10
Problem 8.3.5
"""
import matplotlib.pyplot as plt
from scipy import signal as sig
import math
import numpy as np

Hp = 10**(-.5/20)
Hs = 10**(-20/20)
ws = 5
eps = np.sqrt((1/Hp**2)-1)
print(Hp,eps,1/Hp**2)

n1 = math.acosh((1/eps)*np.sqrt((1/Hs**2)-1))
n2 = 1/math.acosh(ws)
n = n1*n2
n = math.ceil(n)
print('n = ',n)

alpha = (1/eps)+np.sqrt(1+1/eps**2)
a = .5*(alpha**(1/n)-alpha**(-1/n))
b = .5*(alpha**(1/n)+alpha**(-1/n))

Dtheta = 180/n
angle1 = (Dtheta/2)*np.pi/180
angle2 = ((Dtheta/2)+Dtheta)*np.pi/180
#print(angle1*180/np.pi,angle2*180/np.pi)


num = np.zeros(n+1) 

d1 = [1,2*a*np.cos(angle1),(a**2)*(np.cos(angle1)**2)+(b**2)*(np.sin(angle1)**2)]
d2 = [1,2*a*np.cos(angle2),(a**2)*(np.cos(angle2)**2)+(b**2)*(np.sin(angle2)**2)]
den = d1

if(n%2 == 1):    #check if odd
    K = 1
else:          #even
    K = (1/np.sqrt(1+eps**2))*den[n]

num[n] = K


########## Low to High Pass Cheb #############
num[0] = num[2]
num[2] = 0
den[0] = den[2]
den[2] = 1

########## Scaling by 1000 ####################
den[1] = den[1]*1000
den[2] = 10**6

print('num = ',num,'\n','den = ',d1)
system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)



############# Plotting Bode #############################
Ampl = 10**(0.05*Hmag)
#plt.subplot(211)
plt.semilogx(w,Ampl,'k')  # Plot amplitude, not dB.
plt.title('Problem 8.3.5')
#plt.axis([10e3, 10e6, 0, 1])
plt.yticks([0,Hp,Hs, 1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([200,1000])
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
    
    
