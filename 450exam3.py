# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
Exam III
"""

import matplotlib.pyplot as plt
from scipy import signal as sig
import math
import numpy as np
from math import log10

############################ Finding n for each possibility ###############

##### Finding n for Butterworth 
wc = 860
wp = 1/(1000/wc)
Hp = .95
ws = 1/(600/wc)
Hs = .07

na = log10( 1/(Hp**2) - 1)/(2*log10(wp) )
nb = log10( 1/(Hs**2) - 1)/(2*log10(ws) )
print('Butterworth:')
print('np = ',math.ceil(na))
print('ns = ',math.ceil(nb))

####### Finding n for Cheb I 
wp = 1
ws = 1/(600/1000)

eps = np.sqrt((1/Hp**2)-1)
n1 = math.acosh((1/eps)*np.sqrt((1/Hs**2)-1))
n2 = 1/math.acosh(ws)
n = n1*n2
n = math.ceil(n)
print('Chebyshev I:\n','n = ',n)

######## Finding n for Cheb II 
ws = 1
wp = 1/(1000/600)
print(wp,ws)

eps = np.sqrt(Hs**2/(1-Hs**2))
n1 = math.acosh((1/eps)*np.sqrt((1/Hs**2)-1))
n2 = 1/math.acosh(1/wp)
n = n1*n2
n = math.ceil(n)
print('Chebyshev II:\n','n = ',n,'\n\n')


################### Calculating num and den for Cheb I  ###########################
n = 5
wp = 1
ws = 1/(600/1000)

eps = np.sqrt((1/Hp**2)-1)
alpha = (1/eps)+np.sqrt(1+(1/(eps**2)))
a = .5*((alpha**(1/n))-(alpha**(-1/n)))
b = .5*((alpha**(1/n))+(alpha**(-1/n)))

Dtheta = 180/n
angle1 = Dtheta*np.pi/180
angle2 = 2*Dtheta*np.pi/180

d1 = [1, 2*a*np.cos(angle1), (a**2)*((np.cos(angle1))**2)+(b**2)*((np.sin(angle1))**2)]
d2 = [1, 2*a*np.cos(angle2), (a**2)*((np.cos(angle2))**2)+(b**2)*((np.sin(angle2))**2)]
d_odd = [1,a]

d12 = np.convolve(d1,d2)
den = np.convolve(d12,d_odd)

num = np.zeros(np.size(den))
K = den[n]
num[n] = K

print('num = ',num,'\n','den = ',den)

system1 = sig.lti(num,den)
w1, Hmag1, Hphase1 = sig.bode(system1)


########################### Trasnforming from Low to High Pass #############
num = num[::-1]
den = den[::-1]

print('num = ',num,'\n','den = ',den)

system2 = sig.lti(num,den)
w2, Hmag2, Hphase2 = sig.bode(system2)


########################### Scaling to original wp ########################
den = np.array(den)
w = 950
den[1] = (w**1)*den[1]
den[2] = (w**2)*den[2]
den[3] = (w**3)*den[3]
den[4] = (w**4)*den[4]
den[5] = (w**5)*den[5]

print('num = ',num,'\n','den = ',den)

system3 = sig.lti(num,den)
w3, Hmag3, Hphase3 = sig.bode(system3)



############################### Plotting Bode #############################

####### LP CHEB
plt.figure()
Ampl = 10**(0.05*Hmag1)
markers = [600,1000]
plt.semilogx(w1,Ampl,'k',markevery=markers)  # Plot amplitude, not dB.
plt.title('Base Chebyshev I Filter (n=5)')
plt.yticks([0,Hp,Hs,1])
plt.xticks([wp,ws])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.show()

########### HP CHEB
plt.figure()
Ampl = 10**(0.05*Hmag2)
markers = [600,1000]
plt.semilogx(w2,Ampl,'k',markevery=markers)  # Plot amplitude, not dB.
plt.title('High Pass Chebyshev I Filter (n=5)')
plt.yticks([0,Hp,Hs,1])
plt.xticks([wp,1/ws])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.show()


############ Scaled HP CHEB
plt.figure()
Ampl = 10**(0.05*Hmag3)
markers = [600,1000]
plt.semilogx(w3,Ampl,'k',markevery=markers)  # Plot amplitude, not dB.
plt.title('Scaled High Pass Chebyshev I Filter (n=5)')
plt.yticks([0,Hp,Hs,1])
plt.xticks([600,1000])
plt.vlines(1000,0,1)
plt.vlines(600,0,1)
#plt.hlines(.07,0,10000)
#plt.hlines(.95,0,10000)
#plt.axis([100,10000,0,1.1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.show()


############### The sinusodal input ##############################
dt = 0.0001
NN = 2000
DT = dt*1000 
TT = np.arange(0,NN*dt,dt)
y = np.zeros(NN)
f = np.zeros(NN)

A,B,C,D = sig.tf2ss(num,den)
x = np.zeros(np.shape(B))


omega = 600

for n in range(NN):
    f[n] = math.sin(omega*n*dt)

plt.figure()
plt.plot(TT,f,'k')
plt.grid()
plt.ylabel('f')
plt.show()

for m in range(NN):
    x = x + dt*A.dot(x) + dt*B*f[m]
    y[m] = C.dot(x) + D*f[m]

plt.figure() 
plt.plot(TT,y,'k')
plt.grid()
plt.xlabel('T (sec)')
plt.ylabel('y')
plt.show()





