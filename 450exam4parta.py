# -*- coding: utf-8 -*-
"""
Miguel Mares
Ece 450 
Exam IV
Part A
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sig
from math import pi, exp, cos, sin, log, sqrt
from control import margin
from control import tf
from cmath import exp
import numpy as np 
from math import log10

########### Calculating n ##########################
wp = .5
Hp2 = .707

n = log10( 1/Hp2 - 1)/(2*log10(wp) )
print('np = ',n)



ws = 2
Hs2 = .0178

ns = log10( 1/Hs2 - 1)/(2*log10(ws) )
print('ns = ',ns)

####################Plotting Butterworth Filter#########################

num = [1]
d1 = [1,1]
d2 = [1,1,1]
den = np.convolve(d1,d2)

system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)

Ampl = 10**(0.05*Hmag)
#plt.subplot(211)
plt.semilogx(w,Ampl,'k')  # Plot amplitude, not dB.
plt.title('3rd Oder Butterworth')
plt.yticks([0,.707,.0178,1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([.5,1,2])
plt.show()

#plt.subplot(212)
plt.semilogx(w,Hphase,'k')  # Plot amplitude, not dB.
#plt.title('Angle Bode')
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('/_ |H|')
plt.yticks([0,-90,-180,-270])
#plt.xticks([100000,200000,50000])
plt.show()

num = [1e15]
den = [1,2*1e5,2*1e10,1e15]

system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system)

Ampl = 10**(0.05*Hmag)
plt.semilogx(w,Ampl,'k')  # Plot amplitude, not dB.
plt.title('Scaled 3rd Order Butterworth')
plt.yticks([0,.707,.0178, 1])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('|H|')
plt.xticks([100000,200000,50000])
plt.show()





##############Calculating and Plotting Z Bode############################

NN = 5000
phi = np.linspace(0,2*pi,NN)
z = np.zeros(NN, dtype = np.complex)
H = np.zeros(NN, dtype = np.complex)

T=.2
A = np.deg2rad(sqrt(.75)*T)
b = 1   
H1 = ((b/(b-exp(-T))) + (((.577)*exp(-.5*T)*sin(A)*b) - (b**2) + (exp(-.5*T)*cos(A)*b)) / ((b**2) - (2*exp(-.5*T)*cos(A)*b) + exp(-T)))
print("H1 = ",H1)
for n in range(0,NN):
    z = exp(1j*phi[n])
    H[n] = ((z/(z-exp(-T))) + (((.577)*exp(-.5*T)*sin(A)*z) - (z**2) + (exp(-.5*T)*cos(A)*z)) / ((z**2) - (2*exp(-.5*T)*cos(A)*z) + exp(-T)))/H1
    

phi1 = np.rad2deg(wp*T)
phi2 = np.rad2deg(ws*T)
print('phi_pass = ',phi1,'\n','phi_stop = ',phi2)
    
plt.subplot(211)
#plt.plot((180/pi)*phi,abs(H),'k')
plt.semilogx((180/pi)*phi,20*np.log10(H),'k')
plt.axis([1,100, -60, 10])
plt.ylabel('|G| dB')
plt.yticks([ -35,-3,0])
plt.axvline(phi1,color='k')
plt.axvline(phi2,color='k')
plt.text(2,-10,'$\phi1$ = {}'.format(round(phi1,2)),fontsize=12)
plt.text(8,-30,'$\phi2$ = {}'.format(round(phi2,2)),fontsize=12)
plt.title('zbode')
plt.grid(which='both')

aaa = np.angle(H)
#for n in range(NN):
#    if aaa[n] > 0:
#        aaa[n] = aaa[n] - 2*pi

plt.subplot(212)
#plt.plot((180/pi)*phi,(180/pi)*aaa,'k')
plt.semilogx((180/pi)*phi,(180/pi)*aaa,'k')
plt.ylabel('/G (degrees)')
plt.axis([1,100, -180,0])
plt.yticks([-90,-45,-180,0])
plt.axvline(5.7,color='k')
plt.axvline(phi2,color='k')
plt.text(2,-90,'$\phi1$ = {}'.format(round(phi1,2)),fontsize=12)
plt.text(30,-140,'$\phi2$ = {}'.format(round(phi2,2)),fontsize=12)
plt.grid(which='both')
plt.xlabel('$\phi$ (degrees)')
plt.savefig('H_zbode.png',dpi=300)
plt.show()


