# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import signal as sig
from math import pi, exp, cos, sin, log, sqrt
from control import margin
from control import tf
from cmath import exp

NN = 5000
phi = np.linspace(0,2*pi,NN)
z = np.zeros(NN, dtype = np.complex)
H = np.zeros(NN, dtype = np.complex)

a = 0
b = .01
sin_i = sin(a)*math.cosh(b) + 1j*cos(a)*math.sinh(b)
cos_i = cos(a)*math.cos(b) - 1j*sin(a)*math.sinh(b)


for n in range(0,NN):
    z = exp(1j*phi[n])
    H[n] = 6.3 * ( (2.92*(1j)*sin_i*z + (z**2)-.97*cos_i*z) / ((z**2)-1.94*cos_i*z+.94))
    

plt.subplot(211)
#plt.plot((180/pi)*phi,abs(H),'k')
plt.semilogx((180/pi)*phi,20*np.log10(H),'k')
#plt.axis([1,100, -20, 10])
plt.ylabel('|G| dB')
plt.yticks([ -20,-3,0])
plt.axvline(5.7,color='k')
plt.text(2.8,-10,'$\phi$ = {}'.format(round(5.7,1)),fontsize=12)
plt.title('My_zbode')
plt.grid(which='both')

aaa = np.angle(H)
#for n in range(NN):
#    if aaa[n] > 0:
#        aaa[n] = aaa[n] - 2*pi

plt.subplot(212)
#plt.plot((180/pi)*phi,(180/pi)*aaa,'k')
plt.semilogx((180/pi)*phi,(180/pi)*aaa,'k')
plt.ylabel('/G (degrees)')
plt.axis([1,100, -90,0])
plt.yticks([-90,-45,0])
plt.axvline(5.7,color='k')
plt.text(2.8,-70,'$\phi$ = {}'.format(round(5.7,2)),fontsize=12)
plt.grid(which='both')
plt.xlabel('$\phi$ (degrees)')
plt.savefig('H_zbode.png',dpi=300)
plt.show()


print(z)