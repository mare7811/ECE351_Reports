# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 23:05:23 2020

@author: marem
"""

import matplotlib.pyplot as plt
import numpy as np

NN = 5000
phi = np.linspace(0, 2*np.pi, NN)
z = np.zeros(NN, dtype = np.complex)
H = np.zeros(NN, dtype = np.complex)

T = 0.2
alpha = 0.5
omega = 0.866
x = np.deg2rad(T*omega)
k1 = 0.577

for n in range(0, NN):
    z = np.exp(1j*phi[n])
    H[n] = ((z/(z-np.exp(-T))) + ((k1*np.exp(-alpha*T)*np.sin(x)*z) / (z**2-2*np.exp(-alpha*T)*np.cos(x)*z+np.exp(-T))) - ((z**2-np.exp(-alpha*T)*np.cos(x)*z) / (z**2-2*np.exp(-alpha*T)*np.cos(x)*z+np.exp(-T)))) / -4.808



plt.figure(figsize=(6,8))       
plt.subplot(211)
#plt.plot((180/pi)*phi,abs(H),'k')
plt.semilogx((180/np.pi)*phi, 20*np.log10(H), 'k')
plt.title('My_zbode')
#plt.axis([1, 100, -20, 10])
plt.ylabel('|G| dB')
#plt.yticks([-20, -3, 0])
plt.grid(which='both')
plt.axvline(5.7, color='k')
#plt.text(2.8, -10, '$\phi$ = {}'.format(round(5.7, 1)), fontsize=12)

aaa = np.angle(H)
#for n in range(NN):
#    if aaa[n] > 0:
#        aaa[n] = aaa[n] - 2*pi

plt.figure(figsize=(6,8))
plt.subplot(212)
#plt.plot((180/pi)*phi,(180/pi)*aaa,'k')
plt.semilogx((180/np.pi)*phi, (180/np.pi)*aaa, 'k')
#plt.axis([1, 100, -90, 0])
plt.ylabel('/G (degrees)')
#plt.yticks([-90, -45, 0])
plt.xlabel('$\phi$ (degrees)')
plt.grid(which='both')
#plt.axvline(5.7, color='k')
#plt.text(2.8, -70, '$\phi$ = {}'.format(round(5.7, 2)), fontsize=12)
plt.show()
