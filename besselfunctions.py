# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 16:09:49 2020

@author: marem
"""

import matplotlib.pyplot as plt
from scipy import signal as sig
import math
import numpy as np

num0 = [1]
den0 = [1]

num1 = [0,1]
den1 = [1,1]

num2 = [0,0,1]
den2 = [1,3,3]

num3 = [0,0,0,1]
den3 = [1,6,15,15]

system0 = sig.lti(num0,den0)
system1 = sig.lti(num1,den1)
system2 = sig.lti(num2,den2)
system3 = sig.lti(num3,den3)

w0, Hmag0, Hphase0 = sig.bode(system0)
w1, Hmag1, Hphase1 = sig.bode(system1)
w2, Hmag2, Hphase2 = sig.bode(system2)
w3, Hmag3, Hphase3 = sig.bode(system3)

Ampl0 = 10**(0.05*Hmag0)
Ampl1 = 10**(0.05*Hmag1)
Ampl2 = 10**(0.05*Hmag2)
Ampl3 = 10**(0.05*Hmag3)

############# Plotting Bode #############################

plt.semilogx(w0,Ampl0)
plt.semilogx(w1,Ampl1)
plt.semilogx(w2,Ampl2)
plt.semilogx(w3,Ampl3)

plt.title('Problem 8.6.2')
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('/_ |H|')
plt.legend()
plt.show()

plt.figure()
plt.semilogx(w0,Hphase0)
plt.semilogx(w1,Hphase1)
plt.semilogx(w2,Hphase2)
plt.semilogx(w3,Hphase3)

plt.title('Problem 8.6.2')
plt.axis([0,10,-200,0])
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('/_ |H|')
plt.legend()
plt.show()
