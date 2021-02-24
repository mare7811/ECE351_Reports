# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
Exam 4
Part B
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sig
from math import pi, exp, cos, sin, log, sqrt
import pandas as pd

dt = 0.05
x = np.arange(-3,8,dt)
mu = 2.5
sigma = .7
a = 1/(sigma*sqrt(2*pi))

def y(k):
    y = 2*a*exp(-(k-mu)**2/(2*sigma**2))
    return y

H = np.zeros(np.size(x))
for n in range(np.size(x)):
    H[n] = y(x[n])

# for n in range(1,int(np.size(x)/2)-1):
#     H[50-n] = H[n]  
    

plt.plot(x,H)  # Plot amplitude, not dB.
plt.title('Gaussian Function')
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('H(w)')
plt.axis([-0,5,0,1.2])
plt.yticks([0,.15,.9,1])
plt.xticks([1,2,2.5,3,4])
plt.show()


# for n in range(1,100):
#     g[100-n] = g[n]  
    

h = np.fft.ifft(H)

plt.plot(x,h.real)  # Plot amplitude, not dB.
plt.title('Inverse FFT of Guassian')
plt.grid(which='both')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('T (sec)')
plt.ylabel('h[k]')
plt.xlabel('T (sec)')
#plt.axis([-1,-.8,-.2,.3])
#plt.xticks([100000,200000,50000])
plt.show()


M = 8
hh = np.zeros(2*M+1)
k = np.arange(0,2*M+1)

""" Move the filter to the left side """
for n in range(M):
    hh[n+M] = h[n].real
    hh[M-n] = hh[M+n]


plt.plot(hh)
plt.plot(hh,'ok')
plt.axis([0 ,2*M,-.3,.3])
plt.xlabel('Freq (Hz)')
plt.ylabel('hh[k]')
plt.grid()
plt.show()

df = pd.DataFrame({"k":k,"hh":hh})
df.rename(columns={"hh":"k","B":"hh[k]"})
print(df)