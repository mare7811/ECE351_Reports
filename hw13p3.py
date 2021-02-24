# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450 
HW 13
Problem 9.3.1
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy





dt = 0.01
NN = 1000
k = np.arange(0,NN*dt,dt)

x = .5*np.sin(2*np.pi*.1*k) + np.sin(2*np.pi*k) + .5*np.sin(2*np.pi*10*k)

plt.figure()
plt.plot(k,x,'k')
plt.title('input signal')
plt.ylabel('x[k]')
plt.xlabel('sec')
plt.grid()
plt.ylabel('f')


f = scipy.fft.fft(x)
H = 300*np.exp(-((k-1)**2) / (2* .2**2))

plt.figure()
plt.plot(k,np.abs(f),'k')
plt.title('A one sided Gaussian and FFT')
plt.plot(k,H)
plt.ylabel('X{w} / H(w)')
plt.xlabel('freq [Hz]')
plt.grid()
plt.ylabel('f')


h = scipy.fft.ifft(H)

plt.figure()
plt.plot(k,(h))
plt.title('Inverse FFT of H')
plt.ylabel('X{w} / H(w)')
plt.xlabel('freq [Hz]')
plt.grid()
plt.ylabel('f')

hh = np.zeros(NN)
M = 20
for n in range(M):
    hh[n+M] = h[n].real
    hh[M-n] = hh[M+n]

plt.figure()
plt.plot(k,hh)
plt.title('first Few points of h[k]')
#plt.axis([0,.4,-10,10])
plt.ylabel('X{w} / H(w)')
plt.xlabel('freq [Hz]')
plt.grid()
plt.ylabel('f')



l = slice(0,1)
y = np.convolve(x,h[l])

plt.figure()
plt.plot(k,y)
plt.title('Convolution of x and h')
plt.ylabel('X{w} / H(w)')
plt.xlabel('freq [Hz]')
plt.grid()
plt.ylabel('f')


Y = scipy.fft.fft(y)

plt.figure()
plt.plot(k,np.abs(Y))
#plt.axis([0,5,0,100])
plt.title('FFT of y from previous plot')
plt.ylabel('X{w} / H(w)')
plt.xlabel('freq [Hz]')
plt.grid()
plt.ylabel('f')
plt.show()