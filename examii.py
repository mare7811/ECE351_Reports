# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450 
Exam II
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sig
#from math import pi, exp, cos, sin, log, sqrt
from control import margin
#from control import tf
from math import pi, sqrt, exp, cos,sin
#from cmath import exp

# H(s)
K = 400
n1 = [0,0,K*10]
d1 = [1,3,-10]

#Phase Lead
phim = 40 * 0.0174533 #degrees to radians
alpha = (1+sin(phim))/(1-sin(phim))
wm = 100   
wz = wm/sqrt(alpha)
wp = wm*sqrt(alpha)
#n2 = [1,wz]
#d2 = [1,wp]
#k = wp/wz
n2 = [1,2*wz,wz**2]
d2 = [1,2*wp,wp**2]
k = (wp/wz)**2
print(sqrt(k),wz,wp)
n12 = k*np.convolve(n1,n2)
d12 = np.convolve(d1,d2)


#Phase Lag
n3 = [1,1]
d3 = [1,0]
n13 = np.convolve(n1,n3)
d13 = np.convolve(d1,d3)
n123 = np.convolve(n12,n3)
d123 = np.convolve(d12,d3)


#numerator and denominator of G(s)H(s)
num = n123
den = d123

#Calculating values for Bode plots
w = np.linspace(0.01,10000,10000)
system = sig.lti(num,den)
w, Hmag, Hphase = sig.bode(system,w)
gm, pm, wg, wp = margin(Hmag,Hphase,w)
# wp  freq for phase margin at gain crossover (gain = 1)
# pm phase maring


############### Plotting Magnitude Bode ######################
plt.subplot(211)
plt.semilogx(w,Hmag,'k')
plt.semilogx(w,Hmag,'k')
#plt.axis([ 0, 1e2, -25, 50])
#plt.xticks([1,10,30,100,1000])
plt.ylabel('|H| dB',size = 12)
plt.text(1,0,'$\omega$p = {}'.format(round(wp,1)),fontsize=12)
plt.title('Bode Comp')
plt.grid(which='both')

for n in range(100):
    if Hphase[n] > 0:
        Hphase[n] = Hphase[n] - 360


################## Plotting Phase Bode ###########################
plt.subplot(212)
plt.semilogx(w,Hphase,'k')
#plt.axis([ 1, 1e2, -180,0])
#plt.yticks([-180,-90,0])
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('Phase (degrees)',size=12)
plt.text(.1,-150,'pm = {}'.format(round(pm,0)),fontsize=12)
plt.grid(which='both')
plt.show()


################### Time Portion ###############################
dt = 0.001
NN = 10000 
TT = np.arange(0,NN*dt,dt)
step = np.zeros(NN)
ramp = np.zeros(NN)
parabola = np.zeros(NN)
errS = np.zeros(NN)
errR = np.zeros(NN)
errP = np.zeros(NN)

for i in range(NN):
    step[i] = 1.0
    ramp[i] = (dt*i)
    parabola[i] = (dt*i)**(2)
    
denCL = np.add(num,den)

t1, y1, x1 = sig.lsim((num,denCL),step,TT)
t2, y2, x2 = sig.lsim((num,denCL),ramp,TT)
t3, y3, x3 = sig.lsim((num,denCL),parabola,TT)

for i in range(NN):
    errS[i] = step[i] - y1[i]
    errR[i] = ramp[i] - y2[i]
    errP[i] = parabola[i]  - y3[i]    


####### Step Time response #########################
plt.subplot(321)
plt.plot(TT,y1,'b--',label='y1(t)')
plt.plot(TT,step,'k',label='u(t)')
plt.axis([0,1,0,1.5])
plt.ylabel('step')
plt.xlabel('t (sec)')
#plt.yticks([0,.9,1.1,1.5])
#plt.legend()
plt.grid()

####### Step error response #####################
plt.subplot(322)
plt.plot(TT,errS,'k',label='error')
plt.legend()
plt.axis([0,1,-.02,.02])
#plt.yticks([0,0.02,.05,.1])
plt.grid()
plt.savefig('position.png')
plt.show()

########### Ramp Time Response #######################
plt.subplot(321)
plt.plot(TT,y2,'k--',label='y2(t)')
plt.plot(TT,ramp,'k',label='r(t)')
plt.xlabel('t (sec)')
plt.ylabel('ramp')
plt.legend()
plt.grid()

############ Ramp error Response #######################
plt.subplot(322)
plt.plot(TT,errR,'k',label='error')
plt.legend(loc=4)
#plt.text(.3,-.1,'K = {}'.format(round(K,0)),fontsize=12)
plt.xlabel("t (sec)") 
#plt.axis([0,1,-.030,.03])
#plt.yticks([-.5,0,.3,.5])
plt.grid()
plt.savefig('velocity.png')
plt.show()



######## Finding frequency at given H ###################
H_index = 0
H = -10*np.log10(alpha)
for i in range(NN-1):
    if(Hmag[i] > H and Hmag[i+1] < H):
        H_index = i

print('|H| = ',Hmag[H_index],'dB','\n\u03c9 = ',w[H_index],'rad/s')


