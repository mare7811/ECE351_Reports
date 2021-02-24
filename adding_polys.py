# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 22:44:46 2020

@author: marem
"""

import numpy as np

def poly_mult(A,B):
    prod = np.zeros(np.size(A)+np.size(B)-1)   
    
    for i in range(0,np.size(A)):
        for j in range(0,np.size(B)):
            prod[i+j] += A[i] * B[j]
    return prod

A = np.array([1,1])
B = np.array([1,1])
C = poly_mult(A,B)
D = poly_mult(C,A)
E = poly_mult(D,A)
print(E)

s1 = 1
s2 = 1.45*A
s3 = 2.05*C
s4 = 1.39*D
s5 = .529*E
print(s1,s2,s3,s4,s5)

den = np.zeros(5)
den[0] = s1+s2[0]+s3[0]+s4[0]+s5[0]
den[1] = s2[1]+s3[1]+s4[1]+s5[1]
den[2] = s3[2]+s4[2]+s5[2]
den[3] = s4[3]+s5[3]
den[4] = s5[4]
print(den)

num = .514*E
print(num)