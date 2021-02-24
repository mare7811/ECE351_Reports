# -*- coding: utf-8 -*-
"""
Miguel Mares
ECE 450
HW 10
Problem 8.3.3
"""

import math
import numpy as np

Hp = .95
Hs = .07
ws = 1/(600/800)
eps = np.sqrt((1/Hp**2)-1)
print('eps = ',eps)

n1 = math.acosh((1/eps)*np.sqrt((1/Hs**2)-1))
n2 = 1/math.acosh(ws)
n = n1*n2
print('n = ',n)
print('The minimum-order Chebyshev filter must be n=3 order. ')