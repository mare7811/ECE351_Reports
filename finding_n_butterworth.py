# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:07:01 2020

@author: marem
"""

import numpy as np 
from math import log10

wp = 1/1.3
Hp2 = .9

np = log10( 1/Hp2 - 1)/(2*log10(wp) )
print('np = ',np)



ws = 2/1.3
Hs2 = .15

ns = log10( 1/Hs2 - 1)/(2*log10(ws) )
print('ns = ',ns)
