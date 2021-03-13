#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:17:08 2021

@author: acer
"""

import numpy as np
import matplotlib.pyplot as plt
import math

B=1044.5
v=2847227

k_B=1.38*(10**(-23))
T=298.15


h=6.626*(10**(-34))
c=3*(10**(10))
E_0=np.zeros(25)
E_1=np.zeros(25)

for j in range(25):
    E_0[j]=(v*0.5)+(B*j*(j+1))
    E_1[j]=(v*1.5)+(B*j*(j+1))
    print(j)
    
P=np.zeros(25)
R=np.zeros(25)

y_P=np.zeros(25)
y_R=np.zeros(25)

m=n=0
for j in range(25):
    if(j<=23):
        R[m]=E_1[j+1]-E_0[j]
        y_R[m]=((2*j)+1)*math.exp(-(h*c*E_0[j])/(k_B*T))
        print(m)
        m=m+1
    if(j>=1):
        P[n]=E_1[j-1]-E_0[j]
        y_P[m]=((2*j)+1)*math.exp(-(h*c*E_0[j])/(k_B*T))
        
fig = plt.figure()
plt.plot(P,y_P,c='blue')
plt.plot(R,y_R,c='blue')
plt.show()
        
    