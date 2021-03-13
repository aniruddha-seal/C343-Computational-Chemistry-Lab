#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 21:56:52 2021

@author: acer
"""

import numpy as np
import matplotlib.pyplot as plt
import math

B=1.99824
v=0

k_B=1.38*(10**(-23))
T=298.15

h=6.626*(10**(-34))
c=3*(10**(10))

n=25

E=np.zeros(n)

for i in range(n):
    E[i]=(0.5*v)+(B*i*(i+1))
    
st=np.zeros(n)
a_st=np.zeros(n)
y_st=np.zeros(n)
y_a_st=np.zeros(n)

for i in range(n-2):
    st[i]=2*B*((2*i)+3)
    y_st[i]=((2*i)+1)*math.exp(-(E[i]*h*c)/(k_B*T))
    if(i%2==0):
        y_st[i]=y_st[i]*1
    else:
        y_st[i]=y_st[i]*2
    
    if(i>=2):
        a_st[i]=-2*B*((2*i)-1)
        y_a_st[i]=((2*i)+1)*math.exp(-(E[i]*h*c)/(k_B*T))
        if(i%2==0):
            y_a_st[i]=y_a_st[i]*1
        else:
            y_a_st[i]=y_a_st[i]*2

st=st+v
y_st=y_st*1
y_a_st=y_a_st*1
a_st=a_st+v
fig = plt.figure(dpi=600)
#plt.xlim(-700,700)
plt.xlabel("Raman shift (cm$^{-1}$)")
plt.ylabel("Intensity")

# Horizontal Bar Plot 
plt.bar(st[0:24], y_st[0:24], width=1.5, color=  'b',label="Stokes") 
plt.bar(a_st[0:24], y_a_st[0:24], width=1.5, color=  'r',label="Anti-Stokes") 
#plt.bar(P_h[0:24], y_P_h[0:24], width=1.2, color=  'r') 
#plt.bar(R_h[0:24], y_R_h[0:24], width=1.2, color=  'r')
plt.legend(loc=1)
plt.tight_layout() 
plt.savefig("raman.png")
# Show Plot 
plt.show()
        
        
    