#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 15:37:04 2021

@author: acer
"""

from scipy import optimize
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from sys import argv
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from scipy.stats import chisquare

mpl.rcParams['font.family'] = 'FreeSerif'
plt.rcParams['font.size'] = 12

B=10.5933002
v=2990.9248

k_B=1.38*(10**(-23))
T=298.15

h=6.626*(10**(-34))
c=3*(10**(10))

n=20

P=np.zeros(n)
R=np.zeros(n)

y_P=np.zeros(n)
y_R=np.zeros(n)

B_h=10.577394
v_h=2988.6785

P_h=np.zeros(n)
R_h=np.zeros(n)

y_P_h=np.zeros(n)
y_R_h=np.zeros(n)

t=0
s=0

for i in range(n):
    P_h[i] = v_h + (2*B_h*(i+1))
    P[i] = v + (2*B*(i+1))
    R_h[i] = v_h - (2*B_h*(i+1))
    R[i] = v - (2*B*(i+1))

    k = (v*0.5)+(B*i*(i+1))
    k_h = (v_h*0.5)+(B_h*i*(i+1))
    y_P_h[i] = ((2*i)+1)*math.exp((-k_h*h*c)/(k_B*T))*24.23
    y_P[i] = ((2*i)+1)*math.exp((-k*h*c)/(k_B*T))*75.77
    y_R_h[i] = ((2*i)+1)*math.exp((-k_h*h*c)/(k_B*T))*24.23
    y_R[i] = ((2*i)+1)*math.exp((-k*h*c)/(k_B*T))*75.77
    
y_P_h=y_P_h*10
y_P=y_P*10
y_R=y_R*10
y_R_h=y_R_h*10

fig = plt.figure(dpi=650)
plt.xlabel("Frequency (cm$^{-1}$)")
plt.ylabel("Intensity")

# Horizontal Bar Plot 
plt.bar(P[0:24], y_P[0:24], width=1.2, color=  'b',label="H35Cl (75.77%)") 
plt.bar(R[0:24], y_R[0:24], width=1.2, color=  'b') 
plt.bar(P_h[0:24], y_P_h[0:24], width=1.2, color=  'r',label="H37Cl (24.23%)") 
plt.bar(R_h[0:24], y_R_h[0:24], width=1.2, color=  'r') 
plt.legend(loc=1)
plt.tight_layout()
plt.savefig("ro_vib_2.png")
# Show Plot 
plt.show()

'''
plt.scatter(P,y_P,c='blue')
plt.scatter(R,y_R,c='blue')
plt.show()
'''