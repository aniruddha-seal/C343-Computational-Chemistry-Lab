#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 10:08:51 2021

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

E=3.5
n=100

p=np.zeros(1000)
x=np.zeros(1000)
g=np.zeros(1000)

x[0]=-5
dx=0.1
p[0]=0
p[1]=-0.0001
x[1]=x[0]+dx
g[0]=x[0]**2-2*E
g[1]=x[1]**2-2*E


for i in range(1,n):
    x[i+1]=x[i]+dx
    g[i+1]=x[i+1]**2-2*E
    p[i+1]=(2*p[i]-p[i-1]+(5*g[i]*p[i]*(dx**2))/6+(g[i-1]*p[i-1]*(dx**2))/12)/(1-(g[i+1]*(dx**2))/12)
    
fig=plt.figure(dpi=600)
plt.xlabel("$x$")
plt.ylabel("$\Psi(x)$")
#plt.ylim(0.005,27)
plt.scatter(x[0:101],p[0:101])
plt.plot(x[0:101],p[0:101])
plt.tight_layout()
plt.savefig("shm_4.png")
plt.show()