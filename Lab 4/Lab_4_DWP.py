#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 10:44:44 2021

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


a=0.02
b=-3.14882e-05
c=2.84015e-08

def V(x):
    v=a+b*(x**2)+c*(x**4)
    return 2*v


#E=0.01609
E=0.01685
n= 160

p=np.zeros(n)
x=np.zeros(n)
g=np.zeros(n)
pot=np.zeros(n)

x[0]=-40
dx=0.5
p[0]=0
p[1]=0.0001
x[1]=x[0]+dx
g[0]=V(x[0])-2*E
g[1]=V(x[1])-2*E

pot[0]=V(x[0])
pot[1]=V(x[1])

for i in range(1,n-1):
    x[i+1]=x[i]+dx
    g[i+1]=V(x[i+1])-2*E
    p[i+1]=(2*p[i]-p[i-1]+(5*g[i]*p[i]*(dx**2))/6+(g[i-1]*p[i-1]*(dx**2))/12)/(1-(g[i+1]*(dx**2))/12)
    pot[i+1]=V(x[i+1])

fig=plt.figure(dpi=600)
plt.xlabel("$x$")
plt.ylabel("$\Psi(x)$")
#plt.xlim(-40,40)
#plt.scatter(x,pot)
#plt.ylim(0.1,27)
plt.scatter(x,p,c='red')
plt.tight_layout()
#plt.savefig("dwp_2.png")
plt.show()
