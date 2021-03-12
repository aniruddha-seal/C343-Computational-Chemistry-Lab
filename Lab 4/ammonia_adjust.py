#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 17:34:29 2021

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


x,V=np.loadtxt("/home/acer/Downloads/final.txt",unpack=True)


x=x-90

x_p=np.linspace(-40,40,500)
x_chi=np.linspace(-40,40,39)
f_chi=2.84015e-08*x_chi**4 -3.14882e-05*x_chi**2 - 56.374
f=2.84015e-08*x_p**4 -3.14882e-05*x_p**2 - 56.374
MyFile=open('ammonia.txt','w')
k=0
E=0
for element in V:
    er = (f_chi[k]-V[k])**2
    er = er/f_chi[k]
    E=E+er
    MyFile.write(str(x[k])+"      "+str(V[k]))
    MyFile.write('\n')
    k=k+1
    print(k)
MyFile.close()

print(E)
    
fig=plt.figure(dpi=600)
plt.xlabel("x: X-N-H Angle - 90 (in degree)")
plt.ylabel("V(x) (in Hartree)")
plt.scatter(x,V,label="Data points",color='orange')
plt.plot(x_p,f,label="Fitted Function",color='red')
plt.legend(loc=9)
plt.tight_layout()
plt.savefig("fig_1.png")
plt.show()