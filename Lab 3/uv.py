#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 23:16:54 2021

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

x_b,y_b,z_b=np.loadtxt("/home/acer/Downloads/buta_uv.txt",unpack=True)
x_h,y_h,z_h=np.loadtxt("/home/acer/Downloads/hexa_uv.txt",unpack=True)
x_o,y_o,z_o=np.loadtxt("/home/acer/Downloads/octa_uv.txt",unpack=True)

fig=plt.figure(dpi=600)
plt.xlabel("Excitation Energy (nm)")
plt.ylabel("Epsilon")
plt.xlim(0,700)
plt.ylim(0,60000)
plt.plot(x_b,y_b,label="Butadiene")
plt.plot(x_h,y_h,label="Hexatriene")
plt.plot(x_o,y_o,label="Octatetraene")
plt.legend(loc=1)
plt.tight_layout()
plt.savefig("uv.png")
plt.show()