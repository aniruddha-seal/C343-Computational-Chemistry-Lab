import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import support_func_QDP as EF
def Poten_D(x):
    V=K*(C_4*(x**4)-C_2*(x**2)+C_0)
    return V
print("\t THE NUMERICAL SOLUTION FOR SYMMETRIC POTENTIAL WELL\n")
print("The user has to provide an energy range. This program will try to find an eigen value in this range")
print("For using default values enter '$' when asked for the potential form. \nIf not sure about any values please provide '$', whever you want")
h=input("Enter the values of h-bar:")
if h!='$':
    h=float(h)
else:
    h=1
m=input("Enter the value of m:")   
if m!='$':
    m=float(h)
else:
    m=1
print("\nPOTENTAIL FUNCTION DEFENITION")
print("The potentiAL we use is in the form K*(ax^4-bx^2+C)")
print("\tDO YOU WISH TO PROVIDE a,b,c THEN TYPE'C'\n\t\t OR \n\tSPECIFY THE BARRIER HEIGHT AND DISTANCE OF MINIMAS FROM THE BARRIER")
time.sleep(3)
M=input("If you wish to provide the constants then type 'C'\nif you wish the other method type 'B'.\nOR for a default run type '$' :")
if M=='C':
    C_4=input("Enter the value you want for a:")
    C_2=input("Enter the value you want for b:")
    print("The value of C has to be such that the potential is never below zer, if the entered C is insufficient it will be automatically changed")
    C_0=input("Enter the value you want for C:")
    if C_4=='$':
        C_4=4
    else:
        C_4=float(C_4)
    if C_2=='$':
        C_2=40
    else:
        C_2=float(C_2)
    C_0_cal=C_2**2/(4*C_4)
    if C_0=='$':
        C_0=C_2**2/(4*C_4)
    elif float(C_0)<C_0_cal:
        C_0=C_0_cal
    else:
        C_0=float(C_0)
    K=input("Enter the Scaling factor K:")
    if K=='$':
        K=0.1
    else:
        K=float(K)
    print(f"The potential equation is:{K}({round(C_4,1)}x^4-{round(C_2,1)}x^2+{C_0})")
elif M=='B':
    C_0=input("Enter the barrier height:")
    b=input("Enter the distance from barrier to either one of the minimas:")
    if C_0=='$':
        C_0=10
    else:
        C_0=float(C_0)
    if b=='$':
        b=3
    else:
        b=float(b)
    C_4=C_0/(b**4)
    C_2=2*C_0/(b**2)
    print(f"The potential barrier={C_0}/nThe distance from the barrier to minima={b}")
    print(f"The potential equation is:{round(C_4,2)}x^4-{round(C_2,2)}x^2+{C_0}")
    K=float(input("Enter the scaling factor. Enter one if you want absolute scaling:"))
elif M=='$':
    C_4=4
    C_2=40
    C_0=C_2**2/(4*C_4)
    K=0.1
else:
    print("provide a valid input")
    exit()

E_lw=input("\nEnter the lower bound for searching energy(can be zero but not less than that):")
print(f"\nSETTING THE UPER BOUND TO {K*C_0} WILL GIVE YOU THE FIRST FEW EIGEN STATES,\n!!!!ENERGIES GREATER THAN THIS WILL NOT HAVE SPLITTTING AND WILL BE SIMILAR\nTO THE EIGEN STATES OF QUATUM HARMONIC OSSCILATOR.\nyou can of course set the upper bound of energy higher than this.\nbut if this is too high the eigen energies will be too crowded")
E_up=input(f"\nEnter the upper bound for searching energy:")
if E_lw=='$':
    E_lw=0
else:
    E_lw=float(E_lw)
if E_up=='$':
    E_up=K*C_0
else:
    E_up=float(E_up)
Xmin=-((C_2+(C_2**2-4*C_4*(C_0-E_up*10))**.5)/(2*C_4))**(.5)*1.5
Xmax=abs(Xmin)
print(f"\nThe region on X axis starting from {round(Xmin,2)} to {round(Xmax,2)} is considered.\nAnd evaluation of psi is done only in this domain")
print(f"The potential equation is:{K}({round(C_4,1)}x^4-{round(C_2,1)}x^2+{round(C_0,1)})")
print("\n!!PLEASE WAIT!! \n!!CALCULATING THE EIGEN ENERGIES!!")
T_b=time.time()
X=np.linspace(Xmin,Xmax,10**4)
Barrier_poten=K*C_0#the middle potential barrier
if .1>Barrier_poten*.01:
    dE=Barrier_poten/100
else:
    dE=.1
EF.Constant_feeder(h,m,Barrier_poten,Poten_D)
Eigen_E=EF.Eigen_Range_finder(X,E_lw,E_up,10,dE)
T_f=time.time()
print('\nTime taken for computing:',T_f-T_b)
print("The eigen energies are :")
i=1
if len(Eigen_E)>0:
    for e in Eigen_E:
        print(f"Energy  of Eigen state {i}:{e}")
        i+=1
    fig=plt.figure()
    ax=plt.axes(xlabel='X',ylabel='Psi',ylim=(0,Eigen_E[-1]+Barrier_poten*.1))
    X=np.linspace(Xmin,Xmax,10**4)
    EF.Plot_Eq(X,Eigen_E,ax)
    ax.plot(X,Poten_D(X),color='b',label='Potential')
    plt.legend()
    plt.show()
else:
    print("No eigen energies")
