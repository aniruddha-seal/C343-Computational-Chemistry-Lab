import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import support_func_QDP as EF
def Poten_D(x):
    V=K*(C_4*(x**4)+C_2*(x**2)+C_0)
    return V
h,m=1,1
C_4=2.84015e-08
C_2=-3.14882e-05
C_0=0.02
K=1
T_b=time.time()
Xmin,Xmax=-45,45
E_lw,E_up=0.01,.02
dE=.0001
X=np.linspace(Xmin,Xmax,10**4)

'''
fig=plt.figure(dpi=600)
plt.xlabel("x")
plt.ylabel("V'(x)")
plt.plot(X,Poten_D(X))
plt.show()
'''
Barrier_poten=K*C_0#
EF.Constant_feeder(h,m,Barrier_poten,Poten_D)
Eigen_E=EF.Eigen_Range_finder(X,E_lw,E_up,10,dE)
T_f=time.time()

print('\nTime taken for computing:',T_f-T_b)
i=0
print("The eigen energies are :")
if len(Eigen_E)>0:
    for e in Eigen_E:
        print(f"Energy  of Eigen state {i}:{e}")
        i+=1
    fig=plt.figure(dpi=600)
    plt.xlim(Xmin,Xmax)
    plt.ylim(-1,+1)
    ax=plt.axes(xlabel='$\phi$ (X-N-H)',ylabel='$\Psi$',ylim=(-0.01,+0.04))
    X=np.linspace(Xmin,Xmax,10**4)
    EF.Plot_Eq(X,Eigen_E,ax)
    #ax.plot(X,Poten_D(X),color='b',label='Potential')
    plt.legend()
    plt.tight_layout()
    plt.show()
else:
    print("No eigen energies")

