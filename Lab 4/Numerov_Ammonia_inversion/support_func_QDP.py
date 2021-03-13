import matplotlib.pyplot as plt
import numpy as np
import sys
import time
#code for the functions which find the eigen values
h=0
m=0
Barrier_poten=0#An energy constant which is sufficiently greater than the minimum energy but not too high
#The Barrier_poten should be the energy greater then the first few eigen energies. but should be lees than the energy of the 10th or 20th 
Poten=0
def Constant_feeder(h_m,m_m,BP,V):
    global h,m,Barrier_poten,Poten
    h=h_m
    m=m_m
    Barrier_poten=BP
    Poten=V

def F(x,E):
    v=Poten(x)
    return 2*m*(v-E)/h**2

def Numerov(x,q0,q1,psi1,f1,dx,E):
    
    q2=dx*dx*f1*psi1+2*q1-q0
    f2=F(x+dx,E)
    psi2=q2/(1-f2*dx**2/12)
    if psi2>1e250:
        print("!!!!!!!!!!!!\n!!!LIKELY OVERFLOW EXPECTED IN PSI!!!.\nReduce the upper bound for energy or increase the lower bound")
    return q2,f2,psi2

def Initialize(X,E,psi0=1e-5,psi1=1e-5):
    '''
    Initializes all the variables needed for finding the 
    solution using Numerov Method
    q->phi
    f->(V-E)
    psi->Our solution to the schrodinger equation
    '''
    dx=X[1]-X[0]
    f0=F(X[0],E)
    f1=F(X[1],E)
    q0=psi0*(1-dx**2*f0/12)
    q1=psi1*(1-dx**2*f1/12)
    psi=[psi0,psi1]
    f=[f1]
    q=[q0,q1]
    return psi,f,q

def run_eq(X,e,psi,f,q):
    '''
    The function which does the integration of the
    differential equation
    '''
    for i in range(len(X)-2):
        x=X[i+1]
        f1=f[0]
        psi1=psi[-1]
        q1=q[1]
        q0=q[0]
        dx=X[i+2]-X[i+1]
        q2,f2,psi2=Numerov(x,q0,q1,psi1,f1,dx,e)
        q[0]=q1
        q[1]=q2
        f[0]=f2
        psi.append(psi2)
    psi_n=np.copy(psi)
    return psi_n

def run_mult(X,range_E):
    '''
    For running the solutions for a range of Energies
    '''
    data=[]
    for e in range_E:
        psi,f,q=Initialize(X,e)
        psi_n=run_eq(X,e,psi,f,q)
        data.append([e,psi_n])
    return data

def Explo_min_Finder(X,E_range):
    '''
    Finds the local minimum in explosion factor and returns
    the corresponding energy. Local minima steps works given that
    energy array is sorted othewise no meaning. If energies are not
    pre sorted this function sotrts them 
    '''
    E_range=np.copy(E_range)
    V_x=Poten(X)#V of x{V(x)}
    min_V=min(V_x)
    ind_sub0=np.where(E_range<=min_V)#position where energy is less than minimum potential 
    E_range[ind_sub0]=min_V+Barrier_poten*.1
    E_range=np.sort(E_range)
    active_range=np.where(V_x<=max(E_range))
    ind_lw=active_range[0][0]
    ind_up=active_range[0][-1]
    Sol=run_mult(X,E_range)#getting the solutions ofd the diff equation
    Explo_Farray=[]#explosion-factor array
    Optim_E=[]#the energies near to eigen values
    for s in Sol:
        y=s[1]
        Expec_max=max(y[ind_lw:ind_up])#maximum in the expected region
        Expec_min=abs(min(y[ind_lw:ind_up]))
        if Expec_max<Expec_min:
            Expec_max=Expec_min
        Explo_Farray.append(abs(y[-1]/Expec_max))
    for j in range(len(Sol)-2):
        if Explo_Farray[j]>Explo_Farray[j+1]:
             if Explo_Farray[j+1]<Explo_Farray[j+2]:
                Optim_E.append(Sol[j+1][0])
    return Optim_E

def Dec_Eigen_search(X,E_init,Prec=.1,num_div=10):
    '''
    Decimal Eigen search
    This function takes num_div(number of divisions) number of values in both directions
    Where each value differ by Prec(precision). The local minimas are returned
    '''
    if num_div<0:
        num_div=abs(num_div)
    if num_div==0:
        num_div=10
    num_div=int(num_div)
    E_range=[E_init+i*Prec for i in range(-num_div,num_div+1)]
    return Explo_min_Finder(X,E_range)

def Eigen_Range_finder(X,E_min=0,E_max=10,Prec=5,dE=.1):
    '''
    Given an energy range finds all the eigen energy in that range.
    '''
    E_range=np.arange(E_min,E_max+dE,dE)#The basic range with a division length of div 
    E_range=Explo_min_Finder(X,E_range)#local minimas in the superficial range
    E_dyn=np.copy(E_range)#dynamic range which will change for each loop
    E_temp=np.array([])#temperory loop where the new energies will be stored
    l=0#length of E_temp
    c=1#counter
    for i in range(Prec):
        dE=dE/10#precision adjuster
        j=0#position of each e
        for E_m in E_dyn:#E_m E main
            print("|",end="")
            E=Dec_Eigen_search(X,E_m,dE,9)#calling the Dec_eigen multiples times for zooming into each range
            for e in E:#adding each new energy to E_temp
                c+=1
                if j<l:
                    E_temp[j]=e
                    j+=1
                else:
                    E_temp=np.append(E_temp,e)
                    j=j+1
                    l=l+1 
        print("")
        E_dyn=np.copy(E_temp)
    return E_dyn

def Normalize(x,y,norml_Val=1):#function to normalize the function
    '''
    This function normalizes the y value, such that nintegral
    of over entire X will be one
    the process is simple just finds out the area under the graph
    Absolute of areas are taken otherwise the symetric functions will
    be multiplied by a large number since their are is nearly zero
    here we need to normalise psi square so taking the absolute areas dosen't hurt
    '''
    A=0
    for i in range(len(x)-1):
        dx=(x[i+1]-x[i])
        a=abs(dx*(y[i]+y[i+1])/2)
        A=A+a
    norm_y=y/A
    return norm_y

def Run_both(X,e=.5,psi_0=1e-10,psi_1=1e-10):
    '''
    Runs the solutions in both directions and the use merger
    to mergs both the functions.
    '''
    X_b=X[::-1]
    psi,f,q=Initialize(X,e,psi_0,psi_1)
    psi_b,f_b,q_b=Initialize(X_b,e,psi_0,psi_1)
    run_eq(X,e,psi,f,q)
    run_eq(X_b,e,psi_b,f_b,q_b)
    N_psi=Merger(X,psi,X_b,psi_b)
    return N_psi            

def Merger(X,psi,X_b,psi_b):
    '''
    Merges 2 opppositiely run solutions of the equations 
    and merges them both. The merging is not perfect but the effect
    is imperfections in not significant when looked at ovewrall
    function
    '''
    psi_b=np.copy(psi_b[::-1])
    psi_m=np.copy(psi)# the final merged psi
    last_min=0#stores the postion of last minimum in the psi array this indicates the minimum before the explotion
    max_psi=max(psi)
    min_psi=min(psi)
    if abs(min_psi)>max_psi:
        max_psi=min(psi)
    if max_psi==psi[-1]:#detecting an explosion near the far end
        for p in range(len(psi)-2):
            if abs(psi_m[p])>abs(psi_m[p+1]):
                if abs(psi_m[p+1])<abs(psi_m[p+2]):
                    last_min=p+1
        psi_m[last_min:]=psi_b[last_min:]#merging both of them
        psi_n=Normalize(X,psi_m)
    else:#if no explosion the function is fine
        psi_n=Normalize(X,psi_m)
    return psi_n

def Plot_Eq(X,E_range,ax=plt):
    '''
    Plots all the solutions for the given energies
    '''
    print(f"Plottings for {E_range}")
    i=0
    cn=0#color number
    c=['c','r','g','b','deepskyblue','violet','lime','fuchsia','crimson','darkorchid']
    for e in E_range:
        Psi=Run_both(X,e)
        ax.plot(X,Psi+e,'--',color=c[cn],label=f'Eigen State {i}')
        cn+=1
        if cn==10:
        	cn=0
        i+=1