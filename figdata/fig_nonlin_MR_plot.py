import numpy as np
import matplotlib.pyplot as plt
from numpy import log10 as lg
from numpy import pi as pi
from scipy.interpolate import interp1d as sp_interp1d
from scipy.integrate import odeint
from scipy.integrate import ode
import warnings
import timeit
import scipy.optimize as opt
from matplotlib import cm
from astropy import constants as const
from astropy import units as u
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['axes.labelpad'] = 8.0
plt.rcParams['figure.constrained_layout.h_pad'] = 0
plt.rcParams['text.usetex'] = True
plt.rc('text', usetex=True)
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.tick_params(axis='both', which='minor', labelsize=18)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
import warnings
import timeit

import scipy.optimize
from matplotlib import cm
from scipy.integrate import ode as sp_ode

t0 = timeit.time.time()
from shapely.geometry import LineString
from scipy.interpolate import UnivariateSpline
G=const.G.cgs.value
c=const.c.cgs.value
MSUN=const.M_sun.cgs.value
hbar=const.hbar.cgs.value
m_n=const.m_n.cgs.value
KM=10**5
import math

G=const.G.cgs.value
c=const.c.cgs.value
Ms=const.M_sun.cgs.value
hbar=const.hbar.cgs.value
m_n=const.m_n.cgs.value
km=10**5

label=['WFF1','SLy4','AP4', 'MPA1','PAL1']
fig, axs = plt.subplots(2, 4,figsize=(24,12))
colorset=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple'] 
plt.setp(axs, yticks=[0, 0.5, 1,1.5,2,2.5],xticks=[9,11,13,15])

for i in range(2):
    for j in range(4):
        for k in range(5):
            datagr=np.genfromtxt('TOV_5eqs_data'+str(k+1)+'.txt')
            M=datagr[:,1]
            R=datagr[:,2]
            if i==0 and j==0:
              axs[i,j].plot(R/10**5,M/Ms,color=colorset[k],linestyle='--',linewidth=2, label=label[k])
            else:
              axs[i,j].plot(R/10**5,M/Ms,color=colorset[k],linestyle='--',linewidth=2)
            axs[i,j].set_ylim(0,2.6)
            axs[i,j].set_xlim(8.3,15)
            axs[i,j].minorticks_on()
            
# the first plot
for i in range(5):
    data=np.genfromtxt('stgb_solver_pdata2'+str(i+1)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        axs[0,0].plot(b321[ntrimset[ii]:ntrimset[ii+1]]/KM, \
                     b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[i],markersize=5)  

# the second plot
list2=[26,27,28,30]
for i in range(4):
    a=list2[i]
    data=np.genfromtxt('stgb_solver_pdata'+str(list2[i])+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        if a!=30:
            axs[0,1].plot(b321[ntrimset[ii]:ntrimset[ii+1]]/KM, \
            b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[i],markersize=5) 
        else:
            axs[0,1].plot(b321[ntrimset[ii]:ntrimset[ii+1]]/KM, \
            b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[i+1],markersize=5) 

# the third plot

for i in range(5):
    data=np.genfromtxt('stgb_solver_pdata'+str(31+i)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        x=b321[ntrimset[ii]:ntrimset[ii+1]]/KM
        y=b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN
        index=y.argsort()
        ydata=y[index]
        xdata=x[index]
        s1 = UnivariateSpline(ydata, xdata, s=5)
        ys=np.linspace(min(ydata),max(ydata),5)
        xs=s1(ys)
        axs[0,2].plot(xs,ys, color = colorset[i],markersize=5,linewidth=3) 
        
# the fourth plot

for i in range(5):
    data=np.genfromtxt('stgb_solver_pdata'+str(36+i)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        x=b321[ntrimset[ii]:ntrimset[ii+1]]/KM
        y=b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN
        index=y.argsort()
        ydata=y[index]
        xdata=x[index]
        s1 = UnivariateSpline(ydata, xdata, s=5)
        ys=np.linspace(min(ydata),max(ydata),5)
        xs=s1(ys)
        axs[0,3].plot(xs,ys, color = colorset[i],markersize=5,linewidth=3)

# the fifth plot
for i in range(5):
    data=np.genfromtxt('stgb_solver_pdata'+str(i+1)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        axs[1,0].plot(b321[ntrimset[ii]:ntrimset[ii+1]]/KM, \
                     b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[i],markersize=5) 
        
# the sixth plot
for i in range(5):
    data=np.genfromtxt('stgb_solver_pdata'+str(i+6)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        axs[1,1].plot(b321[ntrimset[ii]:ntrimset[ii+1]]/KM, \
                     b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[i],markersize=5) 
        
# the seventh plot
for i in range(5):
    data=np.genfromtxt('stgb_solver_pdata'+str(11+i)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        x=b321[ntrimset[ii]:ntrimset[ii+1]]/KM
        y=b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN
        index=y.argsort()
        ydata=y[index]
        xdata=x[index]
        s1 = UnivariateSpline(ydata, xdata, s=5)
        ys=np.linspace(min(ydata),max(ydata),5)
        xs=s1(ys)
        axs[1,2].plot(xs,ys, color = colorset[i],markersize=5,linewidth=3)
        
# the eighth plot
for i in range(5):
    data=np.genfromtxt('stgb_solver_pdata'+str(16+i)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=data[:, 0], \
    data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], \
    data[:, 6], data[:, 7], data[:, 8], data[:, 9] 
    ecstep = 0.6e+15
    ntrimset = np.array([0])
    for ii in range(0, len(b021)-1):
        testecsep = b021[ii+1]-b021[ii]
        if testecsep > ecstep:     
            ntrimset = np.append(ntrimset, ii+1)
    ntrimset = np.append(ntrimset, len(b021))
    for ii in range(0, len(ntrimset)-1):
        x=b321[ntrimset[ii]:ntrimset[ii+1]]/KM
        y=b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN
        index=y.argsort()
        ydata=y[index]
        xdata=x[index]
        s1 = UnivariateSpline(ydata, xdata, s=5)
        ys=np.linspace(min(ydata),max(ydata),10)
        xs=s1(ys)
        axs[1,3].plot(xs,ys, color = colorset[i],markersize=5,linewidth=3)
        
# plot decorations and labels

axs[0,0].legend(fontsize=13, frameon=False)


fig.text(0.08, 0.5, r'$M\,[{ M_{\odot}}]$', ha='center', fontsize=30,rotation='vertical')
fig.text(0.48, 0.05, r'$R\,[\rm km]$' ,fontsize=30)      
fig.text(0.91, 0.7, r'$a=1$' ,fontsize=30, rotation='-90')          
fig.text(0.91, 0.27, r'$a=0.1$' ,fontsize=30, rotation='-90')     
fig.text(0.18, 0.9, r'$b=-1$' ,fontsize=30)     
fig.text(0.38, 0.9, r'$b=-0.2$' ,fontsize=30)     
fig.text(0.58, 0.9, r'$b=3$' ,fontsize=30)   
fig.text(0.78, 0.9, r'$b=10$' ,fontsize=30)  
plt.savefig("fig_nonlin_MR.pdf", format='pdf', bbox_inches="tight")

plt.show()
