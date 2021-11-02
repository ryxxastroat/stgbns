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

fig, axs = plt.subplots(2, 4,figsize=(24,12),linewidth=2)
colorset=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple'] 
#plt.setp(axs, xticks=[0, 0.5, 1,1.5,2,2.5])

for i in range(2):
    for j in range(4):
        for k in range(5):
            datagr=np.genfromtxt('TOV_5eqs_data'+str(k+1)+'.txt')
            M=datagr[:,1]/MSUN
            I=datagr[:,4]/1e45 
            axs[i,j].plot(M,I,color=colorset[k],linestyle='--',linewidth=0.8)
            #axs[i,j].set_ylim(24,62)
            #axs[i,j].set_xlim(0.4,2.6)
            axs[i,j].minorticks_on()

       
# the first plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+21)+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45
    axs[0,0].plot(x21, y21, color = colorset[i], linewidth=3) 
    
  
# the second plot
list2=[26,27,28,30]
for i in range(4):
    a=list2[i]
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(list2[i])+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45 
    if a!=30:
        axs[0,1].plot(x21, y21, color = colorset[i], linewidth=5) 
   
    else:
        axs[0,1].plot(x21, y21, color = colorset[i+1], linewidth=7)  

# the third plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(31+i)+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45
    #index=x21.argsort()
    #ydata=y21[index]
    #xdata=x21[index]
    #s1 = UnivariateSpline(xdata, ydata, s=5)
    #xs=np.linspace(min(xdata),max(xdata), 20)
    #ys=s1(xs)
    axs[0,2].plot(x21, y21, color = colorset[i], linewidth=3) 
         
# the fourth plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(36+i)+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45
    #index=x21.argsort()
    #ydata=y21[index]
    #xdata=x21[index]
    #s1 = UnivariateSpline(xdata, ydata, s=5)
    #xs=np.linspace(min(xdata),max(xdata),5)
    #ys=s1(xs)
    axs[0,3].plot(x21, y21, color = colorset[i], linewidth=3) 


# the fifth plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(1+i)+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45
    axs[1,0].plot(x21, y21, color = colorset[i], linewidth=3) 
      
# the sixth plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(6+i)+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45
    axs[1,1].plot(x21, y21, color = colorset[i], linewidth=5)  
    
# the seventh plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(11+i)+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45
    #index=x21.argsort()
    #ydata=y21[index]
    #xdata=x21[index]
    #s1 = UnivariateSpline(xdata, ydata, s=5)
    #xs=np.linspace(min(xdata),max(xdata),5)
    #ys=s1(xs)
    axs[1,2].plot(x21, y21, color = colorset[i], linewidth=3) 
      
# the eighth plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(16+i)+'.txt')
    c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, \
    c1021=data[:, 0], data[:, 1], data[:, 2], data[:, 3], \
    data[:, 4], data[:, 5], data[:, 6], data[:, 7], \
    data[:, 8], data[:, 9], data[:, 10] 
    x21 = c321/MSUN 
    y21 = c421/1e45
    #index=x21.argsort()
    #ydata=y21[index]
    #xdata=x21[index]
    #s1 = UnivariateSpline(xdata, ydata, s=5)
    #xs=np.linspace(min(xdata),max(xdata),5)
    #ys=s1(xs)
    axs[1,3].plot(x21, y21, color = colorset[i], linewidth=3)
           
            
       
                      
fig.text(0.08, 0.40,r'$I[\rm 10^{45} \, g\,cm^2]$'  , ha='center', fontsize=30,rotation='vertical')
fig.text(0.48, 0.05, r'$M\,[{ \rm M_{\odot}}]$',fontsize=30)      
fig.text(0.91, 0.7, r'$a=1$' ,fontsize=30, rotation='-90')          
fig.text(0.91, 0.27, r'$a=0.1$' ,fontsize=30, rotation='-90')     
fig.text(0.18, 0.9, r'$b=-1$' ,fontsize=30)     
fig.text(0.38, 0.9, r'$b=-0.2$' ,fontsize=30)     
fig.text(0.58, 0.9, r'$b=3$' ,fontsize=30)   
fig.text(0.78, 0.9, r'$b=10$' ,fontsize=30)  
plt.savefig("fig_nonlin_MI.pdf", format='pdf', bbox_inches="tight")

plt.show()
