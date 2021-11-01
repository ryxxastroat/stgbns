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

plt.tick_params(axis='both', which='minor', labelsize=18)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
import warnings
import timeit
import numpy as np
import scipy.optimize
from numpy import pi
import matplotlib.pyplot as plt
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
plt.subplots_adjust(wspace=0.36)
# plt.setp(axs, xticks=[0, 0.5, 1,1.5,2,2.5])
# from matplotlib import ticker
# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True) 
# formatter.set_powerlimits((-5,-3)) 
# axs[0,2].yaxis.set_major_formatter(formatter)
# axs[0,2].set_yticks([2.5*10**(-4)])



for i in range(2):
    for j in range(4):
        axs[i,j].minorticks_on()
        
# the first plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+21)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
    axs[0,0].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 

list2=[26,27,28,30]
# the second plot
for i in range(4):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(list2[i])+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7]
    if list2[i]!=30:
        axs[0,1].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    else:
        axs[0,1].plot(b321/MSUN, b121, color = colorset[i+1],linewidth=2.5) 
        
# the third plot
for i in range(5):
    """
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+31)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
    axs[0,2].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    """
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+31)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
#     axs[1,2].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    index=b321.argsort()
    ydata=b121[index]
    xdata=b321[index]/MSUN
    s1 = UnivariateSpline(xdata, ydata, s=5)
    xs=np.linspace(min(xdata),max(xdata),20)
    ys=s1(xs)
    axs[0,2].plot(xs, ys, color = colorset[i],linewidth=2.5)    

# the fourth plot
for i in range(5):
    """
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+36)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
    axs[0,3].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    """
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+36)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
#     axs[1,2].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    index=b321.argsort()
    ydata=b121[index]
    xdata=b321[index]/MSUN
    s1 = UnivariateSpline(xdata, ydata, s=5)
    xs=np.linspace(min(xdata),max(xdata),20)
    ys=s1(xs)
    axs[0,3].plot(xs, ys, color = colorset[i],linewidth=2.5)
    
# the fifth plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+1)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
    axs[1,0].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    
    
# the sixth plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(6+i)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7]
    if i!=4:
        axs[1,1].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    else:
        axs[1,1].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5)
        
# the seventh plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+11)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
#     axs[1,2].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    index=b321.argsort()
    ydata=b121[index]
    xdata=b321[index]/MSUN
    s1 = UnivariateSpline(xdata, ydata, s=5)
    xs=np.linspace(min(xdata),max(xdata),20)
    ys=s1(xs)
    axs[1,2].plot(xs, ys, color = colorset[i],linewidth=2.5) 
    
# the eighth plot
for i in range(5):
    data=np.genfromtxt('stgb_tid_v1_comb_data'+str(i+16)+'.txt')
    b021, b121, b221, b321, b421, b521, b621, b721=data[:, 0], data[:, 1], \
    data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7] 
#     axs[1,3].plot(b321/MSUN, b121, color = colorset[i],linewidth=2.5) 
    index=b321.argsort()
    ydata=b121[index]
    xdata=b321[index]/MSUN
    s1 = UnivariateSpline(xdata, ydata, s=5)
    xs=np.linspace(min(xdata),max(xdata),20)
    ys=s1(xs)
    axs[1,3].plot(xs, ys, color = colorset[i],linewidth=2.5)
    
fig.text(0.08, 0.50,r'$\Phi_{\rm c}$'  , ha='center', fontsize=30,rotation='vertical')
fig.text(0.48, 0.05, r'$M\,[{ \rm M_{\odot}}]$',fontsize=30)      
fig.text(0.91, 0.7, r'$a=1$' ,fontsize=30, rotation='-90')          
fig.text(0.91, 0.27, r'$a=0.1$' ,fontsize=30, rotation='-90')     
fig.text(0.18, 0.9, r'$b=-1$' ,fontsize=30)     
fig.text(0.38, 0.9, r'$b=-0.2$' ,fontsize=30)     
fig.text(0.6, 0.9, r'$b=3$' ,fontsize=30)   
fig.text(0.8, 0.9, r'$b=10$' ,fontsize=30)  
plt.savefig("fig_nonlin_Mbphc.pdf", format='pdf', bbox_inches="tight")
plt.show()
