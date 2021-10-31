#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
import numpy as np
import scipy.optimize
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import ode as sp_ode

t0 = timeit.time.time()
from shapely.geometry import LineString


plt.close()


# In[44]:


PI4 = 4.0 * pi
ka = 8.0 * pi
c = 29979245800.0  # cm/s
G = 6.67408e-8  # cm^3/g/s^2

MSUN = 1.98855e33  # g
KM = 1.0e5  # cm
mB = 1.660538921e-24  # g
E_NUCL = 2.0e14  # minimun energy density for NS core; g/cm^3

runit = 10.*KM # length to parametrize quantities

def rdiml(r):
    """ dimensionless mass """
    return r/runit    
def mdiml(m):
    """ dimensionless mass """
    return G*m/c**2/runit
def pdiml(p):
    """ dimensionless pressure """
    return G*p/c**4 * runit**2




# In[42]:


colorset=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']   

fig, (ax1,ax2) = plt.subplots(2, 1,figsize=(10,12),sharex=True)
plt.subplots_adjust(hspace=0.0)
from scipy.interpolate import UnivariateSpline

# ax1.set_xlabel(r'$a$', fontsize=30)
ax1.set_ylabel(r'$b_{+}$', fontsize=30)
ax2.set_xlabel(r'$a$', fontsize=30)
ax2.set_ylabel(r'$b_{-}$', fontsize=30)

 
# ph_r/ph
#ax1.set_ylim([-10, 10])

label=['WFF1','SLy4','AP4', 'MPA1','PAL1']
data1 = np.genfromtxt('stgb_linear_v1_xia_data1.txt')
for j in range(5):
    x=data1[:,0]
    y=data1[:,j+1]
    s1 = UnivariateSpline(x, y, s=5)
    xs=np.linspace(min(x),max(x),20)
    ys=s1(xs)
    ax2.plot( xs,ys , color = colorset[j],linewidth=1.5)

data2 = np.genfromtxt('stgb_linear_v1_xia_data2.txt')
for j in range(5):
    x=data2[:,0]
    y=data2[:,j+1]
    s1 = UnivariateSpline(x, y, s=5)
    xs=np.linspace(min(x),max(x),20)
    ys=s1(xs)
    ax2.plot( xs,ys , color = colorset[j],linewidth=1.5,linestyle='--')
    
    
data3 = np.genfromtxt('stgb_linear_v1_xia_data3.txt')
for j in range(5):
    x=data3[:,0]
    y=data3[:,j+1]
    s1 = UnivariateSpline(x, y, s=5)
    xs=np.linspace(min(x),max(x),20)
    ys=s1(xs)
    ax1.plot( xs,ys , color = colorset[j],linewidth=1.5,label=label[j])
    
data4 = np.genfromtxt('stgb_linear_v1_xia_data4.txt')
for j in range(5):
    x=data4[:,0]
    y=data4[:,j+1]
    s1 = UnivariateSpline(x, y, s=5)
    xs=np.linspace(min(x),max(x),20)
    ys=s1(xs)
    ax1.plot( xs,ys , color = colorset[j],linewidth=1.5,linestyle='--')

ax2.set_ylim(-0.65,0.1)
ax1.set_ylim(-0.5,10.5)
plt.xlim(0,3)

ax1.minorticks_on()
ax2.minorticks_on()
ax1.legend(fontsize=20,frameon=False)


plt.savefig("fig_lin2.pdf", format='pdf', bbox_inches="tight")
# print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show()


