"""
test code: plot M-R
"""

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

def dimllatidal(utidal, ytidal):
    return 16./15. * (1.-2.*utidal)**2 *(2.-ytidal+2.*utidal*(ytidal-1)) / ( 8.*utidal**5*(ytidal+1) + 4.*utidal**4*(3.*ytidal-2.) + 4.*utidal**3*(13.-11.*ytidal) + 6.*utidal**2 *(5.*ytidal-8.) + 2.*utidal*(6.-3.*ytidal) +3.*(1.-2.*utidal)**2*(2.-ytidal+2.*utidal*(ytidal-1))*np.log(1.-2.*utidal) )
          
          
colorset = ['black', 'lightcoral','yellowgreen']

        
ytidalset = np.linspace(1.9, 2.3, 1000)
dimllatidalset1 = dimllatidal(0.01, ytidalset) 
dimllatidalset2 = dimllatidal(0.05, ytidalset) 
dimllatidalset3 = dimllatidal(0.1, ytidalset) 


#print(ytidalset)
#print(dimllatidalset1)
 
plt.rcParams['xtick.labelsize'] = 12  # number size
plt.rcParams['ytick.labelsize'] = 12
#plt.rcParams['xtick.major.size'] = 8 # little segment length
#plt.rcParams['ytick.major.size'] = 8        


plt.plot(ytidalset, dimllatidalset1, linestyle='solid', color=colorset[0], label='$C$ = 0.01')
plt.plot(ytidalset, dimllatidalset2, linestyle='dotted', color=colorset[1], label='$C$ = 0.05')
plt.plot(ytidalset, dimllatidalset3, linestyle='dashed', color=colorset[2], label='$C$ = 0.1')

plt.yscale('symlog')
plt.yticks([-1e8, -1e6, -1e4, -1e2, 0, 1e2, 1e4, 1e6, 1e8]) 
plt.minorticks_off()
plt.ylabel(r'$ \lambda/M^5$',fontsize=15)
plt.xlabel(r'$y$',fontsize=15)

plt.legend(fontsize = 12)



   
plt.savefig("fig_app4.pdf", format='pdf', bbox_inches="tight")
plt.show()          
