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

colorset=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']   


stdata1 = np.genfromtxt('stgb_linear_sol_data1.txt')
stdata2 = np.genfromtxt('stgb_linear_sol_data2.txt')
stdata3 = np.genfromtxt('stgb_linear_sol_data3.txt')
stdata4 = np.genfromtxt('stgb_linear_sol_data4.txt')
stdata5 = np.genfromtxt('stgb_linear_sol_data5.txt')
stdata6 = np.genfromtxt('stgb_linear_sol_data6.txt')

stn1 = len(stdata1)
#stn1, stn2, stn3, stn4, stn5 = len(strot1), len(strot2), len(strot3), len(strot4), len(strot5)
ntrim = 700

# ec, interior ph_r/ph, exterior ph_r/ph
d10, d11, d12 = stdata1[:, 0], stdata1[:, 1], stdata1[:, 2]
d20, d21, d22 = stdata2[:, 0], stdata2[:, 1], stdata2[:, 2]
d30, d31, d32 = stdata3[:, 0], stdata3[:, 1], stdata3[:, 2]
d40, d41, d42 = stdata4[:, 0], stdata4[:, 1], stdata4[:, 2]
d50, d51, d52 = stdata5[:, 0], stdata5[:, 1], stdata5[:, 2]
d60, d61, d62 = stdata6[:, 0], stdata6[:, 1], stdata6[:, 2]


fig = plt.figure(figsize=(6,12))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)




ax1.set_xlabel(r'$\epsilon_0$', fontsize=20)
ax1.set_ylabel(r'$\frac{\Phi^{\prime}}{\Phi}|_s$', fontsize=20)
ax2.set_xlabel(r'$\epsilon_0$', fontsize=20)
ax2.set_ylabel(r'$\frac{\Phi^{\prime}}{\Phi}|_s$', fontsize=20)


 
# ph_r/ph
ax1.set_ylim([-10, 10])
#ax1.set_xlim([0, 10])
ax1.plot( d40[0:ntrim], d41[0:ntrim] , '--', color = colorset[0])
ax1.plot( d40[0:ntrim], d42[0:ntrim] , color = colorset[0])
ax1.plot( d50[0:ntrim], d51[0:ntrim] , '--', color = colorset[1])
ax1.plot( d50[0:ntrim], d52[0:ntrim] , color = colorset[1])
ax1.plot( d60[0:ntrim], d61[0:ntrim] , '--', color = colorset[2])
ax1.plot( d60[0:ntrim], d62[0:ntrim] , color = colorset[2])




ax2.set_ylim([-10, 10])
#ax1.set_xlim([0, 10])
ax2.plot( d10[0:ntrim], d11[0:ntrim] , '--', color = colorset[0])
ax2.plot( d10[0:ntrim], d12[0:ntrim] , color = colorset[0])
ax2.plot( d20[0:ntrim], d21[0:ntrim] , '--', color = colorset[1])
ax2.plot( d20[0:ntrim], d22[0:ntrim] , color = colorset[1])
ax2.plot( d30[0:ntrim], d31[0:ntrim] , '--', color = colorset[2])
ax2.plot( d30[0:ntrim], d32[0:ntrim] , color = colorset[2])




plt.savefig("fig_lin1.pdf", format='pdf', bbox_inches="tight")
print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show()
