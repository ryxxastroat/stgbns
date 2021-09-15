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


stdata1 = np.genfromtxt('stgb_linear_v1_xia_data1.txt')
stdata2 = np.genfromtxt('stgb_linear_v1_xia_data2.txt')
stdata3 = np.genfromtxt('stgb_linear_v1_xia_data3.txt')
stdata4 = np.genfromtxt('stgb_linear_v1_xia_data4.txt')

stn1 = len(stdata1)
#stn1, stn2, stn3, stn4, stn5 = len(strot1), len(strot2), len(strot3), len(strot4), len(strot5)


# ec, interior ph_r/ph, exterior ph_r/ph
d10, d11, d12, d13, d14, d15 = stdata1[:, 0], stdata1[:, 1], stdata1[:, 2], stdata1[:, 3], stdata1[:, 4], stdata1[:, 5]
d20, d21, d22, d23, d24, d25 = stdata2[:, 0], stdata2[:, 1], stdata2[:, 2], stdata2[:, 3], stdata2[:, 4], stdata2[:, 5]
d30, d31, d32, d33, d34, d35 = stdata3[:, 0], stdata3[:, 1], stdata3[:, 2], stdata3[:, 3], stdata3[:, 4], stdata3[:, 5]
d40, d41, d42, d43, d44, d45 = stdata4[:, 0], stdata4[:, 1], stdata4[:, 2], stdata4[:, 3], stdata4[:, 4], stdata4[:, 5]

fig = plt.figure(figsize=(6,12))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)





ax1.set_xlabel(r'$a$', fontsize=20)
ax1.set_ylabel(r'$\xi$', fontsize=20)
ax2.set_xlabel(r'$a$', fontsize=20)
ax2.set_ylabel(r'$\xi$', fontsize=20)

 
# ph_r/ph
#ax1.set_ylim([-10, 10])
ax2.plot( d10, d11 , color = colorset[0])
ax2.plot( d10, d12 , color = colorset[1])
ax2.plot( d10, d13 , color = colorset[2])
ax2.plot( d10, d14 , color = colorset[3])
ax2.plot( d10, d15 , color = colorset[4])

ax2.plot( d20, d21 , '--', color = colorset[0])
ax2.plot( d20, d22 , '--', color = colorset[1])
ax2.plot( d20, d23 , '--', color = colorset[2])
ax2.plot( d20, d24 , '--', color = colorset[3])
ax2.plot( d20, d25 , '--', color = colorset[4])

ax1.plot( d30, d31 , color = colorset[0])
ax1.plot( d30, d32 , color = colorset[1])
ax1.plot( d30, d33 , color = colorset[2])
ax1.plot( d30, d34 , color = colorset[3])
ax1.plot( d30, d35 , color = colorset[4])

ax1.plot( d40, d41 , '--', color = colorset[0])
ax1.plot( d40, d42 , '--', color = colorset[1])
ax1.plot( d40, d43 , '--', color = colorset[2])
ax1.plot( d40, d44 , '--', color = colorset[3])
ax1.plot( d40, d45 , '--', color = colorset[4])


"""
# M vs -nu_s, mu1_s
ax1.plot( c11/MSUN, -0.5*np.log( (1.-2.*mdiml(c11)/rdiml(c12)) ) , color = colorset[0])
ax1.plot( c21/MSUN, -0.5*np.log( (1.-2.*mdiml(c21)/rdiml(c22)) ) , color = colorset[1])
ax1.plot( c31/MSUN, -0.5*np.log( (1.-2.*mdiml(c31)/rdiml(c32)) ) , color = colorset[2])
ax1.plot( c41/MSUN, -0.5*np.log( (1.-2.*mdiml(c41)/rdiml(c42)) ) , color = colorset[3])
ax1.plot( c51/MSUN, -0.5*np.log( (1.-2.*mdiml(c51)/rdiml(c52)) ) , color = colorset[4])
ax1.plot(d13/MSUN, -d14, 'o')
ax1.plot(d23/MSUN, -d24, 'o')
ax1.plot(d33/MSUN, -d34, 'o')
ax1.plot(d43/MSUN, -d44, 'o')
ax1.plot(d53/MSUN, -d54, 'o')
ax1.plot(d13/MSUN, 0.5*np.log(1./(1.-2.*mdiml(d15)/rdiml(d12))), 'v', color = colorset[0])
ax1.plot(d23/MSUN, 0.5*np.log(1./(1.-2.*mdiml(d25)/rdiml(d22))), 'v', color = colorset[1])
ax1.plot(d33/MSUN, 0.5*np.log(1./(1.-2.*mdiml(d35)/rdiml(d32))), 'v', color = colorset[2])
ax1.plot(d43/MSUN, 0.5*np.log(1./(1.-2.*mdiml(d45)/rdiml(d42))), 'v', color = colorset[3])
ax1.plot(d53/MSUN, 0.5*np.log(1./(1.-2.*mdiml(d55)/rdiml(d52))), 'v', color = colorset[4])

"""





plt.savefig("fig_lin2.pdf", format='pdf', bbox_inches="tight")
print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show()
