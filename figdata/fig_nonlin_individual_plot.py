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
#plt.tick_params(axis='both', which='minor', labelsize=18)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
                               
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

def schwarzgtt(r, mass):
   return -1 + 2.0*G*mass/r/c**2

   
colorset=['black', 'lightcoral', 'yellowgreen'] 
linewidthback = 2
linewidthfront = 3
  
grdata1 = np.genfromtxt('TOV_5eqs_sol_data14.txt')
grdata2 = np.genfromtxt('TOV_5eqs_sol_data20.txt')
grdata3 = np.genfromtxt('TOV_5eqs_sol_data22.txt')


# x, gtt, grr, p, e, omega
a10, a11, a12, a13, a14, a15 = grdata1[:, 0], grdata1[:, 1], grdata1[:, 2], grdata1[:, 3], grdata1[:, 4], grdata1[:, 5]
a20, a21, a22, a23, a24, a25 = grdata2[:, 0], grdata2[:, 1], grdata2[:, 2], grdata2[:, 3], grdata2[:, 4], grdata2[:, 5]
a30, a31, a32, a33, a34, a35 = grdata3[:, 0], grdata3[:, 1], grdata3[:, 2], grdata3[:, 3], grdata3[:, 4], grdata3[:, 5]
nptgrext = 1000
x1 = np.linspace(a10[len(a10)-1], 50*KM, nptgrext)
x2 = np.linspace(a20[len(a20)-1], 50*KM, nptgrext)
x3 = np.linspace(a30[len(a30)-1], 50*KM, nptgrext)
y1, y2, y3 = schwarzgtt(x1, 1.4*MSUN), schwarzgtt(x2, 2.0*MSUN), schwarzgtt(x3, 2.2*MSUN)
#print(a10)
#print(x1)


stdata1 = np.genfromtxt('stgb_solver_sol_data1.txt')
stdata2 = np.genfromtxt('stgb_solver_sol_data2.txt')
stdata3 = np.genfromtxt('stgb_solver_sol_data3.txt')
stdata4 = np.genfromtxt('stgb_solver_sol_data4.txt')
stdata5 = np.genfromtxt('stgb_solver_sol_data5.txt')


# r, nu, mu, phi, dphi/dr, p
c10, c11, c12, c13 = stdata1[:, 0], stdata1[:, 1], stdata1[:, 2], stdata1[:, 3]
gtt1, grr1 = np.zeros_like(c11), np.zeros_like(c12)
for i in range(0, len(c10)):
   gtt1[i] = -np.exp(2*c11[i])
   grr1[i] = np.exp(2*c12[i]) 
   
   
c20, c21, c22, c23 = stdata2[:, 0], stdata2[:, 1], stdata2[:, 2], stdata2[:, 3]
gtt2, grr2 = np.zeros_like(c20), np.zeros_like(c20)
for i in range(0, len(c20)):
   gtt2[i] = -np.exp(2*c21[i])
   grr2[i] = np.exp(2*c22[i])
   
   
c30, c31, c32, c33 = stdata3[:, 0], stdata3[:, 1], stdata3[:, 2], stdata3[:, 3]
gtt3, grr3 = np.zeros_like(c30), np.zeros_like(c30)
for i in range(0, len(c30)):
   gtt3[i] = -np.exp(2*c31[i])
   grr3[i] = np.exp(2*c32[i])


c40, c41, c42, c43 = stdata4[:, 0], stdata4[:, 1], stdata4[:, 2], stdata4[:, 3]
gtt4, grr4 = np.zeros_like(c40), np.zeros_like(c40)
for i in range(0, len(c40)):
   gtt4[i] = -np.exp(2*c41[i])
   grr4[i] = np.exp(2*c42[i])
   
   
c50, c51, c52, c53 = stdata5[:, 0], stdata5[:, 1], stdata5[:, 2], stdata5[:, 3]
gtt5, grr5 = np.zeros_like(c50), np.zeros_like(c50)
for i in range(0, len(c50)):
   gtt5[i] = -np.exp(2*c51[i])
   grr5[i] = np.exp(2*c52[i])
   
   
   
     




fig = plt.figure(figsize=(24,24))
ax1 = fig.add_subplot(331)
ax2 = fig.add_subplot(332)
ax3 = fig.add_subplot(333)
ax4 = fig.add_subplot(334)
ax5 = fig.add_subplot(335)
ax6 = fig.add_subplot(336)
ax7 = fig.add_subplot(337)
ax8 = fig.add_subplot(338)
ax9 = fig.add_subplot(339)

fontsizevalue1 = 30
ax1.set_ylabel(r'$g_{tt}$', fontsize=fontsizevalue1)
ax4.set_ylabel(r'$g_{rr}$', fontsize=fontsizevalue1)
ax7.set_ylabel(r'$\Phi$', fontsize=fontsizevalue1)

fig.text(0.2, 0.9, r'$M=1.4M_{\odot}$' ,fontsize=30)     
fig.text(0.47, 0.9, r'$M=2.0M_{\odot}$' ,fontsize=30)     
fig.text(0.74, 0.9, r'$M=2.2M_{\odot}$' ,fontsize=30)  
fig.text(0.48, 0.07, r'$r\,[\rm km]$' ,fontsize=30)



ax1.set_ylim([-1, -0.3])
ax2.set_ylim([-1, -0.1])
ax3.set_ylim([-1, -0.])
"""
ax4.set_xlim([0, 100])
ax5.set_xlim([0, 100])
ax5.set_ylim([1, 3])
ax6.set_xlim([0, 100])
ax6.set_ylim([-0.005, 0.018])

ax1.plot(a10/KM, a11)
ax2.plot(a10/KM, a12)
"""

ntrim1, ntrim2, ntrim3, ntrim4, ntrim5 = 5000, 5000, 5000, 5000, 5000



ax1.plot(c10[:ntrim1]*runit/KM, gtt1[:ntrim1], color = colorset[1], linewidth = linewidthfront, label='$a, b = 0.1, -1$')
ax4.plot(c10[:ntrim1]*runit/KM, grr1[:ntrim1], color = colorset[1], linewidth = linewidthfront)
ax7.plot(c10[:ntrim1]*runit/KM, c13[:ntrim1], color = colorset[1], linewidth = linewidthfront)
ax1.plot(c20[:ntrim2]*runit/KM, gtt2[:ntrim2], color = colorset[2], linewidth = linewidthfront, label='$a, b = 0.1, 3$')
ax4.plot(c20[:ntrim2]*runit/KM, grr2[:ntrim2], color = colorset[2], linewidth = linewidthfront)
ax7.plot(c20[:ntrim2]*runit/KM, c23[:ntrim2], color = colorset[2], linewidth = linewidthfront)
ax1.plot(a10/KM, a11, '--', color = colorset[0], linewidth = linewidthback)
ax4.plot(a10/KM, a12, '--', color = colorset[0], linewidth = linewidthback)
ax1.plot(x1/KM, y1, '--', color = colorset[0], linewidth = linewidthback, label='$a, b = 0, 0$ (GR)')
ax4.plot(x1/KM, -1./y1, '--', color = colorset[0], linewidth = linewidthback)


ax2.plot(c30[:ntrim3]*runit/KM, gtt3[:ntrim3], color = colorset[1], linewidth = linewidthfront, label='$a, b = 1, 3$')
ax5.plot(c30[:ntrim3]*runit/KM, grr3[:ntrim3], color = colorset[1], linewidth = linewidthfront)
ax8.plot(c30[:ntrim3]*runit/KM, c33[:ntrim3], color = colorset[1], linewidth = linewidthfront)
ax2.plot(a20/KM, a21, '--', color = colorset[0], linewidth = linewidthfront, label='$a, b = 0, 0$ (GR)')
ax5.plot(a20/KM, a22, '--', color = colorset[0], linewidth = linewidthfront)
ax2.plot(x2/KM, y2, '--', color = colorset[0], linewidth = linewidthback)
ax5.plot(x2/KM, -1./y2, '--', color = colorset[0], linewidth = linewidthback)


ax3.plot(c40[:ntrim4]*runit/KM, gtt4[:ntrim4], color = colorset[1], linewidth = linewidthfront, label='$a, b = 1, -0.2$')
ax6.plot(c40[:ntrim4]*runit/KM, grr4[:ntrim4], color = colorset[1], linewidth = linewidthfront)
ax9.plot(c40[:ntrim4]*runit/KM, c43[:ntrim4], color = colorset[1], linewidth = linewidthfront)
ax3.plot(c50[:ntrim5]*runit/KM, gtt5[:ntrim5], color = colorset[2], linewidth = linewidthfront, label='$a, b = 1, -1$')
ax6.plot(c50[:ntrim5]*runit/KM, grr5[:ntrim5], color = colorset[2], linewidth = linewidthfront)
ax9.plot(c50[:ntrim5]*runit/KM, c53[:ntrim5], color = colorset[2], linewidth = linewidthfront)
ax3.plot(a30/KM, a31, linestyle='--', color = colorset[0], linewidth = linewidthfront, label='$a, b = 0, 0$ (GR)')
ax6.plot(a30/KM, a32, '--', color = colorset[0], linewidth = linewidthfront)
ax3.plot(x3/KM, y3, '--', color = colorset[0], linewidth = linewidthback)
ax6.plot(x3/KM, -1./y3, '--', color = colorset[0], linewidth = linewidthback)


ax1.legend(fontsize=25, frameon=False)
ax2.legend(fontsize=25, frameon=False)
ax3.legend(fontsize=25, frameon=False)


plt.savefig("fig_nonlin_ind.pdf", format='pdf', bbox_inches="tight")

print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show()
