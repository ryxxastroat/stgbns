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
  
grdata1 = np.genfromtxt('TOV_5eqs_data1.txt')
grdata2 = np.genfromtxt('TOV_5eqs_data2.txt')
grdata3 = np.genfromtxt('TOV_5eqs_data3.txt')
grdata4 = np.genfromtxt('TOV_5eqs_data4.txt')
grdata5 = np.genfromtxt('TOV_5eqs_data5.txt')
# ec, M, R, Mbar, I, J, om_inf, om_s
a10, a11, a12, a13, a14, a15, a16, a17 = grdata1[:, 0], grdata1[:, 1], grdata1[:, 2], grdata1[:, 3], grdata1[:, 4], grdata1[:, 5], grdata1[:, 6], grdata1[:, 7]
a20, a21, a22, a23, a24, a25, a26, a27 = grdata2[:, 0], grdata2[:, 1], grdata2[:, 2], grdata2[:, 3], grdata2[:, 4], grdata2[:, 5], grdata2[:, 6], grdata2[:, 7]
a30, a31, a32, a33, a34, a35, a36, a37 = grdata3[:, 0], grdata3[:, 1], grdata3[:, 2], grdata3[:, 3], grdata3[:, 4], grdata3[:, 5], grdata3[:, 6], grdata3[:, 7]
a40, a41, a42, a43, a44, a45, a46, a47 = grdata4[:, 0], grdata4[:, 1], grdata4[:, 2], grdata4[:, 3], grdata4[:, 4], grdata4[:, 5], grdata4[:, 6], grdata4[:, 7]
a50, a51, a52, a53, a54, a55, a56, a57 = grdata5[:, 0], grdata5[:, 1], grdata5[:, 2], grdata5[:, 3], grdata5[:, 4], grdata5[:, 5], grdata5[:, 6], grdata5[:, 7]



fig = plt.figure(figsize=(24,12))
ax1 = fig.add_subplot(241)
ax1 = fig.add_subplot(241)
ax2 = fig.add_subplot(242)
ax3 = fig.add_subplot(243)
ax4 = fig.add_subplot(244)
ax5 = fig.add_subplot(245)
ax6 = fig.add_subplot(246)
ax7 = fig.add_subplot(247)
ax8 = fig.add_subplot(248)




"""
ax1.set_xlabel(r'$R$', fontsize=20)
ax1.set_ylabel(r'$M$', fontsize=20)
ax2.set_xlabel(r'$x$', fontsize=20)
ax2.set_ylabel(r'$\mu$', fontsize=20)
ax3.set_xlabel(r'$x$', fontsize=20)
ax3.set_ylabel(r'$\Phi$', fontsize=20)
ax4.set_xlabel(r'$x$', fontsize=20)
ax4.set_ylabel(r'$\Phi^{\prime}$', fontsize=20)  

ax5.set_xlabel(r'$R$', fontsize=20)
ax5.set_ylabel(r'$M$', fontsize=20)

ax8.set_xlim([0.5e+15, 1.5e+15])
ax8.set_xlabel(r'$ec$', fontsize=20)
ax8.set_ylabel(r'$\Omega$', fontsize=20)
"""
 
# GR M-R
Rmin, Rmax = 8, 15
#Mmin, Mmax = 0.2, 2.5
ax1.set_xlim([Rmin, Rmax])
#ax1.set_ylim([Mmin, Mmax])
ax1.plot(a12/KM, a11/MSUN, '--', color = colorset[0], label='$a$ = 1')
ax1.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax1.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax1.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax1.plot(a52/KM, a51/MSUN, '--', color = colorset[4])
ax2.set_xlim([Rmin, Rmax])
ax2.plot(a12/KM, a11/MSUN, '--', color = colorset[0])
ax2.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax2.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax2.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax2.plot(a52/KM, a51/MSUN, '--', color = colorset[4])
ax3.set_xlim([Rmin, Rmax])
ax3.plot(a12/KM, a11/MSUN, '--', color = colorset[0])
ax3.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax3.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax3.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax3.plot(a52/KM, a51/MSUN, '--', color = colorset[4])
ax4.set_xlim([Rmin, Rmax])
ax4.plot(a12/KM, a11/MSUN, '--', color = colorset[0])
ax4.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax4.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax4.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax4.plot(a52/KM, a51/MSUN, '--', color = colorset[4])
ax5.set_xlim([Rmin, Rmax])
ax5.plot(a12/KM, a11/MSUN, '--', color = colorset[0], label='$a$ = 0.1')
ax5.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax5.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax5.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax5.plot(a52/KM, a51/MSUN, '--', color = colorset[4])
ax6.set_xlim([Rmin, Rmax])
ax6.plot(a12/KM, a11/MSUN, '--', color = colorset[0])
ax6.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax6.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax6.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax6.plot(a52/KM, a51/MSUN, '--', color = colorset[4])
ax7.set_xlim([Rmin, Rmax])
ax7.plot(a12/KM, a11/MSUN, '--', color = colorset[0])
ax7.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax7.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax7.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax7.plot(a52/KM, a51/MSUN, '--', color = colorset[4])
ax8.set_xlim([Rmin, Rmax])
ax8.plot(a12/KM, a11/MSUN, '--', color = colorset[0])
ax8.plot(a22/KM, a21/MSUN, '--', color = colorset[1])
ax8.plot(a32/KM, a31/MSUN, '--', color = colorset[2])
ax8.plot(a42/KM, a41/MSUN, '--', color = colorset[3])
ax8.plot(a52/KM, a51/MSUN, '--', color = colorset[4])

ax1.legend()
ax5.legend()







# ST M-R
# solver: ec, bphc, R, M_inf, nu_s, M_s, bph_s, psi_s, dimp_s, r_inf
# rot: ec, bphc, R, de R, M_inf, de M_inf,  Mbar, I, J, om_inf
stdata1= np.genfromtxt('stgb_solver_pdata1.txt') 
b01, b11, b21, b31, b41, b51, b61, b71, b81, b91=stdata1[:, 0], stdata1[:, 1], stdata1[:, 2], stdata1[:, 3], stdata1[:, 4], stdata1[:, 5], stdata1[:, 6], stdata1[:, 7], stdata1[:, 8], stdata1[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b01)-1):
   testecsep = b01[ii+1]-b01[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b01))
for ii in range(0, len(ntrimset)-1):
   ax5.plot(b31[ntrimset[ii]:ntrimset[ii+1]]/KM, b41[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[0]) 
   ax5.plot(b31[ntrimset[ii]:ntrimset[ii+1]]/KM, b41[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[0]) 
stdata2= np.genfromtxt('stgb_solver_pdata2.txt') 
b02, b12, b22, b32, b42, b52, b62, b72, b82, b92=stdata2[:, 0], stdata2[:, 1], stdata2[:, 2], stdata2[:, 3], stdata2[:, 4], stdata2[:, 5], stdata2[:, 6], stdata2[:, 7], stdata2[:, 8], stdata2[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b02)-1):
   testecsep = b02[ii+1]-b02[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b02))
for ii in range(0, len(ntrimset)-1):
   ax5.plot(b32[ntrimset[ii]:ntrimset[ii+1]]/KM, b42[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[1]) 
   ax5.plot(b32[ntrimset[ii]:ntrimset[ii+1]]/KM, b42[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[1]) 
stdata3= np.genfromtxt('stgb_solver_pdata3.txt') 
b03, b13, b23, b33, b43, b53, b63, b73, b83, b93=stdata3[:, 0], stdata3[:, 1], stdata3[:, 2], stdata3[:, 3], stdata3[:, 4], stdata3[:, 5], stdata3[:, 6], stdata3[:, 7], stdata3[:, 8], stdata3[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b03)-1):
   testecsep = b03[ii+1]-b03[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b03))
for ii in range(0, len(ntrimset)-1):
   ax5.plot(b33[ntrimset[ii]:ntrimset[ii+1]]/KM, b43[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[2]) 
   ax5.plot(b33[ntrimset[ii]:ntrimset[ii+1]]/KM, b43[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[2])
stdata4= np.genfromtxt('stgb_solver_pdata4.txt') 
b04, b14, b24, b34, b44, b54, b64, b74, b84, b94=stdata4[:, 0], stdata4[:, 1], stdata4[:, 2], stdata4[:, 3], stdata4[:, 4], stdata4[:, 5], stdata4[:, 6], stdata4[:, 7], stdata4[:, 8], stdata4[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b04)-1):
   testecsep = b04[ii+1]-b04[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b04))
for ii in range(0, len(ntrimset)-1):
   ax5.plot(b34[ntrimset[ii]:ntrimset[ii+1]]/KM, b44[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[3]) 
   ax5.plot(b34[ntrimset[ii]:ntrimset[ii+1]]/KM, b44[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[3])
stdata5= np.genfromtxt('stgb_solver_pdata5.txt') 
b05, b15, b25, b35, b45, b55, b65, b75, b85, b95=stdata5[:, 0], stdata5[:, 1], stdata5[:, 2], stdata5[:, 3], stdata5[:, 4], stdata5[:, 5], stdata5[:, 6], stdata5[:, 7], stdata5[:, 8], stdata5[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b05)-1):
   testecsep = b05[ii+1]-b05[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b05))
for ii in range(0, len(ntrimset)-1):
   ax5.plot(b35[ntrimset[ii]:ntrimset[ii+1]]/KM, b45[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[4]) 
   ax5.plot(b35[ntrimset[ii]:ntrimset[ii+1]]/KM, b45[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[4])
   

stdata6= np.genfromtxt('stgb_solver_pdata6.txt') 
b06, b16, b26, b36, b46, b56, b66, b76, b86, b96=stdata6[:, 0], stdata6[:, 1], stdata6[:, 2], stdata6[:, 3], stdata6[:, 4], stdata6[:, 5], stdata6[:, 6], stdata6[:, 7], stdata6[:, 8], stdata6[:, 9]
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b06)-1):
   testecsep = b06[ii+1]-b06[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b06))
for ii in range(0, len(ntrimset)-1):
   ax6.plot(b36[ntrimset[ii]:ntrimset[ii+1]]/KM, b46[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[0]) 
   ax6.plot(b36[ntrimset[ii]:ntrimset[ii+1]]/KM, b46[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[0])  
stdata7= np.genfromtxt('stgb_solver_pdata7.txt') 
b07, b17, b27, b37, b47, b57, b67, b77, b87, b97=stdata7[:, 0], stdata7[:, 1], stdata7[:, 2], stdata7[:, 3], stdata7[:, 4], stdata7[:, 5], stdata7[:, 6], stdata7[:, 7], stdata7[:, 8], stdata7[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b07)-1):
   testecsep = b07[ii+1]-b07[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b07))
for ii in range(0, len(ntrimset)-1):
   ax6.plot(b37[ntrimset[ii]:ntrimset[ii+1]]/KM, b47[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[1]) 
   ax6.plot(b37[ntrimset[ii]:ntrimset[ii+1]]/KM, b47[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[1]) 
stdata8= np.genfromtxt('stgb_solver_pdata8.txt') 
b08, b18, b28, b38, b48, b58, b68, b78, b88, b98=stdata8[:, 0], stdata8[:, 1], stdata8[:, 2], stdata8[:, 3], stdata8[:, 4], stdata8[:, 5], stdata8[:, 6], stdata8[:, 7], stdata8[:, 8], stdata8[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b08)-1):
   testecsep = b08[ii+1]-b08[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b08))
for ii in range(0, len(ntrimset)-1):
   ax6.plot(b38[ntrimset[ii]:ntrimset[ii+1]]/KM, b48[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[2]) 
   ax6.plot(b38[ntrimset[ii]:ntrimset[ii+1]]/KM, b48[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[2]) 
stdata9= np.genfromtxt('stgb_solver_pdata9.txt') 
b09, b19, b29, b39, b49, b59, b69, b79, b89, b99=stdata9[:, 0], stdata9[:, 1], stdata9[:, 2], stdata9[:, 3], stdata9[:, 4], stdata9[:, 5], stdata9[:, 6], stdata9[:, 7], stdata9[:, 8], stdata9[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b09)-1):
   testecsep = b09[ii+1]-b09[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b09))
for ii in range(0, len(ntrimset)-1):
   ax6.plot(b39[ntrimset[ii]:ntrimset[ii+1]]/KM, b49[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[3]) 
   ax6.plot(b39[ntrimset[ii]:ntrimset[ii+1]]/KM, b49[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[3]) 
stdata10= np.genfromtxt('stgb_solver_pdata10.txt') 
b010, b110, b210, b310, b410, b510, b610, b710, b810, b910=stdata10[:, 0], stdata10[:, 1], stdata10[:, 2], stdata10[:, 3], stdata10[:, 4], stdata10[:, 5], stdata10[:, 6], stdata10[:, 7], stdata10[:, 8], stdata10[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b010)-1):
   testecsep = b010[ii+1]-b010[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b010))
for ii in range(0, len(ntrimset)-1):
   ax6.plot(b310[ntrimset[ii]:ntrimset[ii+1]]/KM, b410[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[4]) 
   ax6.plot(b310[ntrimset[ii]:ntrimset[ii+1]]/KM, b410[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[4]) 
   

stdata11= np.genfromtxt('stgb_solver_pdata11.txt') 
b011, b111, b211, b311, b411, b511, b611, b711, b811, b911=stdata11[:, 0], stdata11[:, 1], stdata11[:, 2], stdata11[:, 3], stdata11[:, 4], stdata11[:, 5], stdata11[:, 6], stdata11[:, 7], stdata11[:, 8], stdata11[:, 9] 
ax7.plot(b311/KM, b411/MSUN, color = colorset[0]) 
ax7.plot(b311/KM, b411/MSUN, 'o') 
stdata12= np.genfromtxt('stgb_solver_pdata12.txt') 
b012, b112, b212, b312, b412, b512, b612, b712, b812, b912=stdata12[:, 0], stdata12[:, 1], stdata12[:, 2], stdata12[:, 3], stdata12[:, 4], stdata12[:, 5], stdata12[:, 6], stdata12[:, 7], stdata12[:, 8], stdata12[:, 9] 
ax7.plot(b312/KM, b412/MSUN, color = colorset[1])
ax7.plot(b312/KM, b412/MSUN, 'o') 
stdata13= np.genfromtxt('stgb_solver_pdata13.txt') 
b013, b113, b213, b313, b413, b513, b613, b713, b813, b913=stdata13[:, 0], stdata13[:, 1], stdata13[:, 2], stdata13[:, 3], stdata13[:, 4], stdata13[:, 5], stdata13[:, 6], stdata13[:, 7], stdata13[:, 8], stdata13[:, 9] 
ax7.plot(b313/KM, b413/MSUN, color = colorset[2])
ax7.plot(b313/KM, b413/MSUN, 'o') 
stdata14= np.genfromtxt('stgb_solver_pdata14.txt') 
b014, b114, b214, b314, b414, b514, b614, b714, b814, b914=stdata14[:, 0], stdata14[:, 1], stdata14[:, 2], stdata14[:, 3], stdata14[:, 4], stdata14[:, 5], stdata14[:, 6], stdata14[:, 7], stdata14[:, 8], stdata14[:, 9] 
ax7.plot(b314/KM, b414/MSUN, color = colorset[3]) 
ax7.plot(b314/KM, b414/MSUN, 'o') 
stdata15= np.genfromtxt('stgb_solver_pdata15.txt') 
b015, b115, b215, b315, b415, b515, b615, b715, b815, b915=stdata15[:, 0], stdata15[:, 1], stdata15[:, 2], stdata15[:, 3], stdata15[:, 4], stdata15[:, 5], stdata15[:, 6], stdata15[:, 7], stdata15[:, 8], stdata15[:, 9] 
ax7.plot(b315/KM, b415/MSUN, color = colorset[4]) 
ax7.plot(b315/KM, b415/MSUN, 'o')

stdata16= np.genfromtxt('stgb_solver_pdata16.txt') 
b016, b116, b216, b316, b416, b516, b616, b716, b816, b916=stdata16[:, 0], stdata16[:, 1], stdata16[:, 2], stdata16[:, 3], stdata16[:, 4], stdata16[:, 5], stdata16[:, 6], stdata16[:, 7], stdata16[:, 8], stdata16[:, 9] 
ax8.plot(b316/KM, b416/MSUN, 'o', color = colorset[0]) 
stdata17= np.genfromtxt('stgb_solver_pdata17.txt') 
b017, b117, b217, b317, b417, b517, b617, b717, b817, b917=stdata17[:, 0], stdata17[:, 1], stdata17[:, 2], stdata17[:, 3], stdata17[:, 4], stdata17[:, 5], stdata17[:, 6], stdata17[:, 7], stdata17[:, 8], stdata17[:, 9] 
ax8.plot(b317/KM, b417/MSUN, 'o', color = colorset[1]) 
stdata18= np.genfromtxt('stgb_solver_pdata18.txt') 
b018, b118, b218, b318, b418, b518, b618, b718, b818, b918=stdata18[:, 0], stdata18[:, 1], stdata18[:, 2], stdata18[:, 3], stdata18[:, 4], stdata18[:, 5], stdata18[:, 6], stdata18[:, 7], stdata18[:, 8], stdata18[:, 9] 
ax8.plot(b318/KM, b418/MSUN, 'o', color = colorset[2]) 
stdata19= np.genfromtxt('stgb_solver_pdata19.txt') 
b019, b119, b219, b319, b419, b519, b619, b719, b819, b919=stdata19[:, 0], stdata19[:, 1], stdata19[:, 2], stdata19[:, 3], stdata19[:, 4], stdata19[:, 5], stdata19[:, 6], stdata19[:, 7], stdata19[:, 8], stdata19[:, 9] 
ax8.plot(b319/KM, b419/MSUN, 'o', color = colorset[3]) 
stdata20= np.genfromtxt('stgb_solver_pdata20.txt') 
b020, b120, b220, b320, b420, b520, b620, b720, b820, b920=stdata20[:, 0], stdata20[:, 1], stdata20[:, 2], stdata20[:, 3], stdata20[:, 4], stdata20[:, 5], stdata20[:, 6], stdata20[:, 7], stdata20[:, 8], stdata20[:, 9] 
ax8.plot(b320/KM, b420/MSUN, 'o', color = colorset[4]) 

stdata21= np.genfromtxt('stgb_solver_pdata21.txt') 
b021, b121, b221, b321, b421, b521, b621, b721, b821, b921=stdata21[:, 0], stdata21[:, 1], stdata21[:, 2], stdata21[:, 3], stdata21[:, 4], stdata21[:, 5], stdata21[:, 6], stdata21[:, 7], stdata21[:, 8], stdata21[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b021)-1):
   testecsep = b021[ii+1]-b021[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b021))
for ii in range(0, len(ntrimset)-1):
   ax1.plot(b321[ntrimset[ii]:ntrimset[ii+1]]/KM, b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[0]) 
   ax1.plot(b321[ntrimset[ii]:ntrimset[ii+1]]/KM, b421[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[0])  
stdata22= np.genfromtxt('stgb_solver_pdata22.txt') 
b022, b122, b222, b322, b422, b522, b622, b722, b822, b922=stdata22[:, 0], stdata22[:, 1], stdata22[:, 2], stdata22[:, 3], stdata22[:, 4], stdata22[:, 5], stdata22[:, 6], stdata22[:, 7], stdata22[:, 8], stdata22[:, 9]
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b022)-1):
   testecsep = b022[ii+1]-b022[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b022))
for ii in range(0, len(ntrimset)-1):
   ax1.plot(b322[ntrimset[ii]:ntrimset[ii+1]]/KM, b422[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[1]) 
   ax1.plot(b322[ntrimset[ii]:ntrimset[ii+1]]/KM, b422[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[1])  
stdata23= np.genfromtxt('stgb_solver_pdata23.txt') 
b023, b123, b223, b323, b423, b523, b623, b723, b823, b923=stdata23[:, 0], stdata23[:, 1], stdata23[:, 2], stdata23[:, 3], stdata23[:, 4], stdata23[:, 5], stdata23[:, 6], stdata23[:, 7], stdata23[:, 8], stdata23[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b023)-1):
   testecsep = b023[ii+1]-b023[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b023))
for ii in range(0, len(ntrimset)-1):
   ax1.plot(b323[ntrimset[ii]:ntrimset[ii+1]]/KM, b423[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[2]) 
   ax1.plot(b323[ntrimset[ii]:ntrimset[ii+1]]/KM, b423[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[2]) 
stdata24= np.genfromtxt('stgb_solver_pdata24.txt') 
b024, b124, b224, b324, b424, b524, b624, b724, b824, b924=stdata24[:, 0], stdata24[:, 1], stdata24[:, 2], stdata24[:, 3], stdata24[:, 4], stdata24[:, 5], stdata24[:, 6], stdata24[:, 7], stdata24[:, 8], stdata24[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b024)-1):
   testecsep = b024[ii+1]-b024[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b024))
for ii in range(0, len(ntrimset)-1):
   ax1.plot(b324[ntrimset[ii]:ntrimset[ii+1]]/KM, b424[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[3]) 
   ax1.plot(b324[ntrimset[ii]:ntrimset[ii+1]]/KM, b424[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[3]) 
stdata25= np.genfromtxt('stgb_solver_pdata25.txt') 
b025, b125, b225, b325, b425, b525, b625, b725, b825, b925=stdata25[:, 0], stdata25[:, 1], stdata25[:, 2], stdata25[:, 3], stdata25[:, 4], stdata25[:, 5], stdata25[:, 6], stdata25[:, 7], stdata25[:, 8], stdata25[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b025)-1):
   testecsep = b025[ii+1]-b025[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b025))
for ii in range(0, len(ntrimset)-1):
   ax1.plot(b325[ntrimset[ii]:ntrimset[ii+1]]/KM, b425[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[4]) 
   ax1.plot(b325[ntrimset[ii]:ntrimset[ii+1]]/KM, b425[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[4]) 


stdata26= np.genfromtxt('stgb_solver_pdata26.txt') 
b026, b126, b226, b326, b426, b526, b626, b726, b826, b926=stdata26[:, 0], stdata26[:, 1], stdata26[:, 2], stdata26[:, 3], stdata26[:, 4], stdata26[:, 5], stdata26[:, 6], stdata26[:, 7], stdata26[:, 8], stdata26[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b026)-1):
   testecsep = b026[ii+1]-b026[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b026))
for ii in range(0, len(ntrimset)-1):
   ax2.plot(b326[ntrimset[ii]:ntrimset[ii+1]]/KM, b426[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[0]) 
   ax2.plot(b326[ntrimset[ii]:ntrimset[ii+1]]/KM, b426[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[0])
stdata27= np.genfromtxt('stgb_solver_pdata27.txt') 
b027, b127, b227, b327, b427, b527, b627, b727, b827, b927=stdata27[:, 0], stdata27[:, 1], stdata27[:, 2], stdata27[:, 3], stdata27[:, 4], stdata27[:, 5], stdata27[:, 6], stdata27[:, 7], stdata27[:, 8], stdata27[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b027)-1):
   testecsep = b027[ii+1]-b027[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b027))
for ii in range(0, len(ntrimset)-1):
   ax2.plot(b327[ntrimset[ii]:ntrimset[ii+1]]/KM, b427[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[1]) 
   ax2.plot(b327[ntrimset[ii]:ntrimset[ii+1]]/KM, b427[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[1])
stdata28= np.genfromtxt('stgb_solver_pdata28.txt') 
b028, b128, b228, b328, b428, b528, b628, b728, b828, b928=stdata28[:, 0], stdata28[:, 1], stdata28[:, 2], stdata28[:, 3], stdata28[:, 4], stdata28[:, 5], stdata28[:, 6], stdata28[:, 7], stdata28[:, 8], stdata28[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b028)-1):
   testecsep = b028[ii+1]-b028[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b028))
for ii in range(0, len(ntrimset)-1):
   ax2.plot(b328[ntrimset[ii]:ntrimset[ii+1]]/KM, b428[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[2]) 
   ax2.plot(b328[ntrimset[ii]:ntrimset[ii+1]]/KM, b428[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[2]) 
#stdata29= np.genfromtxt('stgb_solver_pdata29.txt') 
#b029, b129, b229, b329, b429, b529, b629, b729, b829, b929=stdata29[:, 0], stdata29[:, 1], stdata29[:, 2], stdata29[:, 3], stdata29[:, 4], stdata29[:, 5], stdata29[:, 6], stdata29[:, 7], stdata29[:, 8], stdata29[:, 9] 
#ax2.plot(b329/KM, b429/MSUN, 'o') 
stdata30= np.genfromtxt('stgb_solver_pdata30.txt') 
b030, b130, b230, b330, b430, b530, b630, b730, b830, b930=stdata30[:, 0], stdata30[:, 1], stdata30[:, 2], stdata30[:, 3], stdata30[:, 4], stdata30[:, 5], stdata30[:, 6], stdata30[:, 7], stdata30[:, 8], stdata30[:, 9] 
ecstep = 0.6e+15
ntrimset = np.array([0])
for ii in range(0, len(b030)-1):
   testecsep = b030[ii+1]-b030[ii]
   if testecsep > ecstep:     
     #print(ii)
     ntrimset = np.append(ntrimset, ii+1)
ntrimset = np.append(ntrimset, len(b030))
for ii in range(0, len(ntrimset)-1):
   ax2.plot(b330[ntrimset[ii]:ntrimset[ii+1]]/KM, b430[ntrimset[ii]:ntrimset[ii+1]]/MSUN, color = colorset[4]) 
   ax2.plot(b330[ntrimset[ii]:ntrimset[ii+1]]/KM, b430[ntrimset[ii]:ntrimset[ii+1]]/MSUN, 'o', color = colorset[4])
 

stdata31= np.genfromtxt('stgb_solver_pdata31.txt') 
b031, b131, b231, b331, b431, b531, b631, b731, b831, b931=stdata31[:, 0], stdata31[:, 1], stdata31[:, 2], stdata31[:, 3], stdata31[:, 4], stdata31[:, 5], stdata31[:, 6], stdata31[:, 7], stdata31[:, 8], stdata31[:, 9] 
ax3.plot(b331/KM, b431/MSUN, 'o') 
stdata32= np.genfromtxt('stgb_solver_pdata32.txt') 
b032, b132, b232, b332, b432, b532, b632, b732, b832, b932=stdata32[:, 0], stdata32[:, 1], stdata32[:, 2], stdata32[:, 3], stdata32[:, 4], stdata32[:, 5], stdata32[:, 6], stdata32[:, 7], stdata32[:, 8], stdata32[:, 9] 
ax3.plot(b332/KM, b432/MSUN, 'o') 
stdata33= np.genfromtxt('stgb_solver_pdata33.txt') 
b033, b133, b233, b333, b433, b533, b633, b733, b833, b933=stdata33[:, 0], stdata33[:, 1], stdata33[:, 2], stdata33[:, 3], stdata33[:, 4], stdata33[:, 5], stdata33[:, 6], stdata33[:, 7], stdata33[:, 8], stdata33[:, 9] 
ax3.plot(b333/KM, b433/MSUN, 'o') 
stdata34= np.genfromtxt('stgb_solver_pdata34.txt') 
b034, b134, b234, b334, b434, b534, b634, b734, b834, b934=stdata34[:, 0], stdata34[:, 1], stdata34[:, 2], stdata34[:, 3], stdata34[:, 4], stdata34[:, 5], stdata34[:, 6], stdata34[:, 7], stdata34[:, 8], stdata34[:, 9] 
ax3.plot(b334/KM, b434/MSUN, 'o') 
stdata35= np.genfromtxt('stgb_solver_pdata35.txt') 
b035, b135, b235, b335, b435, b535, b635, b735, b835, b935=stdata35[:, 0], stdata35[:, 1], stdata35[:, 2], stdata35[:, 3], stdata35[:, 4], stdata35[:, 5], stdata35[:, 6], stdata35[:, 7], stdata35[:, 8], stdata35[:, 9] 
ax3.plot(b335/KM, b435/MSUN, 'o')
 
stdata36= np.genfromtxt('stgb_solver_pdata36.txt') 
b036, b136, b236, b336, b436, b536, b636, b736, b836, b936=stdata36[:, 0], stdata36[:, 1], stdata36[:, 2], stdata36[:, 3], stdata36[:, 4], stdata36[:, 5], stdata36[:, 6], stdata36[:, 7], stdata36[:, 8], stdata36[:, 9] 
ax4.plot(b336/KM, b436/MSUN, 'o') 
stdata37= np.genfromtxt('stgb_solver_pdata37.txt') 
b037, b137, b237, b337, b437, b537, b637, b737, b837, b937=stdata37[:, 0], stdata37[:, 1], stdata37[:, 2], stdata37[:, 3], stdata37[:, 4], stdata37[:, 5], stdata37[:, 6], stdata37[:, 7], stdata37[:, 8], stdata37[:, 9] 
ax4.plot(b337/KM, b437/MSUN, 'o') 
stdata38= np.genfromtxt('stgb_solver_pdata38.txt') 
b038, b138, b238, b338, b438, b538, b638, b738, b838, b938=stdata38[:, 0], stdata38[:, 1], stdata38[:, 2], stdata38[:, 3], stdata38[:, 4], stdata38[:, 5], stdata38[:, 6], stdata38[:, 7], stdata38[:, 8], stdata38[:, 9] 
ax4.plot(b338/KM, b438/MSUN, 'o') 
stdata39= np.genfromtxt('stgb_solver_pdata39.txt') 
b039, b139, b239, b339, b439, b539, b639, b739, b839, b939=stdata39[:, 0], stdata39[:, 1], stdata39[:, 2], stdata39[:, 3], stdata39[:, 4], stdata39[:, 5], stdata39[:, 6], stdata39[:, 7], stdata39[:, 8], stdata39[:, 9] 
ax4.plot(b339/KM, b439/MSUN, 'o') 
stdata40= np.genfromtxt('stgb_solver_pdata40.txt') 
b040, b140, b240, b340, b440, b540, b640, b740, b840, b940=stdata40[:, 0], stdata40[:, 1], stdata40[:, 2], stdata40[:, 3], stdata40[:, 4], stdata40[:, 5], stdata40[:, 6], stdata40[:, 7], stdata40[:, 8], stdata40[:, 9] 
ax4.plot(b340/KM, b440/MSUN, 'o') 






plt.savefig("fig_nonlin_MR.pdf", format='pdf', bbox_inches="tight")
print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show()
