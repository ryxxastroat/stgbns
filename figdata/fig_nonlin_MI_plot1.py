"""
test code: plot M-I
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
def idiml(i, m):
    """ i/m^3 """
    return i*c**4/m**3/G**2

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
 
# GR M-I
xmin, xmax = 0.7, 2.6
ymin, ymax = 3, 30
ax1.set_xlim([xmin, xmax])
ax2.set_xlim([xmin, xmax])
ax5.set_xlim([xmin, xmax])
ax6.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])
ax2.set_ylim([ymin, ymax])
ax5.set_ylim([ymin, ymax])
ax6.set_ylim([ymin, ymax])



ax1.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0], label='$a$ = 1')
ax1.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax1.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax1.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax1.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])
ax2.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0])
ax2.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax2.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax2.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax2.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])
ax3.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0])
ax3.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax3.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax3.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax3.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])
ax4.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0])
ax4.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax4.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax4.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax4.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])
ax5.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0], label='$a$ = 0.1')
ax5.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax5.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax5.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax5.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])
ax6.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0])
ax6.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax6.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax6.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax6.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])
ax7.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0])
ax7.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax7.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax7.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax7.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])
ax8.plot(a11/MSUN, idiml(a14, a11), '--', color = colorset[0])
ax8.plot(a21/MSUN, idiml(a24, a21), '--', color = colorset[1])
ax8.plot(a31/MSUN, idiml(a34, a31), '--', color = colorset[2])
ax8.plot(a41/MSUN, idiml(a44, a41), '--', color = colorset[3])
ax8.plot(a51/MSUN, idiml(a54, a51), '--', color = colorset[4])

ax1.legend()
ax5.legend()



#ntrim1, ntrim2, ntrim3, ntrim4, ntrim5 = len(b01), len(b02), len(b03), len(b04), len(b05)
#ntrim10 = len(b010)

#ntrim21, ntrim22, ntrim23, ntrim24, ntrim25 = len(b021), len(b022), len(b023), len(b024), len(b025)
#ntrim30 = len(b030)





# ST M-I
stdata1= np.genfromtxt('stgb_tid_v1_comb_t3data1.txt') 
c01, c11, c21, c31, c41, c51, c61, c71, c81, c91, c101=stdata1[:, 0], stdata1[:, 1], stdata1[:, 2], stdata1[:, 3], stdata1[:, 4], stdata1[:, 5], stdata1[:, 6], stdata1[:, 7], stdata1[:, 8], stdata1[:, 9], stdata1[:, 10] 
x1 = c31/MSUN 
y1 = idiml(c41, c31) 
ax5.semilogy(x1, y1, 'o', color = colorset[0]) 
stdata2= np.genfromtxt('stgb_tid_v1_comb_t3data2.txt') 
c02, c12, c22, c32, c42, c52, c62, c72, c82, c92, c102=stdata2[:, 0], stdata2[:, 1], stdata2[:, 2], stdata2[:, 3], stdata2[:, 4], stdata2[:, 5], stdata2[:, 6], stdata2[:, 7], stdata2[:, 8], stdata2[:, 9], stdata2[:, 10] 
x2 = c32/MSUN 
y2 = idiml(c42, c32) 
ax5.semilogy(x2, y2, 'o', color = colorset[1]) 
stdata3= np.genfromtxt('stgb_tid_v1_comb_t3data3.txt') 
c03, c13, c23, c33, c43, c53, c63, c73, c83, c93, c103=stdata3[:, 0], stdata3[:, 1], stdata3[:, 2], stdata3[:, 3], stdata3[:, 4], stdata3[:, 5], stdata3[:, 6], stdata3[:, 7], stdata3[:, 8], stdata3[:, 9], stdata3[:, 10] 
x3 = c33/MSUN 
y3 = idiml(c43, c33) 
ax5.semilogy(x3, y3, 'o', color = colorset[2]) 
stdata4= np.genfromtxt('stgb_tid_v1_comb_t3data4.txt') 
c04, c14, c24, c34, c44, c54, c64, c74, c84, c94, c104=stdata4[:, 0], stdata4[:, 1], stdata4[:, 2], stdata4[:, 3], stdata4[:, 4], stdata4[:, 5], stdata4[:, 6], stdata4[:, 7], stdata4[:, 8], stdata4[:, 9], stdata4[:, 10] 
x4 = c34/MSUN 
y4 = idiml(c44, c34) 
ax5.semilogy(x4, y4, 'o', color = colorset[3]) 
stdata5= np.genfromtxt('stgb_tid_v1_comb_t3data5.txt') 
c05, c15, c25, c35, c45, c55, c65, c75, c85, c95, c105=stdata5[:, 0], stdata5[:, 1], stdata5[:, 2], stdata5[:, 3], stdata5[:, 4], stdata5[:, 5], stdata5[:, 6], stdata5[:, 7], stdata5[:, 8], stdata5[:, 9], stdata5[:, 10] 
x5 = c35/MSUN 
y5 = idiml(c45, c35) 
ax5.semilogy(x5, y5, 'o', color = colorset[4]) 
stdata6= np.genfromtxt('stgb_tid_v1_comb_t3data6.txt') 
c06, c16, c26, c36, c46, c56, c66, c76, c86, c96, c106=stdata6[:, 0], stdata6[:, 1], stdata6[:, 2], stdata6[:, 3], stdata6[:, 4], stdata6[:, 5], stdata6[:, 6], stdata6[:, 7], stdata6[:, 8], stdata6[:, 9], stdata6[:, 10] 
x6 = c36/MSUN 
y6 = idiml(c46, c36) 
ax6.semilogy(x6, y6, 'o', color = colorset[0]) 
stdata7= np.genfromtxt('stgb_tid_v1_comb_t3data7.txt') 
c07, c17, c27, c37, c47, c57, c67, c77, c87, c97, c107=stdata7[:, 0], stdata7[:, 1], stdata7[:, 2], stdata7[:, 3], stdata7[:, 4], stdata7[:, 5], stdata7[:, 6], stdata7[:, 7], stdata7[:, 8], stdata7[:, 9], stdata7[:, 10] 
x7 = c37/MSUN 
y7 = idiml(c47, c37) 
ax6.semilogy(x7, y7, 'o', color = colorset[1]) 
stdata8= np.genfromtxt('stgb_tid_v1_comb_t3data8.txt') 
c08, c18, c28, c38, c48, c58, c68, c78, c88, c98, c108=stdata8[:, 0], stdata8[:, 1], stdata8[:, 2], stdata8[:, 3], stdata8[:, 4], stdata8[:, 5], stdata8[:, 6], stdata8[:, 7], stdata8[:, 8], stdata8[:, 9], stdata8[:, 10] 
x8 = c38/MSUN 
y8 = idiml(c48, c38) 
ax6.semilogy(x8, y8, 'o', color = colorset[2]) 
stdata9= np.genfromtxt('stgb_tid_v1_comb_t3data9.txt') 
c09, c19, c29, c39, c49, c59, c69, c79, c89, c99, c109=stdata9[:, 0], stdata9[:, 1], stdata9[:, 2], stdata9[:, 3], stdata9[:, 4], stdata9[:, 5], stdata9[:, 6], stdata9[:, 7], stdata9[:, 8], stdata9[:, 9], stdata9[:, 10] 
x9 = c39/MSUN 
y9 = idiml(c49, c39) 
ax6.semilogy(x9, y9, 'o', color = colorset[3]) 
stdata10= np.genfromtxt('stgb_tid_v1_comb_t3data10.txt') 
c010, c110, c210, c310, c410, c510, c610, c710, c810, c910, c1010=stdata10[:, 0], stdata10[:, 1], stdata10[:, 2], stdata10[:, 3], stdata10[:, 4], stdata10[:, 5], stdata10[:, 6], stdata10[:, 7], stdata10[:, 8], stdata10[:, 9], stdata10[:, 10] 
x10 = c310/MSUN 
y10 = idiml(c410, c310) 
ax6.semilogy(x10, y10, 'o', color = colorset[4]) 
stdata11= np.genfromtxt('stgb_tid_v1_comb_t3data11.txt') 
c011, c111, c211, c311, c411, c511, c611, c711, c811, c911, c1011=stdata11[:, 0], stdata11[:, 1], stdata11[:, 2], stdata11[:, 3], stdata11[:, 4], stdata11[:, 5], stdata11[:, 6], stdata11[:, 7], stdata11[:, 8], stdata11[:, 9], stdata11[:, 10] 
x11 = c311/MSUN 
y11 = idiml(c411, c311) 
ax7.semilogy(x11, y11, 'o', color = colorset[0]) 
stdata12= np.genfromtxt('stgb_tid_v1_comb_t3data12.txt') 
c012, c112, c212, c312, c412, c512, c612, c712, c812, c912, c1012=stdata12[:, 0], stdata12[:, 1], stdata12[:, 2], stdata12[:, 3], stdata12[:, 4], stdata12[:, 5], stdata12[:, 6], stdata12[:, 7], stdata12[:, 8], stdata12[:, 9], stdata12[:, 10] 
x12 = c312/MSUN 
y12 = idiml(c412, c312) 
ax7.semilogy(x12, y12, 'o', color = colorset[1]) 
stdata13= np.genfromtxt('stgb_tid_v1_comb_t3data13.txt') 
c013, c113, c213, c313, c413, c513, c613, c713, c813, c913, c1013=stdata13[:, 0], stdata13[:, 1], stdata13[:, 2], stdata13[:, 3], stdata13[:, 4], stdata13[:, 5], stdata13[:, 6], stdata13[:, 7], stdata13[:, 8], stdata13[:, 9], stdata13[:, 10] 
x13 = c313/MSUN 
y13 = idiml(c413, c313) 
ax7.semilogy(x13, y13, 'o', color = colorset[2]) 
stdata14= np.genfromtxt('stgb_tid_v1_comb_t3data14.txt') 
c014, c114, c214, c314, c414, c514, c614, c714, c814, c914, c1014=stdata14[:, 0], stdata14[:, 1], stdata14[:, 2], stdata14[:, 3], stdata14[:, 4], stdata14[:, 5], stdata14[:, 6], stdata14[:, 7], stdata14[:, 8], stdata14[:, 9], stdata14[:, 10] 
x14 = c314/MSUN 
y14 = idiml(c414, c314) 
ax7.semilogy(x14, y14, 'o', color = colorset[3]) 
stdata15= np.genfromtxt('stgb_tid_v1_comb_t3data15.txt') 
c015, c115, c215, c315, c415, c515, c615, c715, c815, c915, c1015=stdata15[:, 0], stdata15[:, 1], stdata15[:, 2], stdata15[:, 3], stdata15[:, 4], stdata15[:, 5], stdata15[:, 6], stdata15[:, 7], stdata15[:, 8], stdata15[:, 9], stdata15[:, 10] 
x15 = c315/MSUN 
y15 = idiml(c415, c315) 
ax7.semilogy(x15, y15, 'o', color = colorset[4]) 
stdata16= np.genfromtxt('stgb_tid_v1_comb_t3data16.txt') 
c016, c116, c216, c316, c416, c516, c616, c716, c816, c916, c1016=stdata16[:, 0], stdata16[:, 1], stdata16[:, 2], stdata16[:, 3], stdata16[:, 4], stdata16[:, 5], stdata16[:, 6], stdata16[:, 7], stdata16[:, 8], stdata16[:, 9], stdata16[:, 10] 
x16 = c316/MSUN 
y16 = idiml(c416, c316) 
ax8.semilogy(x16, y16, 'o', color = colorset[0]) 
stdata17= np.genfromtxt('stgb_tid_v1_comb_t3data17.txt') 
c017, c117, c217, c317, c417, c517, c617, c717, c817, c917, c1017=stdata17[:, 0], stdata17[:, 1], stdata17[:, 2], stdata17[:, 3], stdata17[:, 4], stdata17[:, 5], stdata17[:, 6], stdata17[:, 7], stdata17[:, 8], stdata17[:, 9], stdata17[:, 10] 
x17 = c317/MSUN 
y17 = idiml(c417, c317) 
ax8.semilogy(x17, y17, 'o', color = colorset[1]) 
stdata18= np.genfromtxt('stgb_tid_v1_comb_t3data18.txt') 
c018, c118, c218, c318, c418, c518, c618, c718, c818, c918, c1018=stdata18[:, 0], stdata18[:, 1], stdata18[:, 2], stdata18[:, 3], stdata18[:, 4], stdata18[:, 5], stdata18[:, 6], stdata18[:, 7], stdata18[:, 8], stdata18[:, 9], stdata18[:, 10] 
x18 = c318/MSUN 
y18 = idiml(c418, c318) 
ax8.semilogy(x18, y18, 'o', color = colorset[2]) 
stdata19= np.genfromtxt('stgb_tid_v1_comb_t3data19.txt') 
c019, c119, c219, c319, c419, c519, c619, c719, c819, c919, c1019=stdata19[:, 0], stdata19[:, 1], stdata19[:, 2], stdata19[:, 3], stdata19[:, 4], stdata19[:, 5], stdata19[:, 6], stdata19[:, 7], stdata19[:, 8], stdata19[:, 9], stdata19[:, 10] 
x19 = c319/MSUN 
y19 = idiml(c419, c319) 
ax8.semilogy(x19, y19, 'o', color = colorset[3]) 
stdata20= np.genfromtxt('stgb_tid_v1_comb_t3data20.txt') 
c020, c120, c220, c320, c420, c520, c620, c720, c820, c920, c1020=stdata20[:, 0], stdata20[:, 1], stdata20[:, 2], stdata20[:, 3], stdata20[:, 4], stdata20[:, 5], stdata20[:, 6], stdata20[:, 7], stdata20[:, 8], stdata20[:, 9], stdata20[:, 10] 
x20 = c320/MSUN 
y20 = idiml(c420, c320) 
ax8.semilogy(x20, y20, 'o', color = colorset[4]) 
stdata21= np.genfromtxt('stgb_tid_v1_comb_t3data21.txt') 
c021, c121, c221, c321, c421, c521, c621, c721, c821, c921, c1021=stdata21[:, 0], stdata21[:, 1], stdata21[:, 2], stdata21[:, 3], stdata21[:, 4], stdata21[:, 5], stdata21[:, 6], stdata21[:, 7], stdata21[:, 8], stdata21[:, 9], stdata21[:, 10] 
x21 = c321/MSUN 
y21 = idiml(c421, c321) 
ax1.semilogy(x21, y21, 'o', color = colorset[0]) 
stdata22= np.genfromtxt('stgb_tid_v1_comb_t3data22.txt') 
c022, c122, c222, c322, c422, c522, c622, c722, c822, c922, c1022=stdata22[:, 0], stdata22[:, 1], stdata22[:, 2], stdata22[:, 3], stdata22[:, 4], stdata22[:, 5], stdata22[:, 6], stdata22[:, 7], stdata22[:, 8], stdata22[:, 9], stdata22[:, 10] 
x22 = c322/MSUN 
y22 = idiml(c422, c322) 
ax1.semilogy(x22, y22, 'o', color = colorset[1]) 
stdata23= np.genfromtxt('stgb_tid_v1_comb_t3data23.txt') 
c023, c123, c223, c323, c423, c523, c623, c723, c823, c923, c1023=stdata23[:, 0], stdata23[:, 1], stdata23[:, 2], stdata23[:, 3], stdata23[:, 4], stdata23[:, 5], stdata23[:, 6], stdata23[:, 7], stdata23[:, 8], stdata23[:, 9], stdata23[:, 10] 
x23 = c323/MSUN 
y23 = idiml(c423, c323) 
ax1.semilogy(x23, y23, 'o', color = colorset[2]) 
stdata24= np.genfromtxt('stgb_tid_v1_comb_t3data24.txt') 
c024, c124, c224, c324, c424, c524, c624, c724, c824, c924, c1024=stdata24[:, 0], stdata24[:, 1], stdata24[:, 2], stdata24[:, 3], stdata24[:, 4], stdata24[:, 5], stdata24[:, 6], stdata24[:, 7], stdata24[:, 8], stdata24[:, 9], stdata24[:, 10] 
x24 = c324/MSUN 
y24 = idiml(c424, c324) 
ax1.semilogy(x24, y24, 'o', color = colorset[3]) 
stdata25= np.genfromtxt('stgb_tid_v1_comb_t3data25.txt') 
c025, c125, c225, c325, c425, c525, c625, c725, c825, c925, c1025=stdata25[:, 0], stdata25[:, 1], stdata25[:, 2], stdata25[:, 3], stdata25[:, 4], stdata25[:, 5], stdata25[:, 6], stdata25[:, 7], stdata25[:, 8], stdata25[:, 9], stdata25[:, 10] 
x25 = c325/MSUN 
y25 = idiml(c425, c325) 
ax1.semilogy(x25, y25, 'o', color = colorset[4]) 
stdata26= np.genfromtxt('stgb_tid_v1_comb_t3data26.txt') 
c026, c126, c226, c326, c426, c526, c626, c726, c826, c926, c1026=stdata26[:, 0], stdata26[:, 1], stdata26[:, 2], stdata26[:, 3], stdata26[:, 4], stdata26[:, 5], stdata26[:, 6], stdata26[:, 7], stdata26[:, 8], stdata26[:, 9], stdata26[:, 10] 
x26 = c326/MSUN 
y26 = idiml(c426, c326) 
ax2.semilogy(x26, y26, 'o', color = colorset[0]) 
stdata27= np.genfromtxt('stgb_tid_v1_comb_t3data27.txt') 
c027, c127, c227, c327, c427, c527, c627, c727, c827, c927, c1027=stdata27[:, 0], stdata27[:, 1], stdata27[:, 2], stdata27[:, 3], stdata27[:, 4], stdata27[:, 5], stdata27[:, 6], stdata27[:, 7], stdata27[:, 8], stdata27[:, 9], stdata27[:, 10] 
x27 = c327/MSUN 
y27 = idiml(c427, c327) 
ax2.semilogy(x27, y27, 'o', color = colorset[1]) 
stdata28= np.genfromtxt('stgb_tid_v1_comb_t3data28.txt') 
c028, c128, c228, c328, c428, c528, c628, c728, c828, c928, c1028=stdata28[:, 0], stdata28[:, 1], stdata28[:, 2], stdata28[:, 3], stdata28[:, 4], stdata28[:, 5], stdata28[:, 6], stdata28[:, 7], stdata28[:, 8], stdata28[:, 9], stdata28[:, 10] 
x28 = c328/MSUN 
y28 = idiml(c428, c328) 
ax2.semilogy(x28, y28, 'o', color = colorset[2]) 
stdata30= np.genfromtxt('stgb_tid_v1_comb_t3data30.txt') 
c030, c130, c230, c330, c430, c530, c630, c730, c830, c930, c1030=stdata30[:, 0], stdata30[:, 1], stdata30[:, 2], stdata30[:, 3], stdata30[:, 4], stdata30[:, 5], stdata30[:, 6], stdata30[:, 7], stdata30[:, 8], stdata30[:, 9], stdata30[:, 10] 
x30 = c330/MSUN 
y30 = idiml(c430, c330) 
ax2.semilogy(x30, y30, 'o', color = colorset[4]) 
stdata31= np.genfromtxt('stgb_tid_v1_comb_t3data31.txt') 
c031, c131, c231, c331, c431, c531, c631, c731, c831, c931, c1031=stdata31[:, 0], stdata31[:, 1], stdata31[:, 2], stdata31[:, 3], stdata31[:, 4], stdata31[:, 5], stdata31[:, 6], stdata31[:, 7], stdata31[:, 8], stdata31[:, 9], stdata31[:, 10] 
x31 = c331/MSUN 
y31 = idiml(c431, c331) 
ax3.semilogy(x31, y31, 'o', color = colorset[0]) 
stdata32= np.genfromtxt('stgb_tid_v1_comb_t3data32.txt') 
c032, c132, c232, c332, c432, c532, c632, c732, c832, c932, c1032=stdata32[:, 0], stdata32[:, 1], stdata32[:, 2], stdata32[:, 3], stdata32[:, 4], stdata32[:, 5], stdata32[:, 6], stdata32[:, 7], stdata32[:, 8], stdata32[:, 9], stdata32[:, 10] 
x32 = c332/MSUN 
y32 = idiml(c432, c332) 
ax3.semilogy(x32, y32, 'o', color = colorset[1]) 
stdata33= np.genfromtxt('stgb_tid_v1_comb_t3data33.txt') 
c033, c133, c233, c333, c433, c533, c633, c733, c833, c933, c1033=stdata33[:, 0], stdata33[:, 1], stdata33[:, 2], stdata33[:, 3], stdata33[:, 4], stdata33[:, 5], stdata33[:, 6], stdata33[:, 7], stdata33[:, 8], stdata33[:, 9], stdata33[:, 10] 
x33 = c333/MSUN 
y33 = idiml(c433, c333) 
ax3.semilogy(x33, y33, 'o', color = colorset[2]) 
stdata34= np.genfromtxt('stgb_tid_v1_comb_t3data34.txt') 
c034, c134, c234, c334, c434, c534, c634, c734, c834, c934, c1034=stdata34[:, 0], stdata34[:, 1], stdata34[:, 2], stdata34[:, 3], stdata34[:, 4], stdata34[:, 5], stdata34[:, 6], stdata34[:, 7], stdata34[:, 8], stdata34[:, 9], stdata34[:, 10] 
x34 = c334/MSUN 
y34 = idiml(c434, c334) 
ax3.semilogy(x34, y34, 'o', color = colorset[3]) 
stdata35= np.genfromtxt('stgb_tid_v1_comb_t3data35.txt') 
c035, c135, c235, c335, c435, c535, c635, c735, c835, c935, c1035=stdata35[:, 0], stdata35[:, 1], stdata35[:, 2], stdata35[:, 3], stdata35[:, 4], stdata35[:, 5], stdata35[:, 6], stdata35[:, 7], stdata35[:, 8], stdata35[:, 9], stdata35[:, 10] 
x35 = c335/MSUN 
y35 = idiml(c435, c335) 
ax3.semilogy(x35, y35, 'o', color = colorset[4]) 
stdata36= np.genfromtxt('stgb_tid_v1_comb_t3data36.txt') 
c036, c136, c236, c336, c436, c536, c636, c736, c836, c936, c1036=stdata36[:, 0], stdata36[:, 1], stdata36[:, 2], stdata36[:, 3], stdata36[:, 4], stdata36[:, 5], stdata36[:, 6], stdata36[:, 7], stdata36[:, 8], stdata36[:, 9], stdata36[:, 10] 
x36 = c336/MSUN 
y36 = idiml(c436, c336) 
ax4.semilogy(x36, y36, 'o', color = colorset[0]) 
stdata37= np.genfromtxt('stgb_tid_v1_comb_t3data37.txt') 
c037, c137, c237, c337, c437, c537, c637, c737, c837, c937, c1037=stdata37[:, 0], stdata37[:, 1], stdata37[:, 2], stdata37[:, 3], stdata37[:, 4], stdata37[:, 5], stdata37[:, 6], stdata37[:, 7], stdata37[:, 8], stdata37[:, 9], stdata37[:, 10] 
x37 = c337/MSUN 
y37 = idiml(c437, c337) 
ax4.semilogy(x37, y37, 'o', color = colorset[1]) 
stdata38= np.genfromtxt('stgb_tid_v1_comb_t3data38.txt') 
c038, c138, c238, c338, c438, c538, c638, c738, c838, c938, c1038=stdata38[:, 0], stdata38[:, 1], stdata38[:, 2], stdata38[:, 3], stdata38[:, 4], stdata38[:, 5], stdata38[:, 6], stdata38[:, 7], stdata38[:, 8], stdata38[:, 9], stdata38[:, 10] 
x38 = c338/MSUN 
y38 = idiml(c438, c338) 
ax4.semilogy(x38, y38, 'o', color = colorset[2]) 
stdata39= np.genfromtxt('stgb_tid_v1_comb_t3data39.txt') 
c039, c139, c239, c339, c439, c539, c639, c739, c839, c939, c1039=stdata39[:, 0], stdata39[:, 1], stdata39[:, 2], stdata39[:, 3], stdata39[:, 4], stdata39[:, 5], stdata39[:, 6], stdata39[:, 7], stdata39[:, 8], stdata39[:, 9], stdata39[:, 10] 
x39 = c339/MSUN 
y39 = idiml(c439, c339) 
ax4.semilogy(x39, y39, 'o', color = colorset[3]) 
stdata40= np.genfromtxt('stgb_tid_v1_comb_t3data40.txt') 
c040, c140, c240, c340, c440, c540, c640, c740, c840, c940, c1040=stdata40[:, 0], stdata40[:, 1], stdata40[:, 2], stdata40[:, 3], stdata40[:, 4], stdata40[:, 5], stdata40[:, 6], stdata40[:, 7], stdata40[:, 8], stdata40[:, 9], stdata40[:, 10] 
x40 = c340/MSUN 
y40 = idiml(c440, c340) 
ax4.semilogy(x40, y40, 'o', color = colorset[4]) 







plt.savefig("fig_nonlin_MI.pdf", format='pdf', bbox_inches="tight")
print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show()
