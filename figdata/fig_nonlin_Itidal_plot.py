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
from scipy.optimize import curve_fit


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
def Idiml(i):
    """ dimensionless moment of inertia """
    return i*G/c**2 / runit**3
    
    
        

colorset=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple'] 

# ec, M, R, dimllatidal, y
ldata1 = np.genfromtxt('TOV_tidal_v1_data3.txt')
a10, a11, a12, a13, a14 = ldata1[:, 0], ldata1[:, 1], ldata1[:, 2], ldata1[:, 3], ldata1[:, 4]


 
def poly4f(x, a, b, c, d, e):
   return a + b*x + c*x**2 + d*x**3 + e*x**4

ntrim1 = np.argmax( a11 ) 

xset = np.log( a13[0: ntrim1] )
yset = poly4f(xset, 1.496, 0.05951, 0.02238, -6.953e-4, 8.345e-6)     # Yagi and Yunes' fit
yset1 = poly4f(xset, 1.493, 0.06433, 0.02104, -0.0005637, 3.947e-6)   # Rui's fit
#yset2 = poly4f(xset, 1.471, 0.06022, 0.01493, 5.933e-5, -1.007e-5)

xexpset = np.exp(xset)
yexpset = np.exp(yset)
yexpset1 = np.exp(yset1)


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



# GR I-tidal
"""
ax1.plot(xexpset, yexpset, '--', color = colorset[0], label='$a$ = 1')
ax1.plot(xexpset, yexpset1, '--', color = colorset[1])
ax2.plot(xset, yset, '--', color = colorset[0], label='$a$ = 1')
ax2.plot(xset, yset1, '--', color = colorset[1])
ax3.semilogx(xexpset, yexpset, '--', color = colorset[0], label='$a$ = 1')
ax3.semilogx(xexpset, yexpset1, '--', color = colorset[1])
ax4.loglog(xexpset, yexpset, '--', color = colorset[0], label='$a$ = 1')
ax4.loglog(xexpset, yexpset1, '--', color = colorset[1])
"""
ax1.loglog(xexpset, yexpset1, '--', color = 'black', label='$a$ = 1')
ax2.loglog(xexpset, yexpset1, '--', color = 'black')
ax3.loglog(xexpset, yexpset1, '--', color = 'black')
ax4.loglog(xexpset, yexpset1, '--', color = 'black')
ax5.loglog(xexpset, yexpset1, '--', color = 'black', label='$a$ = 0.1')
ax6.loglog(xexpset, yexpset1, '--', color = 'black')
ax7.loglog(xexpset, yexpset1, '--', color = 'black')
ax8.loglog(xexpset, yexpset1, '--', color = 'black')


ax1.legend()
ax5.legend()





# ec, bphc, M, ytidal, dimllatidal, dimllatidal1, dimllatidal2, dimllatidal3, sigma1, sigma2, sigma3, R, deR, r_inf, I
stdata1= np.genfromtxt('stgb_tid_test_data1.txt') 
c01, c11, c21, c31, c41, c51, c61, c71, c141=stdata1[:, 0], stdata1[:, 1], stdata1[:, 2], stdata1[:, 3], stdata1[:, 4], stdata1[:, 5], stdata1[:, 6], stdata1[:, 7], stdata1[:, 14] 
x1 = c51 
y1 = Idiml(c141)/(mdiml(c21))**3 
ax5.loglog(x1, y1, 'o') 
stdata2= np.genfromtxt('stgb_tid_test_data2.txt') 
c02, c12, c22, c32, c42, c52, c62, c72, c142=stdata2[:, 0], stdata2[:, 1], stdata2[:, 2], stdata2[:, 3], stdata2[:, 4], stdata2[:, 5], stdata2[:, 6], stdata2[:, 7], stdata2[:, 14] 
x2 = c52 
y2 = Idiml(c142)/(mdiml(c22))**3 
ax5.loglog(x2, y2, 'o') 
stdata3= np.genfromtxt('stgb_tid_test_data3.txt') 
c03, c13, c23, c33, c43, c53, c63, c73, c143=stdata3[:, 0], stdata3[:, 1], stdata3[:, 2], stdata3[:, 3], stdata3[:, 4], stdata3[:, 5], stdata3[:, 6], stdata3[:, 7], stdata3[:, 14] 
x3 = c53 
y3 = Idiml(c143)/(mdiml(c23))**3 
ax5.loglog(x3, y3, 'o') 
stdata4= np.genfromtxt('stgb_tid_test_data4.txt') 
c04, c14, c24, c34, c44, c54, c64, c74, c144=stdata4[:, 0], stdata4[:, 1], stdata4[:, 2], stdata4[:, 3], stdata4[:, 4], stdata4[:, 5], stdata4[:, 6], stdata4[:, 7], stdata4[:, 14] 
x4 = c54 
y4 = Idiml(c144)/(mdiml(c24))**3 
ax5.loglog(x4, y4, 'o') 
stdata5= np.genfromtxt('stgb_tid_test_data5.txt') 
c05, c15, c25, c35, c45, c55, c65, c75, c145=stdata5[:, 0], stdata5[:, 1], stdata5[:, 2], stdata5[:, 3], stdata5[:, 4], stdata5[:, 5], stdata5[:, 6], stdata5[:, 7], stdata5[:, 14] 
x5 = c55 
y5 = Idiml(c145)/(mdiml(c25))**3 
ax5.loglog(x5, y5, 'o') 
stdata6= np.genfromtxt('stgb_tid_test_data6.txt') 
c06, c16, c26, c36, c46, c56, c66, c76, c146=stdata6[:, 0], stdata6[:, 1], stdata6[:, 2], stdata6[:, 3], stdata6[:, 4], stdata6[:, 5], stdata6[:, 6], stdata6[:, 7], stdata6[:, 14] 
x6 = c56 
y6 = Idiml(c146)/(mdiml(c26))**3 
ax6.loglog(x6, y6, 'o') 
stdata7= np.genfromtxt('stgb_tid_test_data7.txt') 
c07, c17, c27, c37, c47, c57, c67, c77, c147=stdata7[:, 0], stdata7[:, 1], stdata7[:, 2], stdata7[:, 3], stdata7[:, 4], stdata7[:, 5], stdata7[:, 6], stdata7[:, 7], stdata7[:, 14] 
x7 = c57 
y7 = Idiml(c147)/(mdiml(c27))**3 
ax6.loglog(x7, y7, 'o') 
stdata8= np.genfromtxt('stgb_tid_test_data8.txt') 
c08, c18, c28, c38, c48, c58, c68, c78, c148=stdata8[:, 0], stdata8[:, 1], stdata8[:, 2], stdata8[:, 3], stdata8[:, 4], stdata8[:, 5], stdata8[:, 6], stdata8[:, 7], stdata8[:, 14] 
x8 = c58 
y8 = Idiml(c148)/(mdiml(c28))**3 
ax6.loglog(x8, y8, 'o') 
stdata9= np.genfromtxt('stgb_tid_test_data9.txt') 
c09, c19, c29, c39, c49, c59, c69, c79, c149=stdata9[:, 0], stdata9[:, 1], stdata9[:, 2], stdata9[:, 3], stdata9[:, 4], stdata9[:, 5], stdata9[:, 6], stdata9[:, 7], stdata9[:, 14] 
x9 = c59 
y9 = Idiml(c149)/(mdiml(c29))**3 
ax6.loglog(x9, y9, 'o') 
stdata10= np.genfromtxt('stgb_tid_test_data10.txt') 
c010, c110, c210, c310, c410, c510, c610, c710, c1410=stdata10[:, 0], stdata10[:, 1], stdata10[:, 2], stdata10[:, 3], stdata10[:, 4], stdata10[:, 5], stdata10[:, 6], stdata10[:, 7], stdata10[:, 14] 
x10 = c510 
y10 = Idiml(c1410)/(mdiml(c210))**3 
ax6.loglog(x10, y10, 'o') 
stdata11= np.genfromtxt('stgb_tid_test_data11.txt') 
c011, c111, c211, c311, c411, c511, c611, c711, c1411=stdata11[:, 0], stdata11[:, 1], stdata11[:, 2], stdata11[:, 3], stdata11[:, 4], stdata11[:, 5], stdata11[:, 6], stdata11[:, 7], stdata11[:, 14] 
x11 = c511 
y11 = Idiml(c1411)/(mdiml(c211))**3 
ax7.loglog(x11, y11, 'o') 
stdata12= np.genfromtxt('stgb_tid_test_data12.txt') 
c012, c112, c212, c312, c412, c512, c612, c712, c1412=stdata12[:, 0], stdata12[:, 1], stdata12[:, 2], stdata12[:, 3], stdata12[:, 4], stdata12[:, 5], stdata12[:, 6], stdata12[:, 7], stdata12[:, 14] 
x12 = c512 
y12 = Idiml(c1412)/(mdiml(c212))**3 
ax7.loglog(x12, y12, 'o') 
stdata13= np.genfromtxt('stgb_tid_test_data13.txt') 
c013, c113, c213, c313, c413, c513, c613, c713, c1413=stdata13[:, 0], stdata13[:, 1], stdata13[:, 2], stdata13[:, 3], stdata13[:, 4], stdata13[:, 5], stdata13[:, 6], stdata13[:, 7], stdata13[:, 14] 
x13 = c513 
y13 = Idiml(c1413)/(mdiml(c213))**3 
ax7.loglog(x13, y13, 'o') 
stdata14= np.genfromtxt('stgb_tid_test_data14.txt') 
c014, c114, c214, c314, c414, c514, c614, c714, c1414=stdata14[:, 0], stdata14[:, 1], stdata14[:, 2], stdata14[:, 3], stdata14[:, 4], stdata14[:, 5], stdata14[:, 6], stdata14[:, 7], stdata14[:, 14] 
x14 = c514 
y14 = Idiml(c1414)/(mdiml(c214))**3 
ax7.loglog(x14, y14, 'o') 
stdata15= np.genfromtxt('stgb_tid_test_data15.txt') 
c015, c115, c215, c315, c415, c515, c615, c715, c1415=stdata15[:, 0], stdata15[:, 1], stdata15[:, 2], stdata15[:, 3], stdata15[:, 4], stdata15[:, 5], stdata15[:, 6], stdata15[:, 7], stdata15[:, 14] 
x15 = c515 
y15 = Idiml(c1415)/(mdiml(c215))**3 
ax7.loglog(x15, y15, 'o') 
stdata16= np.genfromtxt('stgb_tid_test_data16.txt') 
c016, c116, c216, c316, c416, c516, c616, c716, c1416=stdata16[:, 0], stdata16[:, 1], stdata16[:, 2], stdata16[:, 3], stdata16[:, 4], stdata16[:, 5], stdata16[:, 6], stdata16[:, 7], stdata16[:, 14] 
x16 = c516 
y16 = Idiml(c1416)/(mdiml(c216))**3 
ax8.loglog(x16, y16, 'o') 
stdata17= np.genfromtxt('stgb_tid_test_data17.txt') 
c017, c117, c217, c317, c417, c517, c617, c717, c1417=stdata17[:, 0], stdata17[:, 1], stdata17[:, 2], stdata17[:, 3], stdata17[:, 4], stdata17[:, 5], stdata17[:, 6], stdata17[:, 7], stdata17[:, 14] 
x17 = c517 
y17 = Idiml(c1417)/(mdiml(c217))**3 
ax8.loglog(x17, y17, 'o') 
stdata18= np.genfromtxt('stgb_tid_test_data18.txt') 
c018, c118, c218, c318, c418, c518, c618, c718, c1418=stdata18[:, 0], stdata18[:, 1], stdata18[:, 2], stdata18[:, 3], stdata18[:, 4], stdata18[:, 5], stdata18[:, 6], stdata18[:, 7], stdata18[:, 14] 
x18 = c518 
y18 = Idiml(c1418)/(mdiml(c218))**3 
ax8.loglog(x18, y18, 'o') 
stdata19= np.genfromtxt('stgb_tid_test_data19.txt') 
c019, c119, c219, c319, c419, c519, c619, c719, c1419=stdata19[:, 0], stdata19[:, 1], stdata19[:, 2], stdata19[:, 3], stdata19[:, 4], stdata19[:, 5], stdata19[:, 6], stdata19[:, 7], stdata19[:, 14] 
x19 = c519 
y19 = Idiml(c1419)/(mdiml(c219))**3 
ax8.loglog(x19, y19, 'o') 
stdata20= np.genfromtxt('stgb_tid_test_data20.txt') 
c020, c120, c220, c320, c420, c520, c620, c720, c1420=stdata20[:, 0], stdata20[:, 1], stdata20[:, 2], stdata20[:, 3], stdata20[:, 4], stdata20[:, 5], stdata20[:, 6], stdata20[:, 7], stdata20[:, 14] 
x20 = c520 
y20 = Idiml(c1420)/(mdiml(c220))**3 
ax8.loglog(x20, y20, 'o') 
stdata21= np.genfromtxt('stgb_tid_test_data21.txt') 
c021, c121, c221, c321, c421, c521, c621, c721, c1421=stdata21[:, 0], stdata21[:, 1], stdata21[:, 2], stdata21[:, 3], stdata21[:, 4], stdata21[:, 5], stdata21[:, 6], stdata21[:, 7], stdata21[:, 14] 
x21 = c521 
y21 = Idiml(c1421)/(mdiml(c221))**3 
ax1.loglog(x21, y21, 'o') 
stdata22= np.genfromtxt('stgb_tid_test_data22.txt') 
c022, c122, c222, c322, c422, c522, c622, c722, c1422=stdata22[:, 0], stdata22[:, 1], stdata22[:, 2], stdata22[:, 3], stdata22[:, 4], stdata22[:, 5], stdata22[:, 6], stdata22[:, 7], stdata22[:, 14] 
x22 = c522 
y22 = Idiml(c1422)/(mdiml(c222))**3 
ax1.loglog(x22, y22, 'o') 
stdata23= np.genfromtxt('stgb_tid_test_data23.txt') 
c023, c123, c223, c323, c423, c523, c623, c723, c1423=stdata23[:, 0], stdata23[:, 1], stdata23[:, 2], stdata23[:, 3], stdata23[:, 4], stdata23[:, 5], stdata23[:, 6], stdata23[:, 7], stdata23[:, 14] 
x23 = c523 
y23 = Idiml(c1423)/(mdiml(c223))**3 
ax1.loglog(x23, y23, 'o') 
stdata24= np.genfromtxt('stgb_tid_test_data24.txt') 
c024, c124, c224, c324, c424, c524, c624, c724, c1424=stdata24[:, 0], stdata24[:, 1], stdata24[:, 2], stdata24[:, 3], stdata24[:, 4], stdata24[:, 5], stdata24[:, 6], stdata24[:, 7], stdata24[:, 14] 
x24 = c524 
y24 = Idiml(c1424)/(mdiml(c224))**3 
ax1.loglog(x24, y24, 'o') 
stdata25= np.genfromtxt('stgb_tid_test_data25.txt') 
c025, c125, c225, c325, c425, c525, c625, c725, c1425=stdata25[:, 0], stdata25[:, 1], stdata25[:, 2], stdata25[:, 3], stdata25[:, 4], stdata25[:, 5], stdata25[:, 6], stdata25[:, 7], stdata25[:, 14] 
x25 = c525 
y25 = Idiml(c1425)/(mdiml(c225))**3 
ax1.loglog(x25, y25, 'o') 
stdata26= np.genfromtxt('stgb_tid_test_data26.txt') 
c026, c126, c226, c326, c426, c526, c626, c726, c1426=stdata26[:, 0], stdata26[:, 1], stdata26[:, 2], stdata26[:, 3], stdata26[:, 4], stdata26[:, 5], stdata26[:, 6], stdata26[:, 7], stdata26[:, 14] 
x26 = c526 
y26 = Idiml(c1426)/(mdiml(c226))**3 
ax2.loglog(x26, y26, 'o') 
stdata27= np.genfromtxt('stgb_tid_test_data27.txt') 
c027, c127, c227, c327, c427, c527, c627, c727, c1427=stdata27[:, 0], stdata27[:, 1], stdata27[:, 2], stdata27[:, 3], stdata27[:, 4], stdata27[:, 5], stdata27[:, 6], stdata27[:, 7], stdata27[:, 14] 
x27 = c527 
y27 = Idiml(c1427)/(mdiml(c227))**3 
ax2.loglog(x27, y27, 'o') 
stdata28= np.genfromtxt('stgb_tid_test_data28.txt') 
c028, c128, c228, c328, c428, c528, c628, c728, c1428=stdata28[:, 0], stdata28[:, 1], stdata28[:, 2], stdata28[:, 3], stdata28[:, 4], stdata28[:, 5], stdata28[:, 6], stdata28[:, 7], stdata28[:, 14] 
x28 = c528 
y28 = Idiml(c1428)/(mdiml(c228))**3 
ax2.loglog(x28, y28, 'o') 
#stdata29= np.genfromtxt('stgb_tid_test_data29.txt') 
#c029, c129, c229, c329, c429, c529, c629, c729, c1429=stdata29[:, 0], stdata29[:, 1], stdata29[:, 2], stdata29[:, 3], stdata29[:, 4], stdata29[:, 5], stdata29[:, 6], stdata29[:, 7], stdata29[:, 14] 
#x29 = c529 
#y29 = Idiml(c1429)/(mdiml(c229))**3 
#ax2.loglog(x29, y29, 'o') 
stdata30= np.genfromtxt('stgb_tid_test_data30.txt') 
c030, c130, c230, c330, c430, c530, c630, c730, c1430=stdata30[:, 0], stdata30[:, 1], stdata30[:, 2], stdata30[:, 3], stdata30[:, 4], stdata30[:, 5], stdata30[:, 6], stdata30[:, 7], stdata30[:, 14] 
x30 = c530 
y30 = Idiml(c1430)/(mdiml(c230))**3 
ax2.loglog(x30, y30, 'o') 
stdata31= np.genfromtxt('stgb_tid_test_data31.txt') 
c031, c131, c231, c331, c431, c531, c631, c731, c1431=stdata31[:, 0], stdata31[:, 1], stdata31[:, 2], stdata31[:, 3], stdata31[:, 4], stdata31[:, 5], stdata31[:, 6], stdata31[:, 7], stdata31[:, 14] 
x31 = c531 
y31 = Idiml(c1431)/(mdiml(c231))**3 
ax3.loglog(x31, y31, 'o') 
stdata32= np.genfromtxt('stgb_tid_test_data32.txt') 
c032, c132, c232, c332, c432, c532, c632, c732, c1432=stdata32[:, 0], stdata32[:, 1], stdata32[:, 2], stdata32[:, 3], stdata32[:, 4], stdata32[:, 5], stdata32[:, 6], stdata32[:, 7], stdata32[:, 14] 
x32 = c532 
y32 = Idiml(c1432)/(mdiml(c232))**3 
ax3.loglog(x32, y32, 'o') 
stdata33= np.genfromtxt('stgb_tid_test_data33.txt') 
c033, c133, c233, c333, c433, c533, c633, c733, c1433=stdata33[:, 0], stdata33[:, 1], stdata33[:, 2], stdata33[:, 3], stdata33[:, 4], stdata33[:, 5], stdata33[:, 6], stdata33[:, 7], stdata33[:, 14] 
x33 = c533 
y33 = Idiml(c1433)/(mdiml(c233))**3 
ax3.loglog(x33, y33, 'o') 
stdata34= np.genfromtxt('stgb_tid_test_data34.txt') 
c034, c134, c234, c334, c434, c534, c634, c734, c1434=stdata34[:, 0], stdata34[:, 1], stdata34[:, 2], stdata34[:, 3], stdata34[:, 4], stdata34[:, 5], stdata34[:, 6], stdata34[:, 7], stdata34[:, 14] 
x34 = c534 
y34 = Idiml(c1434)/(mdiml(c234))**3 
ax3.loglog(x34, y34, 'o') 
stdata35= np.genfromtxt('stgb_tid_test_data35.txt') 
c035, c135, c235, c335, c435, c535, c635, c735, c1435=stdata35[:, 0], stdata35[:, 1], stdata35[:, 2], stdata35[:, 3], stdata35[:, 4], stdata35[:, 5], stdata35[:, 6], stdata35[:, 7], stdata35[:, 14] 
x35 = c535 
y35 = Idiml(c1435)/(mdiml(c235))**3 
ax3.loglog(x35, y35, 'o') 

stdata36= np.genfromtxt('stgb_tid_test_data36.txt') 
c036, c136, c236, c336, c436, c536, c636, c736, c1436=stdata36[:, 0], stdata36[:, 1], stdata36[:, 2], stdata36[:, 3], stdata36[:, 4], stdata36[:, 5], stdata36[:, 6], stdata36[:, 7], stdata36[:, 14] 
x36 = c736 
y36 = Idiml(c1436)/(mdiml(c236))**3 
ax4.loglog(x36, y36, 'o') 
stdata37= np.genfromtxt('stgb_tid_test_data37.txt') 
c037, c137, c237, c337, c437, c537, c637, c737, c1437=stdata37[:, 0], stdata37[:, 1], stdata37[:, 2], stdata37[:, 3], stdata37[:, 4], stdata37[:, 5], stdata37[:, 6], stdata37[:, 7], stdata37[:, 14] 
x37 = c737 
y37 = Idiml(c1437)/(mdiml(c237))**3 
ax4.loglog(x37, y37, 'o') 
stdata38= np.genfromtxt('stgb_tid_test_data38.txt') 
c038, c138, c238, c338, c438, c538, c638, c738, c1438=stdata38[:, 0], stdata38[:, 1], stdata38[:, 2], stdata38[:, 3], stdata38[:, 4], stdata38[:, 5], stdata38[:, 6], stdata38[:, 7], stdata38[:, 14] 
x38 = c738 
y38 = Idiml(c1438)/(mdiml(c238))**3 
ax4.loglog(x38, y38, 'o') 
stdata39= np.genfromtxt('stgb_tid_test_data39.txt') 
c039, c139, c239, c339, c439, c539, c639, c739, c1439=stdata39[:, 0], stdata39[:, 1], stdata39[:, 2], stdata39[:, 3], stdata39[:, 4], stdata39[:, 5], stdata39[:, 6], stdata39[:, 7], stdata39[:, 14] 
x39 = c739 
y39 = Idiml(c1439)/(mdiml(c239))**3 
ax4.loglog(x39, y39, 'o') 
stdata40= np.genfromtxt('stgb_tid_test_data40.txt') 
c040, c140, c240, c340, c440, c540, c640, c740, c1440=stdata40[:, 0], stdata40[:, 1], stdata40[:, 2], stdata40[:, 3], stdata40[:, 4], stdata40[:, 5], stdata40[:, 6], stdata40[:, 7], stdata40[:, 14] 
x40 = c740 
y40 = Idiml(c1440)/(mdiml(c240))**3 
ax4.loglog(x40, y40, 'o')





plt.savefig("fig_nonlin_Itidal.pdf", format='pdf', bbox_inches="tight")
print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show() 
