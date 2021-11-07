"""
solve the moment of inertia of the scalarized neutron stars
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
from STGB_eqsa import STGB_eqs as STGBeqs

from EOSa import EOS

PI4 = 4.0 * pi
ka = 8.0 * pi
c = 29979245800.0  # cm/s
G = 6.67408e-8  # cm^3/g/s^2

MSUN = 1.98855e33  # g
KM = 1.0e5  # cm
mB = 1.660538921e-24  # g
E_NUCL = 2.0e14  # minimun energy density for NS core; g/cm^3

runit = 10.*KM # length to parametrize quantities
dimlmB = G*mB/c**2/runit

def rdiml(r):
    """ dimensionless length """
    return r/runit
    
def mdiml(m):
    """ dimensionless mass """
    return G*m/c**2/runit
    
def ediml(e):
    """ dimensionless energy density """
    return G*e/c**2 * runit**2
    
def pdiml(p):
    """ dimensionless pressure """
    return G*p/c**4 * runit**2

def ndiml(n):
    """ dimensionless pressure """
    return n * runit**3
    
def csqdiml(csq):
    """ dimensionless speed of sound """
    return csq/c**2 
        
def mdim(m):
    """ dimensionful mass """
    return m*c**2 * runit/G

def edim(e):
    """ dimensionful energy density """
    return c**2 * e/G/runit**2
    
def pdim(p):
    """ dimensionful pressure """
    return c**4 *p /G / runit**2

def Idim(i):
    """ dimensionful moment of inertia """
    return i*c**2 /G * runit**3
    
    
class STGB_solver(object):

    def __init__(self, EOS_name='AP4', theory=STGBeqs(xi=4.4,msq=1) ):
        """ Initialize EOS and theoretical parameters """
        self.EOS_name = EOS_name
        self.EOS = EOS(EOS_name)
        self.STGBeqs = theory


    def odesin(self, r, y):
        """ y = [mu1, bph, bphp, p, Mbar, om, omp] 
            equation arguments: (r, nu, mu1, bph, psi, p, e) 
        """
        dydr = [self.STGBeqs.dmu1dr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )), 
                self.STGBeqs.dphdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )), 
                self.STGBeqs.ddphdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )), 
                self.STGBeqs.dpdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )),
                self.STGBeqs.dMbardr( r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) ), ndiml( self.EOS.p2n( pdim(y[3]) ) ) ),
                self.STGBeqs.domdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) ), y[5], y[6] ),
                self.STGBeqs.ddomdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) ), y[5], y[6] )]
        return dydr

    def odesex(self, r, y):
        """ y = [mu1, bph, bphp, om, omp] 
            equation arguments: (r, nu, mu1, bph, psi, p, e) 
        """
        dydr = [self.STGBeqs.dmu1dr(r, 0., y[0], y[1], y[2], 0., 0.), 
                self.STGBeqs.dphdr(r, 0., y[0], y[1], y[2], 0., 0.), 
                self.STGBeqs.ddphdr(r, 0., y[0], y[1], y[2], 0., 0.),
                self.STGBeqs.domdr(r, 0., y[0], y[1], y[2], 0., 0., y[3], y[4] ),
                self.STGBeqs.ddomdr(r, 0., y[0], y[1], y[2], 0., 0., y[3], y[4] )]
        return dydr
        
        
    def bph2coef(self, ec, bphc):
        """ coefficients for the quartic equation at the center"""
        df0 = self.STGBeqs.dfxif(bphc)
        u0 = self.STGBeqs.Umf(bphc)
        du0 = self.STGBeqs.dUmf(bphc)
        a4 = 18432.*df0**3
        a3 = -768.*df0**2*(9.+2.*df0*du0)
        a2 = 288.*df0*(3.+2.*df0*du0)
        a1 = -12.*( 3. + 6.*df0*du0 + 2.*df0**2 * (2*ka*pdiml(self.EOS.e2p(ec)) + u0) * (2*ka*ediml(ec) + u0) )
        a0 = 3.*du0 + 2.*df0*( 2*ka*ediml(ec) + u0 ) * ( ka*ediml(ec) + 3.*ka*pdiml(self.EOS.e2p(ec)) - u0 )
        return [a4, a3, a2, a1, a0 ]           
        
        
    def solout(self, r, y):
        """ Stop condition for ode """
        if y[3] <= pdiml(self.EOS.min_p):
            return -1
        else:
            return 0
            
    def solout1(self, r, y):
        """ Stop condition for ode """
        if y[3] <= 1.e-10:
            return -1
        else:
            return 0 
            
    def ode_solver(self, ec, bphc, Rmax = 3.0e6, rerror=1.e-7, verbose=False):
        dimlpc = pdiml(self.EOS.e2p(ec))
        #print( self.EOS.p2n( pdim(dimlpc) )  )
        y0 = [ 0., bphc, 0., pdiml( self.EOS.e2p(ec) ), 0., 1., 0. ]
        xmax = rdiml(Rmax)
        xini = 2.e-3 # does not work for 1.e-4, 1.e-5, ...
        
        [nu2, mu12, bph2, p2, check] = self.STGBeqs.coef2( dimlpc, ediml(ec), bphc)
        e2 = p2 / csqdiml(self.EOS.p2csq( self.EOS.e2p(ec) ))
        [nu4, mu14, bph4, p4] = self.STGBeqs.coef4( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2)
        om2 = self.STGBeqs.rotcoef2( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, 1.)        
        om4 = self.STGBeqs.rotcoef4( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2, nu4, mu14, bph4, 1.)            
        """
        y0[0] = mu12*xini**2 
        y0[1] = bphc + bph2*xini**2 
        y0[2] = 2.*bph2*xini 
        y0[3] = dimlpc + p2*xini**2 
        """
        y0[0] = mu12*xini**2 + mu14*xini**4
        y0[1] = bphc + bph2*xini**2 + bph4*xini**4
        y0[2] = 2.*bph2*xini + 4*bph4*xini**3
        y0[3] = dimlpc + p2*xini**2 + p4*xini**4
        y0[4] = 4./3. *np.pi *xini**3 * dimlmB * ndiml( self.EOS.e2n( ec ) )
        y0[6] = 2.*om2*xini + 4.*om4*xini**3

        #print(y0, xini)
        solver = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
        solver.set_solout(self.solout)  # stop condition
        solver.set_initial_value(y0, xini)
        warnings.filterwarnings("ignore", category=UserWarning)
        solver.integrate( xmax )
        codein = sp_ode.get_return_code(solver)
        """
        1--successful; 2--successful and interrupted by solout; -1--input not consistent; -2--larger nstep needed; -3--step size too small; -4--stiff system
        """
        warnings.resetwarnings()

        
        R, M_s , p_s, x_s = 0, 0, 0, 0
        y_s = y0
        y_inft = [0, 0, 0, 0 ]

        y_s = solver.y
        x_s = solver.t       
        R, M_s , Mbar, p_s = runit*x_s, mdim(x_s*( 1. - np.exp(-2.*y_s[1]) )/2.), mdim(y_s[4]), pdim(y_s[4])
        #print(x_s, y_s)

 
        xinft = 20./np.sqrt(self.STGBeqs.msq)        #  cannot be small; large is fine                    
        nex = 1000
        xexset = np.linspace(x_s+(xinft-x_s)/nex, xinft, nex)
        mu1exset, bphexset, psiexset, omexset, ompexset = np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset)
        y0ex = [y_s[0], y_s[1], y_s[2], y_s[5], y_s[6]]
   
           
        solverex = sp_ode(self.odesex).set_integrator('dopri5', rtol=rerror)
        solverex.set_initial_value(y0ex, x_s)
        warnings.filterwarnings("ignore", category=UserWarning)
        for i in range(0, nex):
           solverex.integrate( xexset[i] )
           yf = solverex.y
           mu1exset[i], bphexset[i], psiexset[i], omexset[i], ompexset[i] = yf[0], yf[1], yf[2], yf[3], yf[4]
                   
        i_infmu = np.argmin( abs(mu1exset) )
        i_infbph = np.argmin( abs(bphexset) ) # suitable numerical infinity
        i_inf = i_infbph
        #print(i_infmu, i_infbph)
        dimlM_inf = xexset[i_inf]*( 1. - np.exp(-2.*mu1exset[i_inf]) )/2.
        M_inf = mdim( dimlM_inf )
        dimlJ_inf = ompexset[i_inf]*xexset[i_inf]**4 /6.
        dimlom_inf = omexset[i_inf]
        dimlI_inf = dimlJ_inf/dimlom_inf
        I_inf = Idim( dimlI_inf ) 
        #M_inf = M_inf/MSUN
        #print(dimlM_inf)

        if verbose: 
           nin = 100
           xinset = np.linspace(xini, x_s, nin)
           mu1inset, bphinset, psiinset, pinset, mbarinset, ominset, ompinset = np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset)              
           mu1inset[0], bphinset[0], psiinset[0], pinset[0], mbarinset[0], ominset[0], ompinset[0] = y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], y0[6]
              
           solverin = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
           solverin.set_solout(self.solout)  # stop condition
           solverin.set_initial_value(y0, xini)
           warnings.filterwarnings("ignore", category=UserWarning)
           for i in range(1, nin): 
              solverin.integrate( xinset[i] )
              yinf = solverin.y
              mu1inset[i], bphinset[i], psiinset[i], pinset[i], mbarinset[i], ominset[i], ompinset[i] = yinf[0], yinf[1], yinf[2], yinf[3], yinf[4], yinf[5], yinf[6]
                                                          
           f = open("stgb_5eqs_sol_data.txt", 'w+')              
           for i in range(0, nin):            
              f.write( ('%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f' % (xinset[[i]], mu1inset[i], bphinset[i], psiinset[i], pinset[i], mbarinset[i], ominset[i], ompinset[i]  ) ) + '\n' )                                    
           for i in range(0, nex): 
              f.write( ('%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f' % (xexset[[i]], mu1exset[i], bphexset[i], psiexset[i], 0., 0.,  omexset[i], ompexset[i]  ) ) + '\n' )                                                
           f.close

           f1 = open("stgb_5eqs_rot_data.txt", 'w+') 
           Jexset, Iexset = np.zeros_like(xexset), np.zeros_like(xexset)
           for i in range(0, nex):
              Jexset[i] = ompexset[i]*xexset[i]**4 /6.      
              Iexset[i] = Jexset[i]/omexset[i]       
              f1.write( ('%.9f %.9f' % (xexset[[i]], Iexset[i]  ) ) + '\n' )   
           f1.close  
              
              
                                                     
           fig = plt.figure(figsize=(18,12))
           ax1 = fig.add_subplot(241)
           ax2 = fig.add_subplot(242)
           ax3 = fig.add_subplot(243)
           ax4 = fig.add_subplot(244)
           ax5 = fig.add_subplot(245)
           ax6 = fig.add_subplot(246)
           ax7 = fig.add_subplot(247)
           ax8 = fig.add_subplot(248)
           ax1.set_xlabel(r'$x$', fontsize=20)
           ax1.set_ylabel(r'$\nu$', fontsize=20)
           ax2.set_xlabel(r'$x$', fontsize=20)
           ax2.set_ylabel(r'$\mu$', fontsize=20)
           ax3.set_xlabel(r'$x$', fontsize=20)
           ax3.set_ylabel(r'$\Phi$', fontsize=20)
           ax4.set_xlabel(r'$x$', fontsize=20)
           ax4.set_ylabel(r'$\Phi^{\prime}$', fontsize=20)
                     
           ax1.plot(xinset, mu1inset, 'o')
           ax2.plot(xinset, bphinset)
           ax3.plot(xinset, psiinset)
           ax4.plot(xinset, pinset)
           
           ax5.plot(xinset, mbarinset)          
           ax6.plot(xinset, ominset, 'o')
           ax7.plot(xexset, Jexset)
           ax8.plot(xexset, Iexset)
                           
              
                
        res = {'R': R, 'M_inf': M_inf, 'Mbar': Mbar, 'mu_s': y_s[0], 'bph_s' : y_s[1], 'psi_s' : y_s[2], 'dimp_s': y_s[4], 'codein': codein, 'r_inf': xexset[i_inf], 'mu_inf': mu1exset[i_inf], 'bph_inf' : bphexset[i_inf], 'psi_inf' : psiexset[i_inf], 'dimlJ_inf':dimlJ_inf, 'dimlom_inf':dimlom_inf, 'I':I_inf, 'I_sph': 0.4*M_inf*R**2}  
                     
        return res

        
    def ode_solsearchprep(self, ec, verbose=False):        # preliminary search

        bphcmin = 0.001
        bphcstep = 0.001  # for xi < 0
        """
        bphcmin = 0.00001
        bphcstep = 0.00001 # for xi > 0 
        """        
        bphctest = bphcmin               
        conroots = np.roots( self.bph2coef(ec, bphctest) )
        sol = self.ode_solverprep(ec, bphctest)

        #f = open("stgb_solver_search_data1.txt", 'w+')        
        codein = sol['codein']  
        bph_inft = sol['bph_inft'] 
        bph_inft1 = bph_inft         
        print( bphctest, conroots[3], bph_inft )
        
        signcount = 0
        bphclower, bphcupper = 0., 0.
        while codein==2 or codein==-2 and np.imag(conroots[3]) == 0.0 and bphctest<1. :
           bphctest = bphctest + bphcstep
           conroots = np.roots( self.bph2coef(ec, bphctest) )
           sol = self.ode_solverprep(ec, bphctest)
           codein = sol['codein']  
           bph_inft = sol['bph_inft'] 
           print( bphctest, conroots[3], codein, bph_inft )
           #f.write( ('%.6f %.6f' % (bphctest, bph_inft ) ) + '\n' ) 
           
           if bph_inft * bph_inft1 < 0.0:
              signcount = signcount + 1
              bphclower, bphcupper = bphctest-bphcstep, bphctest
              #f.write( ('%.6f, %.6f' % (bphctest-bphcstep, bphctest ) ) + '\n' ) 
           bph_inft1 = bph_inft
           
        #f.close()
           
        return {'ec': ec, 'signcount': signcount, 'codein': codein, 'bphclower': bphclower, 'bphcupper': bphcupper}
                
        
    def ode_solsearch(self, ec, bphcl, bphch, verbose=False):       # search for high accuracy; bisection method fails for accuracy at 1e-7 
        bphcmin = bphcl
        bphcstep = 1.e-5 # for xi < 0 
        """
        bphcmin = bphcl
        bphcstep = 1.e-7 # for xi > 0 
        """
        bphctest = bphcmin                
        conroots = np.roots( self.bph2coef(ec, bphctest) )
        sol = self.ode_solverprep(ec, bphctest)

        bph_inft = sol['bph_inft'] 
        bph_inft1 = bph_inft         
        print( bphctest, conroots[3], bph_inft )
        
        signcount = 0
        bphclower, bphcupper = 0., 0.
        while np.imag(conroots[3]) == 0.0 and bphctest<bphch :
           bphctest = bphctest + bphcstep
           conroots = np.roots( self.bph2coef(ec, bphctest) )
           sol = self.ode_solverprep(ec, bphctest)
           
           bph_inft = sol['bph_inft'] 
           print( bphctest, conroots[3], bph_inft )
           #f.write( ('%.6f %.6f' % (bphctest, bph_inft ) ) + '\n' ) 
           
           if bph_inft * bph_inft1 < 0.0:
              signcount = signcount + 1
              bphclower, bphcupper = bphctest-bphcstep, bphctest
              break
               
           bph_inft1 = bph_inft
           
           
        return {'ec': ec, 'signcount': signcount, 'bphclower': bphclower, 'bphcupper': bphcupper}
    
    
    def bisec_bphc(self, ec, bphclower, bphcupper):
        bphctest = bphclower
        conroots = np.roots( self.bph2coef(ec, bphctest) )
        sol = self.ode_solverprep(ec, bphctest)
         
        bph_inft = sol['bph_inft'] 
        bph_inft1 = bph_inft 
          
        niter = 0      
        while abs(bphctest - (bphclower + bphcupper)/2.) > 1.e-7 and niter < 100:
           bphctest = (bphclower + bphcupper)/2.
           conroots = np.roots( self.bph2coef(ec, bphctest) )
           sol = self.ode_solverprep(ec, bphctest)
           bph_inft = sol['bph_inft'] 
           
           #print(bphclower, bphcupper, bph_inft )
           
           if bph_inft * bph_inft1 > 0.0:
              bphclower = bphctest
           else:
              bphcupper = bphctest 
              
           bph_inft1 = bph_inft  
           niter = niter + 1  
           
        return (bphclower + bphcupper)/2.
    
        
    def ode_solec(self, ec0, ec1):        
        nec = 100
        ecset = np.linspace(ec0, ec1, nec)

        f = open("stgb_solver_search_data1.txt", 'w+')        
        for i in range(0, nec):
           bphcsol = self.ode_solsearch(ecset[i])
           f.write( ('%.9e %.9f %.9f' % (ecset[i], bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' ) 
           
           if bphcsol['bphcupper'] == 0:
              break
                      
        f.close()
        
        return {i, ecset[i]}
        
        
        
        
        
        
t0 = timeit.time.time()

eos_nameset = ['WFF1','SLy4','AP4','MPA1','PAL1']
xiset = [3., 5., 10., -0.2, -0.5, -1.]
msqset = [0.01**2, 0.1**2, 1.**2]

ieos = 2
msq = 1.**2.
xi = -1.
eos_name = eos_nameset[ieos]
zz = STGB_solver(EOS_name=eos_name, theory=STGBeqs(xi, msq) )

switch = 1

if switch: # single run for testing use
  ec = 2.319000E+15
  bphc = 0.01388
  zzz = zz.ode_solver(ec, bphc)
  print( zzz )
  #ec = 3.8377045727049740e+15
  #bphc1 = 0.01389  
  #zzz1 = zz.ode_solver(ec, bphc1)
  #print(('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (ec, bphc, zzz['R'], zzz['M_inf']/MSUN, zzz['Mbar']/MSUN, zzz['I'], zzz['dimlJ_inf'], zzz['dimlom_inf'] ) ))  
  #print(('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (ec, bphc1, zzz1['R'], zzz1['M_inf']/MSUN, zzz1['Mbar']/MSUN, zzz1['I'], zzz1['dimlJ_inf'], zzz1['dimlom_inf'] ) ))  
  
  
  
if not switch:  # multiple runs for generating the data
  f = open("stgb_rot_data1.txt", 'w+') 
  f1 = open("stgb_rot_error1.txt", 'w+')   
  seadata = np.genfromtxt('stgb_solver_pdata1.txt')
  if seadata.size !=0:
    c0, c1, c2, c3, c4 = seadata[:, 0], seadata[:, 1], seadata[:, 2], seadata[:, 3], seadata[:, 4]
    nec = len(seadata)       
    for i in range(0, nec):
       zzz = zz.ode_solver(c0[i], c1[i])
       zzz1 = zz.ode_solver(c0[i], c2[i])
       f.write( ('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (c0[i], c1[i], c2[i], c3[i], (zzz1['R']+zzz['R'])/2.-c3[i], c4[i], (zzz1['M_inf']+zzz['M_inf'])/2.-c4[i], (zzz1['Mbar']+zzz['Mbar'])/2., (zzz1['I']+zzz['I'])/2., (zzz1['dimlJ_inf']+zzz['dimlJ_inf'])/2., (zzz1['dimlom_inf']+zzz['dimlom_inf'])/2. ) ) + '\n' ) 
       f1.write( ('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (c0[i], c1[i], c2[i], (zzz1['R']-zzz['R'])/2., (zzz1['M_inf']-zzz['M_inf'])/2., (zzz1['Mbar']-zzz['Mbar'])/2., (zzz1['I']-zzz['I'])/2., (zzz1['dimlJ_inf']-zzz['dimlJ_inf'])/2., (zzz1['dimlom_inf']-zzz['dimlom_inf'])/2. ) ) + '\n' )       
                        
  f.close()
  f1.close()




print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))


plt.show()
