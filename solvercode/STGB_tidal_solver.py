"""
solve the tidal deformability of scalarized neutron stars
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
from scipy import interpolate
from scipy.optimize import curve_fit

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

def rdiml(r):
    """ dimensionless length """
    return r/runit
    
def mdiml(m):
    """ dimensionless energy density """
    return G*m/c**2/runit
    
def ediml(e):
    """ dimensionless energy density """
    return G*e/c**2 * runit**2
    
def pdiml(p):
    """ dimensionless pressure """
    return G*p/c**4 * runit**2
    
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

class STGB_solver(object):

    def __init__(self, EOS_name='AP4', theory=STGBeqs(xi=4.4,msq=1) ):
        """ Initialize EOS and theoretical parameters """
        self.EOS_name = EOS_name
        self.EOS = EOS(EOS_name)
        self.STGBeqs = theory


    def odesin(self, r, y):
        """ y = [ mu1, bph, bph', p, H_0, K, debph, debph']  
        """
        dydr = [self.STGBeqs.dmu1dr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )), 
                self.STGBeqs.dphdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )), 
                self.STGBeqs.ddphdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )), 
                self.STGBeqs.dpdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) )),
                self.STGBeqs.deh0dr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) ), csqdiml(self.EOS.p2csq( pdim(y[3]) )), y[4], y[5], y[6], y[7]),
                self.STGBeqs.dekdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) ), csqdiml(self.EOS.p2csq( pdim(y[3]) )), y[4], y[5], y[6], y[7]),
                self.STGBeqs.ddebphdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) ), csqdiml(self.EOS.p2csq( pdim(y[3]) )), y[4], y[5], y[6], y[7]),
                self.STGBeqs.dddebphdr(r, 0., y[0], y[1], y[2], y[3], ediml( self.EOS.p2e( pdim(y[3]) ) ), csqdiml(self.EOS.p2csq( pdim(y[3]) )), y[4], y[5], y[6], y[7])]
        return dydr

    def odesex(self, r, y):
        """ y = [ mu1, bph, bph', H_0, K, debph, debph'] 
        """
        dydr = [self.STGBeqs.dmu1dr(r, 0., y[0], y[1], y[2], 0., 0.), 
                self.STGBeqs.dphdr(r, 0., y[0], y[1], y[2], 0., 0.), 
                self.STGBeqs.ddphdr(r, 0., y[0], y[1], y[2], 0., 0.), 
                self.STGBeqs.deh0dr(r, 0., y[0], y[1], y[2], 0., 0., 1., y[3], y[4], y[5], y[6]),
                self.STGBeqs.dekdr(r, 0., y[0], y[1], y[2], 0., 0., 1., y[3], y[4], y[5], y[6]),
                self.STGBeqs.ddebphdr(r, 0., y[0], y[1], y[2], 0., 0., 1., y[3], y[4], y[5], y[6]),
                self.STGBeqs.dddebphdr(r, 0., y[0], y[1], y[2], 0., 0., 1., y[3], y[4], y[5], y[6])]
        return dydr
        
  
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

    def odein_solver(self, eh00, debph0, ec, bphc, Rmax = 3.e6, rerror=1.e-7, verbose=False): # rerror does not work at 1e-8; 10 mins with 100-point plot
        dimlpc = pdiml(self.EOS.e2p(ec))
        csq0 = csqdiml(self.EOS.p2csq( self.EOS.e2p(ec) ))
        y0sol1 = [ 0., bphc, 0., pdiml( self.EOS.e2p(ec) ), 0, 0, 0, 0 ]
        xmax = rdiml(Rmax)
        xini = 2.e-3 # does not work for 1.e-4, 1.e-5, ...       
        [nu2, mu12, bph2, p2, check] = self.STGBeqs.coef2( dimlpc, ediml(ec), bphc)
        e2 = p2 / csq0
        [nu4, mu14, bph4, p4] = self.STGBeqs.coef4( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2)
        
        #eh00, debph0 = 1., 0.
        ek0 = (1-16.*self.STGBeqs.xi*bphc*bph2)*eh00 + 32*self.STGBeqs.xi*nu2*bphc*debph0 
        [eh02, debph2, ek2] = self.STGBeqs.tidcoef2( dimlpc, ediml(ec), csq0, bphc, nu2, mu12, bph2, nu4, mu14, bph4, eh00, debph0 )
        
        #print(nu2, mu12, bph2, p2)
        #print(nu4, mu14, bph4, p4)
        #print(ek0)
        #print(eh02, debph2, ek2)
        """
        y0sol1[2] = 2.*bph2*xini  
        y0sol1[4] = xini**2 * eh00       
        y0sol1[5] = xini**2 * ek0          
        y0sol1[6] = xini**2 * debph0          
        y0sol1[7] = 2.*debph0*xini        
        """
        y0sol1[0] = mu12*xini**2 + mu14*xini**4
        y0sol1[1] = bphc + bph2*xini**2 + bph4*xini**4
        y0sol1[2] = 2.*bph2*xini + 4*bph4*xini**3
        y0sol1[3] = dimlpc + p2*xini**2 + p4*xini**4
        y0sol1[4] = xini**2 * (eh00 + eh02*xini**2)       
        y0sol1[5] = xini**2 * (ek0 + ek2*xini**2)         
        y0sol1[6] = xini**2 * (debph0 + debph2*xini**2)         
        y0sol1[7] = 2.*debph0*xini + 4.*debph2*xini**3         

        print( y0sol1 )         
        solver = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
        solver.set_solout(self.solout)  # stop condition
        solver.set_initial_value(y0sol1, xini)
        warnings.filterwarnings("ignore", category=UserWarning)
        solver.integrate( xmax )
        codein = sp_ode.get_return_code(solver)
        """
        1--successful; 2--successful and interrupted by solout; -1--input not consistent; -2--larger nstep needed; -3--step size too small; -4--stiff system
        """
        warnings.resetwarnings()
        y_s = solver.y      
        x_s = solver.t   
        print(codein)
        print(x_s, y_s)
         
        res = { 'codein': codein }  
        return res
        
    def ode_solver(self, eh00, debph0, ec, bphc, Rmax = 3.e6, rerror=1.e-7, verbose=False): # rerror does not work at 1e-8; 10 mins with 100-point plot
        dimlpc = pdiml(self.EOS.e2p(ec))
        csq0 = csqdiml(self.EOS.p2csq( self.EOS.e2p(ec) ))
        y0sol1 = [ 0., bphc, 0., pdiml( self.EOS.e2p(ec) ), 0, 0, 0, 0 ]
        xmax = rdiml(Rmax)
        xini = 2.e-3 # does not work for 1.e-4, 1.e-5, ...       
        [nu2, mu12, bph2, p2, check] = self.STGBeqs.coef2( dimlpc, ediml(ec), bphc)
        e2 = p2 / csq0
        [nu4, mu14, bph4, p4] = self.STGBeqs.coef4( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2)
        
        #eh00, debph0 = 1., 0.
        ek0 = (1-16.*self.STGBeqs.xi*bphc*bph2)*eh00 + 32*self.STGBeqs.xi*nu2*bphc*debph0 
        [eh02, debph2, ek2] = self.STGBeqs.tidcoef2( dimlpc, ediml(ec), csq0, bphc, nu2, mu12, bph2, nu4, mu14, bph4, eh00, debph0 )
        
        #print(nu2, mu12, bph2, p2)
        #print(nu4, mu14, bph4, p4)
        #print(ek0)
        #print(eh02, debph2, ek2)
        """
        y0sol1[2] = 2.*bph2*xini  
        y0sol1[4] = xini**2 * eh00       
        y0sol1[5] = xini**2 * ek0          
        y0sol1[6] = xini**2 * debph0          
        y0sol1[7] = 2.*debph0*xini        
        """
        y0sol1[0] = mu12*xini**2 + mu14*xini**4
        y0sol1[1] = bphc + bph2*xini**2 + bph4*xini**4
        y0sol1[2] = 2.*bph2*xini + 4*bph4*xini**3
        y0sol1[3] = dimlpc + p2*xini**2 + p4*xini**4
        y0sol1[4] = xini**2 * (eh00 + eh02*xini**2)       
        y0sol1[5] = xini**2 * (ek0 + ek2*xini**2)         
        y0sol1[6] = xini**2 * (debph0 + debph2*xini**2)         
        y0sol1[7] = 2.*debph0*xini + 4.*debph2*xini**3         

        print( y0sol1 )         
        solver = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
        solver.set_solout(self.solout)  # stop condition
        solver.set_initial_value(y0sol1, xini)
        warnings.filterwarnings("ignore", category=UserWarning)
        solver.integrate( xmax )
        codein = sp_ode.get_return_code(solver)
        """
        1--successful; 2--successful and interrupted by solout; -1--input not consistent; -2--larger nstep needed; -3--step size too small; -4--stiff system
        """
        warnings.resetwarnings()
        y_s = solver.y      
        x_s = solver.t   
        print(codein)
        print(x_s, y_s)
 
        nex = 100        
        xinft = 20./np.sqrt(self.STGBeqs.msq) 
        xexset = np.linspace(x_s + (xinft-x_s)/nex, xinft, nex)
        mu1exset1, bphexset1, psiexset1, eh0exset1, ekexset1, debphexset1, depsiexset1 = np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset)               
        yex0sol = [y_s[0], y_s[1], y_s[2], y_s[4], y_s[5], y_s[6], y_s[7]]
        solverex = sp_ode(self.odesex).set_integrator('dopri5', rtol=rerror) 
        solverex.set_initial_value(yex0sol, x_s)
        warnings.filterwarnings("ignore", category=UserWarning)
        for i in range(0, nex):
           solverex.integrate( xexset[i] )
           yexsol = solverex.y
           mu1exset1[i], bphexset1[i], psiexset1[i], eh0exset1[i], ekexset1[i], debphexset1[i], depsiexset1[i] = yexsol[0], yexsol[1], yexsol[2], yexsol[3], yexsol[4], yexsol[5], yexsol[6]

        i_inf = np.argmin( abs(bphexset1) )
        dimlM_inf = 0.5* xexset[i_inf]*( 1.0-np.exp(-2.* mu1exset1[i_inf]) )                      
        M_inf = mdim(dimlM_inf)
        

        if verbose:
          nin = 100
          xinset = np.linspace(xini, x_s, nin)
          mu1inset1, bphinset1, psiinset1, pinset1, eh0inset1, ekinset1, debphinset1, depsiinset1 = np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset)
          mu1inset1[0], bphinset1[0], psiinset1[0], pinset1[0], eh0inset1[0], ekinset1[0], debphinset1[0], depsiinset1[0] = y0sol1[0], y0sol1[1], y0sol1[2], y0sol1[3], y0sol1[4], y0sol1[5], y0sol1[6], y0sol1[7]                          
          solverin = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
          solverin.set_initial_value(y0sol1, xini)
          warnings.filterwarnings("ignore", category=UserWarning)
          for i in range(1, nin):
             solverin.integrate( xinset[i] )
             yin = solverin.y
             mu1inset1[i], bphinset1[i], psiinset1[i], pinset1[i] = yin[0], yin[1], yin[2], yin[3]
             eh0inset1[i], ekinset1[i], debphinset1[i], depsiinset1[i] = yin[4], yin[5], yin[6], yin[7]

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
          ax1.set_ylabel(r'$\mu$', fontsize=20)
          ax2.set_xlabel(r'$x$', fontsize=20)
          ax2.set_ylabel(r'$\Phi$', fontsize=20)
          ax3.set_xlabel(r'$x$', fontsize=20)
          ax3.set_ylabel(r'$H_0$', fontsize=20)
          ax4.set_xlabel(r'$x$', fontsize=20)
          ax4.set_ylabel(r'$\delta \Phi$', fontsize=20)


          ax1.plot(xinset, mu1inset1) 
          ax1.plot(xexset, mu1exset1)       
          ax5.plot(xinset, bphinset1) 
          ax5.plot(xexset, bphexset1) 
                
          ax2.plot(xinset, eh0inset1)              
          ax6.plot(xexset, eh0exset1)

          ax3.plot(xinset, ekinset1)              
          ax7.plot(xexset, ekexset1)
                  
          ax4.plot(xinset, debphinset1) 
          ax8.plot(xexset, debphexset1)          
                 
          ax8.plot(xinset, depsiinset1)               
          ax8.plot(xexset, depsiexset1)          
          
                  
        res = {'mu_s': y_s[0], 'bph_s': y_s[1], 'psi_s': y_s[2], 'dimp_s': y_s[3], 'eh0_s': y_s[4], 'ek_s': y_s[5], 'debph_s': y_s[6], 'depsi_s': y_s[7], 'mu_inf': mu1exset1[i_inf], 'bph_inf': bphexset1[i_inf], 'psi_inf': psiexset1[i_inf], 'eh0_inf': eh0exset1[i_inf], 'ek_inf': ekexset1[i_inf], 'debph_inf': debphexset1[i_inf], 'depsi_inf': depsiexset1[i_inf], 'r_inf': xexset[i_inf], 'M_inf': M_inf, 'codein': codein  }
        return res

        
    def tidal_solver(self, ec, bphc, Rmax = 3.e6, rerror=1.e-7, verbose=False):        
        dimlpc = pdiml(self.EOS.e2p(ec))
        csq0 = csqdiml(self.EOS.p2csq( self.EOS.e2p(ec) ))
        y0sol1 = [ 0., bphc, 0., pdiml( self.EOS.e2p(ec) ), 0, 0, 0, 0 ]
        xmax = rdiml(Rmax)
        xini = 2.e-3 # does not work for 1.e-4, 1.e-5, ...       
        [nu2, mu12, bph2, p2, check] = self.STGBeqs.coef2( dimlpc, ediml(ec), bphc)
        e2 = p2 / csq0
        [nu4, mu14, bph4, p4] = self.STGBeqs.coef4( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2)

        # 1st solution        
        eh00, debph0 = 1., 0.
        """
        1: (1,0)(1,0.4)      2: (1,0.01)(1,0.4) 
        3: (0.1,0)(0.1,0.4)  4:(0.1,0.01)(0.1,0.4) 
        5:(0.1,0)(0.1,0.1)   6:(0.1,0.01)(0.1,0.1)
        1, 3, 5 better than 2, 4, 6 
        """
        ek0 = (1-16.*self.STGBeqs.xi*bphc*bph2)*eh00 + 32*self.STGBeqs.xi*nu2*bphc*debph0 
        [eh02, debph2, ek2] = self.STGBeqs.tidcoef2( dimlpc, ediml(ec), csq0, bphc, nu2, mu12, bph2, nu4, mu14, bph4, eh00, debph0 )
        
        #print(nu2, mu12, bph2, p2)
        #print(nu4, mu14, bph4, p4)
        #print(ek0)
        #print(eh02, debph2, ek2)
        """
        y0sol1[2] = 2.*bph2*xini  
        y0sol1[4] = xini**2 * eh00       
        y0sol1[5] = xini**2 * ek0          
        y0sol1[6] = xini**2 * debph0          
        y0sol1[7] = 2.*debph0*xini        
        """
        y0sol1[0] = mu12*xini**2 + mu14*xini**4
        y0sol1[1] = bphc + bph2*xini**2 + bph4*xini**4
        y0sol1[2] = 2.*bph2*xini + 4*bph4*xini**3
        y0sol1[3] = dimlpc + p2*xini**2 + p4*xini**4
        y0sol1[4] = xini**2 * (eh00 + eh02*xini**2)       
        y0sol1[5] = xini**2 * (ek0 + ek2*xini**2)         
        y0sol1[6] = xini**2 * (debph0 + debph2*xini**2)         
        y0sol1[7] = 2.*debph0*xini + 4.*debph2*xini**3         

        #print( y0sol1 )         
        solver = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
        solver.set_solout(self.solout)  # stop condition
        solver.set_initial_value(y0sol1, xini)
        warnings.filterwarnings("ignore", category=UserWarning)
        solver.integrate( xmax )
        codein1 = sp_ode.get_return_code(solver)
        """
        1--successful; 2--successful and interrupted by solout; -1--input not consistent; -2--larger nstep needed; -3--step size too small; -4--stiff system
        """
        warnings.resetwarnings()
        y1_s = solver.y
        x1_s = solver.t 
        R1 = x1_s *runit                 
        print(codein1)
        print(x1_s, y1_s)        

        # 2nd solution
        eh00, debph0 = 1., 1.
        ek0 = (1-16.*self.STGBeqs.xi*bphc*bph2)*eh00 + 32*self.STGBeqs.xi*nu2*bphc*debph0 
        [eh02, debph2, ek2] = self.STGBeqs.tidcoef2( dimlpc, ediml(ec), csq0, bphc, nu2, mu12, bph2, nu4, mu14, bph4, eh00, debph0 )
        y0sol1[0] = mu12*xini**2 + mu14*xini**4
        y0sol1[1] = bphc + bph2*xini**2 + bph4*xini**4
        y0sol1[2] = 2.*bph2*xini + 4*bph4*xini**3
        y0sol1[3] = dimlpc + p2*xini**2 + p4*xini**4
        y0sol1[4] = xini**2 * (eh00 + eh02*xini**2)       
        y0sol1[5] = xini**2 * (ek0 + ek2*xini**2)         
        y0sol1[6] = xini**2 * (debph0 + debph2*xini**2)         
        y0sol1[7] = 2.*debph0*xini + 4.*debph2*xini**3         

        #print( y0sol1 )         
        solver = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
        solver.set_solout(self.solout)  # stop condition
        solver.set_initial_value(y0sol1, xini)
        warnings.filterwarnings("ignore", category=UserWarning)
        solver.integrate( xmax )
        codein2 = sp_ode.get_return_code(solver)
        warnings.resetwarnings()
        y2_s = solver.y
        x2_s = solver.t
        R2 = x2_s *runit         
        print(codein2)
        print(x2_s, y2_s)

        codein = min(codein1, codein2)            
        res = {'codein': codein, 'mu_inf': 0, 'bph_inf': 0, 'psi_inf': 0, 'eh0_inf': 0, 'ek_inf': 0, 'debph_inf': 0, 'depsi_inf': 0, 'r_inf': 0, 'M_inf': 0, 'eh0p_inf': 0, 'ytidal': 0, 'dimllatidal': 0, 'dimllatidal1': 0, 'dimllatidal2': 0, 'dimllatidal3': 0, 'dimllatidal4': 0, 'sigma1': 0, 'sigma2': 0, 'sigma3': 0, 'sigma4': 0, 'R': 0, 'deR': 0 }
              
        if min(x1_s, x2_s) > xini*400:
          if np.sqrt(self.STGBeqs.msq) >= 1.: 
            nex = 1000
          else:
            nex = 1000           
          
          x_s, y_s = x1_s, y1_s  
          xinft = 10./np.sqrt(self.STGBeqs.msq)          
          xexset = np.linspace(x_s + (xinft-x_s)/nex, xinft, nex)
          mu1exset1, bphexset1, psiexset1, eh0exset1, ekexset1, debphexset1, depsiexset1, mu1exset2, bphexset2, psiexset2, eh0exset2, ekexset2, debphexset2, depsiexset2 = np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset)
          yex0sol = [y_s[0], y_s[1], y_s[2], y_s[4], y_s[5], y_s[6], y_s[7]]
          solverex = sp_ode(self.odesex).set_integrator('dopri5', rtol=rerror) 
          solverex.set_initial_value(yex0sol, x_s)
          warnings.filterwarnings("ignore", category=UserWarning)
          for i in range(0, nex):
             solverex.integrate( xexset[i] )
             yexsol = solverex.y
             mu1exset1[i], bphexset1[i], psiexset1[i], eh0exset1[i], ekexset1[i], debphexset1[i], depsiexset1[i] = yexsol[0], yexsol[1], yexsol[2], yexsol[3], yexsol[4], yexsol[5], yexsol[6]

          y_s = y2_s
          yex0sol = [y_s[0], y_s[1], y_s[2], y_s[4], y_s[5], y_s[6], y_s[7]]
          solverex = sp_ode(self.odesex).set_integrator('dopri5', rtol=rerror) 
          solverex.set_initial_value(yex0sol, x_s)
          warnings.filterwarnings("ignore", category=UserWarning)
          for i in range(0, nex):
             solverex.integrate( xexset[i] )
             yexsol = solverex.y
             mu1exset2[i], bphexset2[i], psiexset2[i], eh0exset2[i], ekexset2[i], debphexset2[i], depsiexset2[i] = yexsol[0], yexsol[1], yexsol[2], yexsol[3], yexsol[4], yexsol[5], yexsol[6]

          i_inf1 = np.argmin( abs(bphexset1) )
          i_inf2 = np.argmin( abs(bphexset2) )
          i_inf = min(i_inf1, i_inf2)
          dimlM_inf = 0.5* xexset[i_inf]*( 1.0-np.exp(-2.* mu1exset1[i_inf]) )
          M_inf = mdim(dimlM_inf)
          print(i_inf1, i_inf2, xexset[i_inf], dimlM_inf )
                         
          y1oy2 = debphexset1[i_inf]/debphexset2[i_inf]       
          eh0exset = eh0exset1 - y1oy2*eh0exset2 
          ekexset = ekexset1 - y1oy2*ekexset2
          debphexset = debphexset1 - y1oy2*debphexset2
          depsiexset = depsiexset1 - y1oy2*depsiexset2  
                 
          if verbose:
            f = open("stgb_tidal_sol_v1_data1.txt", 'w+')        
            f1 = open("stgb_tidal_sol1_v1_data1.txt", 'w+')
            f2 = open("stgb_tidal_sol2_v1_data1.txt", 'w+')
            for i in range (0, nex):
               f1.write( ('%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f' % ( xexset[i], mu1exset1[i], bphexset1[i], psiexset1[i], eh0exset1[i], ekexset1[i], debphexset1[i], depsiexset1[i] ) ) + '\n' )               
               f2.write( ('%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f' % ( xexset[i], mu1exset2[i], bphexset2[i], psiexset2[i], eh0exset2[i], ekexset2[i], debphexset2[i], depsiexset2[i] ) ) + '\n' )
              
               utidal = dimlM_inf/xexset[i]
               eh0p = self.STGBeqs.deh0dr(xexset[i], 0., mu1exset1[i], bphexset1[i], psiexset1[i], 0., 0., 1., eh0exset[i], ekexset[i], debphexset[i], depsiexset[i] )
               ytidal = xexset[i]*eh0p/eh0exset[i]
               dimllatidal = 16./15. * (1.-2.*utidal)**2 *(2.-ytidal+2.*utidal*(ytidal-1)) / ( 8.*utidal**5*(ytidal+1) + 4.*utidal**4*(3.*ytidal-2.) + 4.*utidal**3*(13.-11.*ytidal) + 6.*utidal**2 *(5.*ytidal-8.) + 2.*utidal*(6.-3.*ytidal) +3.*(1.-2.*utidal)**2*(2.-ytidal+2.*utidal*(ytidal-1))*np.log(1.-2.*utidal) )   
               f.write( ('%.9f %.9f %.9f %.9f %.9f %.9f %.9f' % ( xexset[i], eh0exset[i], ekexset[i], debphexset[i], depsiexset[i], ytidal, dimllatidal ) ) + '\n' )                                                
            f.close()  
            f1.close()
            f2.close()

          # calculate tidal deformability using Legendre function
          utidal = dimlM_inf/xexset[i_inf]
          eh0p = self.STGBeqs.deh0dr(xexset[i_inf], 0., mu1exset1[i_inf], bphexset1[i_inf], psiexset1[i_inf], 0., 0., 1., eh0exset[i_inf], ekexset[i_inf], debphexset[i_inf], depsiexset[i_inf] )
          ytidal = xexset[i_inf]*eh0p/eh0exset[i_inf]
          dimllatidal = 16./15. * (1.-2.*utidal)**2 *(2.-ytidal+2.*utidal*(ytidal-1)) / ( 8.*utidal**5*(ytidal+1) + 4.*utidal**4*(3.*ytidal-2.) + 4.*utidal**3*(13.-11.*ytidal) + 6.*utidal**2 *(5.*ytidal-8.) + 2.*utidal*(6.-3.*ytidal) +3.*(1.-2.*utidal)**2*(2.-ytidal+2.*utidal*(ytidal-1))*np.log(1.-2.*utidal) )         

          # calculate tidal deformability by fitting
          def asyh0(x, a2, a1, am3):
             return a2*x**2+a1*x+am3/x**3

          def removediag(amatrix):
             bmatrix = np.zeros_like(amatrix)
             for i in range(0, bmatrix.shape[0]):
                for j in range(0, bmatrix.shape[1]):
                   if i!=j:
                     bmatrix[i, j] = amatrix[i,j]
             return bmatrix 
                
          intvwid = int( i_inf/5 )          
          nfit4 = i_inf - intvwid
          nfit3 = i_inf - 2*intvwid
          nfit2 = i_inf - 3*intvwid
          nfit1 = i_inf - 4*intvwid 
                       
          xdata = xexset[nfit1:nfit2]   
          ydata = eh0exset[nfit1:nfit2]   
          popt, pcov = curve_fit(asyh0, xdata, ydata)
          print(popt, pcov)
          dimllatidal1 = popt[2]/popt[0]/3./dimlM_inf**5.
          std = np.sqrt(np.diag(pcov))
          pcovmod=removediag(pcov)          
          relerror1, correl1 = abs(std[0]/popt[0])+abs(std[2]/popt[2]), np.max(abs(pcovmod))          

          xdata = xexset[nfit2:nfit3]   
          ydata = eh0exset[nfit2:nfit3]   
          popt, pcov = curve_fit(asyh0, xdata, ydata)
          print(popt, pcov)
          dimllatidal2 = popt[2]/popt[0]/3./dimlM_inf**5.
          std = np.sqrt(np.diag(pcov))
          pcovmod=removediag(pcov)          
          relerror2, correl2 = abs(std[0]/popt[0])+abs(std[2]/popt[2]), np.max(abs(pcovmod))

          xdata = xexset[nfit3:nfit4]   
          ydata = eh0exset[nfit3:nfit4]   
          popt, pcov = curve_fit(asyh0, xdata, ydata)
          print(popt, pcov)
          dimllatidal3 = popt[2]/popt[0]/3./dimlM_inf**5.
          std = np.sqrt(np.diag(pcov))
          pcovmod=removediag(pcov)          
          relerror3, correl3 = abs(std[0]/popt[0])+abs(std[2]/popt[2]), np.max(abs(pcovmod)) 

          xdata = xexset[nfit4:i_inf]   
          ydata = eh0exset[nfit4:i_inf]   
          popt, pcov = curve_fit(asyh0, xdata, ydata)
          print(popt, pcov)
          dimllatidal4 = popt[2]/popt[0]/3./dimlM_inf**5.
          std = np.sqrt(np.diag(pcov))
          pcovmod=removediag(pcov)          
          relerror4, correl4 = abs(std[0]/popt[0])+abs(std[2]/popt[2]), np.max(abs(pcovmod))                    
         
          res = {'codein': codein, 'mu_inf': mu1exset1[i_inf], 'bph_inf': bphexset1[i_inf], 'psi_inf': psiexset1[i_inf], 'eh0_inf': eh0exset[i_inf], 'ek_inf': ekexset[i_inf], 'debph_inf': debphexset[i_inf], 'depsi_inf': depsiexset[i_inf], 'r_inf': xexset[i_inf], 'M_inf': M_inf, 'eh0p_inf': eh0p, 'ytidal': ytidal, 'dimllatidal': dimllatidal, 'dimllatidal1': dimllatidal1, 'dimllatidal2': dimllatidal2, 'dimllatidal3': dimllatidal3, 'dimllatidal4': dimllatidal4, 'relerror1': relerror1, 'relerror2': relerror2, 'relerror3': relerror3, 'relerror4': relerror4, 'correl1': correl1, 'correl2': correl2, 'correl3': correl3, 'correl4': correl4, 'R': R1, 'deR': R2-R1 }
        return res        
        
        
        
        
        
        
        
t0 = timeit.time.time()

eos_nameset = ['WFF1','SLy4','AP4','MPA1','PAL1']
xiset = [3., 5., 10., -0.2, -0.5, -1.]
msqset = [0.01**2, 0.1**2, 1.**2]

ieos = 2
msq = 1.**2
xi = 3.
eos_name = eos_nameset[ieos]

zz = STGB_solver(EOS_name=eos_name, theory=STGBeqs(xi, msq) )

switch = 1
sswitch = 1

if switch:
  if sswitch: # single run for testing use
    ec = 1.6859273485299200e+15 
    bphc = 2.682000000e-04
    zzz = zz.odein_solver(1., 0.4, ec, bphc)
    #zzz = zz.tidal_solver(ec, bphc)
    #print( zzz )
    #bphc1 = 2.683000000e-04
    #zzz1 = zz.tidal_solver(ec, bphc1)    
    #print( (zzz['dimllatidal1']+zzz1['dimllatidal1'])/2., (zzz1['dimllatidal1']-zzz['dimllatidal1'])/2., (zzz['dimllatidal2']+zzz1['dimllatidal2'])/2., (zzz1['dimllatidal2']-zzz['dimllatidal2'])/2., (zzz['dimllatidal3']+zzz1['dimllatidal3'])/2., (zzz1['dimllatidal3']-zzz['dimllatidal3'])/2., (zzz['dimllatidal4']+zzz1['dimllatidal4'])/2., (zzz1['dimllatidal4']-zzz['dimllatidal4'])/2., (zzz['dimllatidal']+zzz1['dimllatidal'])/2., (zzz1['dimllatidal']-zzz['dimllatidal'])/2. )
    

  if not sswitch: # single run for testing use 
    tidprepdata = np.genfromtxt('stgb_tidprep_data31.txt' )
    if tidprepdata.size !=0:
      linech = 13
      ec, bphc = tidprepdata[linech-1, 0], tidprepdata[linech-1, 1]
      print(ec, bphc)
      zzz = zz.tidal_solver( ec, bphc, verbose = 0)
      print(zzz)    




if not switch:
  if sswitch: # unnecessary part for testing use
    tidprepdata = np.genfromtxt('stgb_tidprepinput_data1.txt' )
    f = open('stgb_tidprep_v1_data1.txt', 'w+')
    if tidprepdata.size !=0:
      istep = int(len(tidprepdata) /15)+1
      print(istep)
      i = 0
      while i < len(tidprepdata):
         ec, bphc = tidprepdata[i, 0], tidprepdata[i, 1]
         print( ec, bphc )
         zzz1 = zz.odein_solver(1., 0., ec, bphc)
         zzz2 = zz.odein_solver(1., 0.5, ec, bphc)
         print(zzz1, zzz2)    
         i = i + istep
         if zzz1['codein'] == 2 and zzz2['codein'] == 2:
           f.write( ('%.16e %.9e' % (ec, bphc) ) + '\n' )
           
    f.close()
    
  if not sswitch: # calculate the tidal deformability of neutron stars in the output of STGB_solver_5eqs; run time is 40 mins for one neutron star configuration
    tidprepdata = np.genfromtxt('stgb_rot_data1.txt' )  
    f = open('stgb_tid_v1_data1.txt', 'w+')
    f1 = open('stgb_tid_v1_error1.txt', 'w+')    
    if tidprepdata.size !=0:
      #print( len(tidprepdata) )
      istep = int(len(tidprepdata)/16)+1
      irun = 0
      while irun < len(tidprepdata):
         ec, bphcl, bphcu, M, moi = tidprepdata[irun, 0], tidprepdata[irun, 1], tidprepdata[irun, 2], tidprepdata[irun, 5], tidprepdata[irun, 8]
         print(irun, ec, bphcl, bphcu, M, moi)
         zzz = zz.tidal_solver(ec, bphcl)
         zzz1 = zz.tidal_solver(ec, bphcu)                  
         irun = irun + istep
         f.write( ('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (ec, bphcl, bphcu, M, moi, (zzz1['ytidal']+zzz['ytidal'])/2., (zzz1['dimllatidal']+zzz['dimllatidal'])/2., (zzz1['dimllatidal1']+zzz['dimllatidal1'])/2., (zzz1['dimllatidal2']+zzz['dimllatidal2'])/2., (zzz1['dimllatidal3']+zzz['dimllatidal3'])/2., (zzz1['dimllatidal4']+zzz['dimllatidal4'])/2., (zzz1['M_inf']+zzz['M_inf'])/2.-M, (zzz1['r_inf']+zzz['r_inf'])/2.  ) ) + '\n' ) 
         f1.write( ('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (ec, bphcl, bphcu, (zzz1['ytidal']-zzz['ytidal'])/2., (zzz1['dimllatidal']-zzz['dimllatidal'])/2., (zzz1['dimllatidal1']-zzz['dimllatidal1'])/2., (zzz1['dimllatidal2']-zzz['dimllatidal2'])/2., (zzz1['dimllatidal3']-zzz['dimllatidal3'])/2., (zzz1['dimllatidal4']-zzz['dimllatidal4'])/2.  ) ) + '\n' )         
    f.close()
    f1.close()
      
      
      


print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))

plt.show()
