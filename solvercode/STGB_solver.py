"""
solve sperical neutron stars in the scalar-tensor theory with the Gauss-Bonnet term
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
        """ y = [nu, mu1, bph, bphp, p] 
            equation arguments: (r, nu, mu1, bph, psi, p, e) 
        """
        dydr = [self.STGBeqs.dnudr(r, y[0], y[1], y[2], y[3], y[4], ediml( self.EOS.p2e( pdim(y[4]) ) )), 
                self.STGBeqs.dmu1dr(r, y[0], y[1], y[2], y[3], y[4], ediml( self.EOS.p2e( pdim(y[4]) ) )), 
                self.STGBeqs.dphdr(r, y[0], y[1], y[2], y[3], y[4], ediml( self.EOS.p2e( pdim(y[4]) ) )), 
                self.STGBeqs.ddphdr(r, y[0], y[1], y[2], y[3], y[4], ediml( self.EOS.p2e( pdim(y[4]) ) )), 
                self.STGBeqs.dpdr(r, y[0], y[1], y[2], y[3], y[4], ediml( self.EOS.p2e( pdim(y[4]) ) ))]
        return dydr

    def odesex(self, r, y):
        """ y = [nu, mu1, bph, bphp] 
            equation arguments: (r, nu, mu1, bph, psi, p, e) 
        """
        dydr = [self.STGBeqs.dnudr(r, y[0], y[1], y[2], y[3], 0., 0.), 
                self.STGBeqs.dmu1dr(r, y[0], y[1], y[2], y[3], 0., 0.), 
                self.STGBeqs.dphdr(r, y[0], y[1], y[2], y[3], 0., 0.), 
                self.STGBeqs.ddphdr(r, y[0], y[1], y[2], y[3], 0., 0.)]
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
        if y[4] <= pdiml(self.EOS.min_p):
            return -1
        else:
            return 0
            
    def solout1(self, r, y):
        """ Stop condition for ode """
        if y[4] <= 1.e-10:
            return -1
        else:
            return 0 
            
    def ode_solverprep(self, ec, bphc, Rmax = 3.0e6, rerror=1.e-7, verbose=False):
        dimlpc = pdiml(self.EOS.e2p(ec))
        y0 = [0., 0., bphc, 0., pdiml( self.EOS.e2p(ec) ) ]
        xmax = rdiml(Rmax)
        xini = 2.e-3 # does not work for 1.e-4, 1.e-5, ...

        R, M_s, p_s, codeex = 0, 0, 0, 0
        y_s = y0
        y_inft = [0, 0, 0, 0]
        
        [nu2, mu12, bph2, p2, check] = self.STGBeqs.coef2( dimlpc, ediml(ec), bphc)
        if check != 0:
          codein = -5 
        else:
          e2 = p2 / csqdiml(self.EOS.p2csq( self.EOS.e2p(ec) ))
          [nu4, mu14, bph4, p4] = self.STGBeqs.coef4( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2)
          """
          y0[0] = nu2*xini**2 
          y0[1] = mu12*xini**2 
          y0[2] = bphc + bph2*xini**2 
          y0[3] = 2.*bph2*xini 
          y0[4] = dimlpc + p2*xini**2 
          """
          y0[0] = nu2*xini**2 + nu4*xini**4
          y0[1] = mu12*xini**2 + mu14*xini**4
          y0[2] = bphc + bph2*xini**2 + bph4*xini**4
          y0[3] = 2.*bph2*xini + 4*bph4*xini**3
          y0[4] = dimlpc + p2*xini**2 + p4*xini**4

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

          x_s = solver.t
          y_s = solver.y
          R, M_s , p_s = runit*x_s, mdim(x_s*( 1. - np.exp(-2.*y_s[1]) )/2.), pdim(y_s[4])
          if x_s < 100*xini:
            codein = -4
 
          y_inft = [0, 0, 0, 0 ]      
          codeex = 0
          if codein==2 or codein==-2: 
            xinft = 100./np.sqrt(self.STGBeqs.msq)         # cannot be small; large is fine
            #print(xinft)
            """
            a = 0.0: 1000
            a = 0.1:  100
            a = 1: 50             
            """
            y0ex = [y_s[0], y_s[1], y_s[2], y_s[3]]
           
            solverex = sp_ode(self.odesex).set_integrator('dopri5', rtol=rerror)
            solverex.set_initial_value(y0ex, x_s)
            warnings.filterwarnings("ignore", category=UserWarning)
            solverex.integrate( xinft )
            codeex = sp_ode.get_return_code(solverex)
            warnings.resetwarnings()           
            y_inft = solverex.y            
                     
          else:
            if verbose:
              print('\n (-_-) Integration is not ended with solout, ' 'for ec = %.2e \n' % ec)
                                
        res = {'R': R, 'M_s': M_s, 'p_s': p_s, 'nu_s': y_s[0], 'mu_s': y_s[1], 'bph_s' : y_s[2], 'psi_s' : y_s[3], 'dimp_s': y_s[4], 'codein': codein, 'nu_inft': y_inft[0], 'mu_inft': y_inft[1], 'bph_inft' : y_inft[2], 'psi_inft' : y_inft[3], 'codeex': codeex }  
                     
        return res     
     
    def ode_solver(self, ec, bphc, Rmax = 3.0e6, rerror=1.e-7, verbose=False): # 1.e-8 is the best accuracy; 1.e-9 works but codein = -2 occurs frequently 
        dimlpc = pdiml(self.EOS.e2p(ec))
        y0 = [0., 0., bphc, 0., pdiml( self.EOS.e2p(ec) ) ]
        xmax = rdiml(Rmax)
        xini = 2.e-3 # does not work for 1.e-4, 1.e-5, ...
        
        [nu2, mu12, bph2, p2, check] = self.STGBeqs.coef2( dimlpc, ediml(ec), bphc)
        print(check)
        e2 = p2 / csqdiml(self.EOS.p2csq( self.EOS.e2p(ec) ))
        [nu4, mu14, bph4, p4] = self.STGBeqs.coef4( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2)
        #print( dimlpc, ediml(ec), bphc, nu2, mu12, bph2, p2, e2)
        #print( nu4, mu14, bph4, p4 )
        """
        y0[0] = nu2*xini**2 
        y0[1] = mu12*xini**2 
        y0[2] = bphc + bph2*xini**2 
        y0[3] = 2.*bph2*xini 
        y0[4] = dimlpc + p2*xini**2 
        """
        y0[0] = nu2*xini**2 + nu4*xini**4        
        y0[1] = mu12*xini**2 + mu14*xini**4
        y0[2] = bphc + bph2*xini**2 + bph4*xini**4
        y0[3] = 2.*bph2*xini + 4*bph4*xini**3
        y0[4] = dimlpc + p2*xini**2 + p4*xini**4

        print(y0, xini)
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
        #print(codein)

        x_s = solver.t
        y_s = solver.y
        print(x_s, y_s)
        R, M_s , p_s = runit*x_s, mdim(x_s*( 1. - np.exp(-2.*y_s[1]) )/2.), pdim(y_s[4])
        if x_s < 100*xini:   # stiff system
          codein = -4
                  
        if codein==2 or codein==-2 or codein==-4:
           
           xinft = 20/np.sqrt(self.STGBeqs.msq)            #  cannot be small; large is fine
           #print(xinfty)
           """
           a = 0.01: 1000
           a = 0.1:  100
           a = 1: 50
           """                     
           nex = 20000
           xexset = np.linspace(x_s+(xinft-x_s)/nex, xinft, nex)
           y0exset, y1exset, y2exset, y3exset = np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset)
           y0ex = [y_s[0], y_s[1], y_s[2], y_s[3]]
   
           
           solverex = sp_ode(self.odesex).set_integrator('dopri5', rtol=rerror)
           solverex.set_initial_value(y0ex, x_s)
           warnings.filterwarnings("ignore", category=UserWarning)
           for i in range(0, nex):
              solverex.integrate( xexset[i] )
              yf = solverex.y
              y0exset[i], y1exset[i], y2exset[i], y3exset[i] = yf[0], yf[1], yf[2], yf[3]
           
           
           i_infmu = np.argmin( abs(y1exset) )
           i_infbph = np.argmin( abs(y2exset) ) # suitable numerical infinity
           i_inf = i_infbph
           #print(i_infmu, i_infbph)
           dimlM_inf = xexset[i_inf]*( 1. - np.exp(-2.*y1exset[i_inf]) )/2.
           M_inf = mdim( dimlM_inf )
           #M_inf = M_inf/MSUN
           y0c_inf, y1c_inf = y0exset[i_inf], y1exset[i_inf]
           print(xexset[i_inf], dimlM_inf)

           if verbose: 
              nin = 1000
              xinset = np.linspace(xini, x_s, nin)
              y0inset, y1inset, y2inset, y3inset, y4inset, minset = np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset), np.zeros_like(xinset)             
              y0inset[0], y1inset[0], y2inset[0], y3inset[0], y4inset[0] = y0[0], y0[1], y0[2], y0[3], y0[4]
              
              solverin = sp_ode(self.odesin).set_integrator('dopri5', rtol=rerror) 
              solverin.set_solout(self.solout)  # stop condition
              solverin.set_initial_value(y0, xini)
              warnings.filterwarnings("ignore", category=UserWarning)
              for i in range(1, nin): 
                 solverin.integrate( xinset[i] )
                 yinf = solverin.y
                 y0inset[i], y1inset[i], y2inset[i], y3inset[i], y4inset[i] = yinf[0], yinf[1], yinf[2], yinf[3], yinf[4]
                 minset[i] = mdim(xinset[i]*( 1. - np.exp(-2.*y1inset[i]) )/2.)
                                               

              f = open("stgb_solver_sol_data.txt", 'w+')              
              for i in range(0, nin):            
                 y0inset[i] = y0inset[i] - y0c_inf-y1c_inf 
                 f.write( ('%.9f %.9f %.9f %.9f %.9f %.9f' % (xinset[[i]], y0inset[i], y1inset[i], y2inset[i], y3inset[i], y4inset[i]  ) ) + '\n' )                   
                 
              mexset = np.zeros_like(xexset)
              for i in range(0, nex): 
                 #print(i)
                 #print(y0exset[i], - y0exset[i_inf] - y1exset[i_inf], -y0c_inf-y1c_inf )           
                 y0exset[i] = y0exset[i] - y0c_inf-y1c_inf 
                 mexset[i] = mdim(xexset[i]*( 1. - np.exp(-2.*y1exset[i]) )/2.)
                 f.write( ('%.9f %.9f %.9f %.9f %.9f %.9f' % (xexset[[i]], y0exset[i], y1exset[i], y2exset[i], y3exset[i], 0.  ) ) + '\n' )                 
                 
                   
              f.close()
              
                                         

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
                     
              ax1.plot(xinset, y0inset, 'o')
              ax2.plot(xinset, y1inset)
              ax3.plot(xinset, y2inset)
              ax4.plot(xinset, y3inset)
              ax5.plot(xinset, y4inset)
              """                
              ax1.plot(xexset[0:i_inf+1], y2exset[0:i_inf+1])
              ax3.plot(xexset[0:i_inf+1], y0exset[0:i_inf+1])
              ax2.plot(xexset[0:i_inf+1], mexset[0:i_inf+1])
              ax4.plot(xexset[0:i_inf+1], y1exset[0:i_inf+1])
              """
             
              ax1.plot(xexset, y0exset, 'o')
              ax2.plot(xexset, y1exset)
              ax3.plot(xexset, y2exset)
              ax4.plot(xexset, y3exset)
              
              
              
              
              """
              # solve exterior eqs using nin-1 interior values as the boundary condition
              y0exset1, y1exset1, y2exset1, y3exset1 = np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset), np.zeros_like(xexset)
              y0exset1[0], y1exset1[0], y2exset1[0], y3exset1[0] = y0inset[nin-1], y1inset[nin-1], y2inset[nin-1], y3inset[nin-1]
              y0ex1 = [y0inset[nin-1], y1inset[nin-1], y2inset[nin-1], y3inset[nin-1]] 
           
              solverex1 = sp_ode(self.odesex).set_integrator('dopri5')
              solverex1.set_initial_value(y0ex1, x_s)
              warnings.filterwarnings("ignore", category=UserWarning)
              for i in range(1, nex):
                 solverex1.integrate( xexset[i] )
                 yf1 = solverex1.y
                 y0exset1[i], y1exset1[i], y2exset1[i], y3exset1[i] = yf1[0], yf1[1], yf1[2], yf1[3]
                  
              i_inf1 = np.argmin( abs(y2exset1) ) 
              M_inf1 = mdim(xexset[i_inf1]*( 1. - np.exp(-2.*y1exset1[i_inf1]) )/2.)
              print(xinset[nin-1], y0inset[nin-1], y1inset[nin-1], y2inset[nin-1], y3inset[nin-1], y4inset[nin-1])
              print(x_s, y_s)
              
              print(xexset[i_inf1], y0exset1[i_inf1], y1exset1[i_inf1], y2exset1[i_inf1], y3exset1[i_inf1], M_inf1 )
              """
                     
        else:
           if verbose:
              print('\n (-_-) Integration is not ended with solout, ' 'for ec = %.6e \n' % ec)
              
                
        res = {'R': R, 'M_s': M_s, 'p_s': p_s, 'nu_s': y_s[0]-y0c_inf-y1c_inf, 'mu_s': y_s[1], 'bph_s' : y_s[2], 'psi_s' : y_s[3], 'dimp_s': y_s[4], 'codein': codein, 'r_inf': xexset[i_inf], 'nu_inf': -y1exset[i_inf], 'mu_inf': y1exset[i_inf], 'bph_inf' : y2exset[i_inf], 'psi_inf' : y3exset[i_inf], 'M_inf': M_inf}  
                     
        return res

    def ode_solsearchprep(self, ec, verbose=False):        # preliminary search
        #print(self.STGBeqs.xi) 
        if self.STGBeqs.xi < 0:
           bphcmin = 0.001
           bphcstep = 0.001  # for xi < 0
        else:
           bphcmin = 0.00001
           bphcstep = 0.00001 # for xi > 0 
             
        bphctest = bphcmin               
        conroots = np.roots( self.bph2coef(ec, bphctest) )
        sol = self.ode_solverprep(ec, bphctest)
        
        codein = sol['codein']  
        bph_inft = sol['bph_inft'] 
        bph_inft1 = bph_inft         
        print( bphctest, conroots[3], bph_inft )
        
        signcount = 0
        bphclower, bphcupper = 0., 0.
        while (codein==2 or codein==-2) and np.imag(conroots[3]) == 0.0 and bphctest<1. and signcount<1:
           bphctest = bphctest + bphcstep
           conroots = np.roots( self.bph2coef(ec, bphctest) )
           sol = self.ode_solverprep(ec, bphctest)
           codein = sol['codein']  
           bph_inft = sol['bph_inft'] 
           print( bphctest, conroots[3], codein, bph_inft ) 

           if np.imag(conroots[3]) != 0.0:
              break
                         
           if bph_inft * bph_inft1 < 0.0:
              signcount = signcount + 1
              if signcount == 1:  
                bphclower, bphcupper = bphctest-bphcstep, bphctest
           bph_inft1 = bph_inft
           
          
        return {'ec': ec, 'signcount': signcount, 'codein': codein, 'bphclower': bphclower, 'bphcupper': bphcupper}
                
        
    def ode_solsearch(self, ec, bphcl, bphch, verbose=False):       # search for high accuracy; bisection method fails for accuracy at 1e-7 
        if self.STGBeqs.xi < 0:    
           bphcmin = bphcl
           bphcstep = 1.e-5 # for xi < 0 
        else:
           bphcmin = bphcl
           bphcstep = 1.e-7 # for xi > 0 

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
    
        
    def search_onsetec(self, ec1, ec2, verbose=False):        
        if self.STGBeqs.xi < 0:    
          bphcmin = 1.e-5 # for xi < 0 
        else:
          bphcmin = 1.e-7 # for xi > 0  
                    
        conroots = np.roots( self.bph2coef(ec1, bphcmin) )
        sol = self.ode_solverprep(ec1, bphcmin)
        codein, bph_inft1 = sol['codein'], sol['bph_inft'] 
        ec, bph_inft = ec1, bph_inft1         
        print(ec1, codein, conroots[3], bph_inft1 )
        conroots = np.roots( self.bph2coef(ec2, bphcmin) )
        sol = self.ode_solverprep(ec2, bphcmin)
        codein, bph_inft1 = sol['codein'], sol['bph_inft']         
        print(ec2, codein, conroots[3], bph_inft1 ) 
        
        if bph_inft * bph_inft1 < 0.0:
          print('input okay')
          while (codein==2 or codein==-2) and np.imag(conroots[3]) == 0.0 and abs(ec - (ec1 + ec2)/2.)/ec1 > 1.e-6:
             ec = (ec1 + ec2)/2.
             conroots = np.roots( self.bph2coef(ec, bphcmin) )
             sol = self.ode_solverprep(ec, bphcmin)
             codein, bph_inft1 = sol['codein'], sol['bph_inft']
             print( ec, codein, conroots[3], bph_inft1 )     
                           
             if bph_inft * bph_inft1 > 0.0:        
               ec1 = ec
             else:
               ec2 = ec 
          print('bisection search finished')   
               
          ec, ecstep = ec2, (ec2-ec1)/10.
          bphcmax = 2*bphcmin        
          conroots = np.roots( self.bph2coef(ec2, bphcmin) )
          sol = self.ode_solverprep(ec2, bphcmin)
          codein, bph_inft1 = sol['codein'], sol['bph_inft']         
          print(ec2, codein, conroots[3], bph_inft1 ) 
          conroots = np.roots( self.bph2coef(ec2, bphcmax) )
          sol = self.ode_solverprep(ec2, bphcmax)
          codein, bph_inft2 = sol['codein'], sol['bph_inft']         
          print(ec2, codein, conroots[3], bph_inft2 ) 

          while bph_inft1*bph_inft2 > 0 and ec>ec1:
             ec = ec - ecstep
             conroots = np.roots( self.bph2coef(ec, bphcmin) )
             sol = self.ode_solverprep(ec, bphcmin)
             codein, bph_inft1 = sol['codein'], sol['bph_inft']         
             print(ec, codein, conroots[3], bph_inft1 ) 
             conroots = np.roots( self.bph2coef(ec, bphcmax) )
             sol = self.ode_solverprep(ec, bphcmax)
             codein, bph_inft2 = sol['codein'], sol['bph_inft']         
             print(ec, codein, conroots[3], bph_inft2 )  
        
        else:
          print('invalid input')
                                                
        return {'ec1':ec1, 'ec2':ec2, 'ec':ec}    
                 
        
        
        
        
        
        
t0 = timeit.time.time()

eos_nameset = ['WFF1','SLy4','AP4','MPA1','PAL1']
xiset = [3., 5., 10., -0.2, -0.5, -1.]
msqset = [0.01**2, 0.1**2, 1.**2]
ieos = 2
msq = 1.**2
xi = -1.
eos_name = eos_nameset[ieos]


switch = 1           
sswitch = 1

zz = STGB_solver(EOS_name=eos_name, theory=STGBeqs(xi, msq) )

if switch:
  if sswitch: # single run for testing use
    ec = 2.319000E+15
    bphc = 0.01388
    #ec = 3.8377045727049740e+15
    #bphc = 4.620000000e-03
    #zzz = zz.ode_solverprep(ec, bphc)
    zzz = zz.ode_solver(ec, bphc, verbose=True)
    #print( zzz )
    #print(('%.16e %.9e %.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (ec, bphc, zzz['R'], zzz['M_inf'], zzz['nu_s'], zzz['M_s'], zzz['bph_s'], zzz['psi_s'], zzz['dimp_s'], zzz['r_inf'] ) ))
    #print( zz.ode_solsearchprep(ec) )
    #zzz = zz.ode_solsearch(ec, 0.00001, 0.00002)
    #sol = zz.ode_solver(ec, zzz['bphclower'])
    #print(sol['M_inf']/MSUN)

  if not sswitch: # single run; find the onset ec
    ec1, ec2 = 2.7113449507015170e+15, 2.7123449507015170e+15
    zzz = zz.search_onsetec(ec1, ec2)
    print(zzz)
  
  
if not switch:
  if sswitch:  # find the brief bounds for suitable values of the central scalar; multiple runs for hpc use
    ecstepset = [ 0.1e+15, 0.01e+15, 0.001e+15, 0.0001e+15, 0.00001e+15] 
      
    f = open("stgb_solverprep_data13.txt", 'w+')     
    stadata = np.genfromtxt('stgb_solversta_data13.txt', usecols = (0) )
    if stadata.size !=0:
      if stadata.size == 1:
        for (j, ecstep) in enumerate(ecstepset):
           print(j)
           ectest = stadata - ecstep
           bphcsol = zz.ode_solsearchprep(ectest)            
           if bphcsol['bphcupper'] != 0:
             break 
           else:
             ectest = stadata + ecstep
             bphcsol = zz.ode_solsearchprep(ectest)
             if bphcsol['bphcupper'] != 0:
               break
                           
        if bphcsol['bphcupper'] != 0:
          f.write( ('%.16e %.9f %.9f' % (ectest, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' )             
          ecstep1 = ecstep /10.
          ec = ectest           
          while bphcsol['bphcupper'] != 0:
             ec = ec - ecstep1
             bphcsol = zz.ode_solsearchprep(ec)
             f.write( ('%.16e %.9f %.9f' % (ec, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' )    
          ec = ectest + ecstep1
          bphcsol = zz.ode_solsearchprep(ec)                     
          f.write( ('%.16e %.9f %.9f' % (ec, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' ) 
          while bphcsol['bphcupper'] != 0:         
             ec = ec + ecstep1
             bphcsol = zz.ode_solsearchprep(ec)           
             f.write( ('%.16e %.9f %.9f' % (ec, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' ) 
                        
      else:
        for (i, ec0) in enumerate(stadata):
           print(i)
           for (j, ecstep) in enumerate(ecstepset):
              print(j)
              ectest = ec0 - ecstep
              bphcsol = zz.ode_solsearchprep(ectest)            
              if bphcsol['bphcupper'] != 0:
                break 
              else:
                ectest = ec0 + ecstep
                bphcsol = zz.ode_solsearchprep(ectest)
                if bphcsol['bphcupper'] != 0:
                  break

           if bphcsol['bphcupper'] != 0:
             f.write( ('%.16e %.9f %.9f' % (ectest, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' )                    
             ecstep1 = ecstep /5.
             ec = ectest           
             while bphcsol['bphcupper'] != 0:
                ec = ec - ecstep1
                bphcsol = zz.ode_solsearchprep(ec)
                f.write( ('%.16e %.9f %.9f' % (ec, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' )    
             ec = ectest + ecstep1
             bphcsol = zz.ode_solsearchprep(ec)   
             f.write( ('%.16e %.9f %.9f' % (ec, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' )                   
             while bphcsol['bphcupper'] != 0:         
                ec = ec + ecstep1
                bphcsol = zz.ode_solsearchprep(ec)
                f.write( ('%.16e %.9f %.9f' % (ec, bphcsol['bphclower'], bphcsol['bphcupper'] ) ) + '\n' ) 
                                                              
    f.close()
    
    
  if not sswitch: # use the brief bounds on the central scalar to find solutions of spontaneous scalarization
    f = open("stgb_solver_data2.txt", 'w+') 
    f1 = open("stgb_solver_error2.txt", 'w+')     
    seadata = np.genfromtxt('stgb_solverprep_data2.txt')
    if seadata.size !=0:
      c0, c1, c2 = seadata[:, 0], seadata[:, 1], seadata[:, 2]
      nec = len(seadata)       
      for i in range(0, nec):
         if c2[i] != 0: 
           print(i)         
           zzz = zz.ode_solsearch(c0[i], c1[i], c2[i])
           odesol1 = zz.ode_solver(c0[i], zzz['bphclower'])
           odesol2 = zz.ode_solver(c0[i], zzz['bphcupper'])
           f.write( ('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (c0[i], zzz['bphclower'], zzz['bphcupper'], (odesol1['R']+odesol2['R'])/2., (odesol1['M_inf']+odesol2['M_inf'])/2., (odesol1['nu_s']+odesol2['nu_s'])/2., (odesol1['M_s']+odesol2['M_s'])/2., (odesol1['bph_s']+odesol2['bph_s'])/2., (odesol1['psi_s']+odesol2['psi_s'])/2., (odesol1['dimp_s']+odesol2['dimp_s'])/2., (odesol1['r_inf']+odesol2['r_inf'])/2. ) ) + '\n' ) 
           f1.write( ('%.16e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e' % (c0[i], zzz['bphclower'], zzz['bphcupper'], (odesol2['R']-odesol1['R'])/2., (odesol2['M_inf']-odesol1['M_inf'])/2., (odesol2['nu_s']-odesol1['nu_s'])/2., (odesol2['M_s']-odesol1['M_s'])/2., (odesol2['bph_s']-odesol1['bph_s'])/2., (odesol2['psi_s']-odesol1['psi_s'])/2., (odesol2['dimp_s']-odesol1['dimp_s'])/2., (odesol2['r_inf']-odesol1['r_inf'])/2. ) ) + '\n' )                                  
    f.close()
    f1.close()    




print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))


plt.show()
