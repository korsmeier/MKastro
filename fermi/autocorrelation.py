import numpy as np
import os
from scipy import integrate

import MKastro.fermi.observables        as fermiObs

CDIR = os.path.dirname(os.path.realpath(__file__))

data__ebins  = np.genfromtxt(CDIR+'/data/Cp/Cov_table__Ebins.dat.txt'  )
data__Cp     = np.genfromtxt(CDIR+'/data/Cp/Cov_table.dat.txt'         )
data__Cp_err = np.genfromtxt(CDIR+'/data/Cp/Cov_table__err.dat.txt'    )

n_bins = len(data__ebins[:,0])

def n():
    return n_bins
def get_Cp( ):
    return data__ebins, data__Cp, data__Cp_err

## Transformtions of the Cp matrix to a 1D array
# - Use the symmetry of the matrix nxn -> n(n+1)/2
# - Index k refers to the 1D array
# - Indices i,j refer to the 2D array

def _get_k( i, j ):
    if j>i:
        return get_k(j,i)
    return i*(i+1)/2 + j
get_k  = np.vectorize(_get_k)

def _get_ij( k ):
    i = int( -0.5 + np.sqrt( 0.25+2*k )+0.0000001  )
    j = k - i*(i+1)/2
    return i,j
get_ij = np.vectorize(_get_ij)



def powerlaw_integral(xmin, xmax, alpha):
    if np.fabs(1+alpha) < 1e-5:
        return np.ln( xmax/xmin )
    return 1./(1+alpha) * ( np.power( xmax, 1+alpha ) - np.power( xmin, 1+alpha ) )

######### Functions for simple Cp model
#
#   dN_ph/dE ~ E^-Gamma
#   dN/dS    = A (S/S_0)^-beta
#
#      S: Flux in the energy intervall 1GV - 100 GeV
#      S_0 = 1e-10 cm^2 s^-1
#
def _get_Cp_simple_model( Emin, Emax, E2min, E2max,  A, beta, Gamma, k=1.  ):
    S_0 = 1e-10
    F1000_cut_4FGL = fermiObs.F_thresh__Gamma_4FGL( Gamma )
    F1000_cut_3FHL = fermiObs.F_thresh__Gamma_3FHL( Gamma )
    
    F1000_cut_3FHL = fermiObs.F_E__from__F_Efrom( F_Efrom=F1000_cut_3FHL, Gamma=Gamma, Emin=1, Emax=100, Emin_from=10, Emax_from=1000 )
    
    cut1 = F1000_cut_4FGL
    cut2 = F1000_cut_4FGL
    if Emin  > 14.:
        cut1 = np.amin( [F1000_cut_4FGL, F1000_cut_3FHL] )
    if Emin  > 120.:
        cut1 = F1000_cut_3FHL
    if E2min > 14.:
        cut2 = np.amin( [F1000_cut_4FGL, F1000_cut_3FHL] )
    if E2min > 120.:
        cut2 = F1000_cut_3FHL
    F_cut = np.amin( [cut1, cut2] )

    integral_E1 = powerlaw_integral(  Emin,  Emax,          -Gamma       )
    integral_E2 = powerlaw_integral(  E2min, E2max,         -Gamma       )
    integral_E  = powerlaw_integral(  1.,    100. ,         -Gamma       )
    integral_F  = powerlaw_integral(  0.,    F_cut*k,       -beta+2    )
    return A/S_0**-beta * integral_E1*integral_E2/integral_E**2 * integral_F
get_Cp_simple_model = np.vectorize(_get_Cp_simple_model)


######### Functions for simple Cp model
#
#   dN_ph/dE ~ E^-Gamma
#   dN/dS    = A (S/S_0)^-beta * 1/sqrt(2pi)/sigma * exp( - (Gamma-mu)^2/2/sigam^2 )
#
#      S: Flux in the energy intervall 1GV - 100 GeV
#      S_0 = 1e-10 cm^2 s^-1
#
def _get_Cp_simple_model_GammaDistr( Emin, Emax, E2min, E2max,  A, beta, mu, sigma, k=1.  ):
    Gamma_min = np.amax([1.0, mu-5*sigma])
    Gamma_max = np.amin([3.0, mu+5*sigma])
    return integrate.quad(lambda Gamma: _get_Cp_simple_model(Emin, Emax, E2min, E2max,  A, Gamma=Gamma, beta=beta, k=k)*1./(sigma * np.sqrt(2 * np.pi)) * np.exp( - (Gamma - mu)**2 / (2 * sigma**2) ), Gamma_min, Gamma_max)[0]
get_Cp_simple_model_GammaDistr = np.vectorize(_get_Cp_simple_model_GammaDistr)
