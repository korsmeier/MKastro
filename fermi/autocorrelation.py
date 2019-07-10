import numpy as np
import os

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



######### Functions for simple Cp model
#
#
#   A E^-beta S^-gamma
#
#

def powerlaw_integral(xmin, xmax, alpha):
    if np.fabs(1+alpha) < 1e-5:
        return np.ln( xmax/xmin )
    return 1./(1+alpha) * ( np.power( xmax, 1+alpha ) - np.power( xmin, 1+alpha ) )

def _get_Cp_simple_model( Emin, Emax, E2min, E2max,  A, beta, gamma, k=1.  ):
    F1000_cut_4FGL = fermiObs.F_thresh__Gamma_4FGL( beta )
    F1000_cut_3FHL = fermiObs.F_thresh__Gamma_3FHL( beta )
    
    F1000_cut_3FHL = fermiObs.F_E__from__F_Efrom( F_Efrom=F1000_cut_3FHL, Gamma=gamma, Emin=1, Emax=100, Emin_from=10, Emax_from=1000 )
    
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

    integral_E1 = powerlaw_integral(  Emin,  Emax,          -beta       )
    integral_E2 = powerlaw_integral(  E2min, E2max,         -beta       )
    integral_E  = powerlaw_integral(  1.,    100. ,         -beta       )
    integral_F  = powerlaw_integral(  0.,    F_cut*k,       -gamma+2    )
    return A*integral_E1*integral_E2/integral_E**2 * integral_F

get_Cp_simple_model = np.vectorize(_get_Cp_simple_model)
