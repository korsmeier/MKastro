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
    return int( i*(i+1)/2 + j + 0.0000001 )
get_k  = np.vectorize(_get_k)

def _get_ij( k ):
    i = int( -0.5 + np.sqrt( 0.25+2*k )+0.0000001  )
    j = k - (i*(i+1))/2
    return int(i+0.0000001),int(j+0.0000001)
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

    #print('%e    %e            %e' % (F1000_cut_4FGL,F1000_cut_3FHL,F_cut))

    integral_E1 = powerlaw_integral(  Emin,  Emax,          -Gamma       )
    integral_E2 = powerlaw_integral(  E2min, E2max,         -Gamma       )
    integral_E  = powerlaw_integral(  1.,    100. ,         -Gamma       )
    integral_F  = powerlaw_integral(  0.,    F_cut*k,       -beta+2    )
    return A/S_0**-beta * integral_E1*integral_E2/integral_E**2 * integral_F
get_Cp_simple_model = np.vectorize(_get_Cp_simple_model)


######### Functions for simple Cp model
#
#   dN_ph/dE      ~ E^-Gamma
#   dN/(dSdGamma) = A (S/S_0)^-beta * 1/sqrt(2pi)/sigma * exp( - (Gamma-mu)^2/2/sigam^2 )
#
#      S: Flux in the energy intervall 1GV - 100 GeV
#      S_0 = 1e-10 cm^2 s^-1
#
def _get_Cp_simple_model_GammaDistr( Emin, Emax, E2min, E2max,  A, beta, mu, sigma, k=1.  ):
    Gamma_min = np.amax([1.0, mu-5*sigma])
    Gamma_max = np.amin([3.0, mu+5*sigma])
    return integrate.quad(lambda Gamma: _get_Cp_simple_model(Emin, Emax, E2min, E2max,  A, Gamma=Gamma, beta=beta, k=k)*1./(sigma * np.sqrt(2 * np.pi)) * np.exp( - (Gamma - mu)**2 / (2 * sigma**2) ), Gamma_min, Gamma_max)[0]
get_Cp_simple_model_GammaDistr = np.vectorize(_get_Cp_simple_model_GammaDistr)



#
#   integrand ~ x^alpha * exp(-x/x_b)
#
def powerlaw_cutoff_integral(xmin, xmax, alpha, x_b):
    res = integrate.quad(lambda logx: np.power(np.exp(logx),1+alpha) * np.exp( -np.exp(logx)/x_b ), np.log(xmin), np.log(xmax))
    return res[0]

######### Functions for simple Cp model with exponential cutoff
#
#   dN_ph/dE ~ E^-Gamma * exp(-E/E_b)
#   dN/dS    = A (S/S_0)^-beta
#
#      S: Flux in the energy intervall 1GV - 100 GeV
#      S_0 = 1e-10 cm^2 s^-1
#
def _get_Cp_simple_cutoff_model( Emin, Emax, E2min, E2max,  A, beta, Gamma, E_b, k=1.  ):
    S_0 = 1e-10
    F1000_cut_4FGL  = fermiObs.F_thresh__Gamma_4FGL( Gamma )
    F10000_cut_3FHL = fermiObs.F_thresh__Gamma_3FHL( Gamma )
    
    integral_1_100    = powerlaw_cutoff_integral(   1.,  100., -Gamma, E_b)
    integral_10_1000  = powerlaw_cutoff_integral(  10., 1000., -Gamma, E_b)
    
    F1000_cut_3FHL =  integral_1_100 / integral_10_1000 * F10000_cut_3FHL
    
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

    #print('%e    %e            %e' % (F1000_cut_4FGL,F1000_cut_3FHL,F_cut))

    integral_E1 = powerlaw_cutoff_integral(  Emin,  Emax,   -Gamma, E_b  )
    integral_E2 = powerlaw_cutoff_integral(  E2min, E2max,  -Gamma, E_b  )
    integral_E  = integral_1_100

    integral_F  = powerlaw_integral(  0.,    F_cut*k,       -beta+2    )

    return A/S_0**-beta * integral_E1*integral_E2/integral_E**2 * integral_F
get_Cp_simple_cutoff_model = np.vectorize(_get_Cp_simple_cutoff_model)

######### Functions for simple Cp model with exponential cutoff
#
#   dN_ph/dE      ~ E^-Gamma * exp(-E/E_b)
#   dN/(dSdGamma) = A (S/S_0)^-beta * 1/sqrt(2pi)/sigma * exp( - (Gamma-mu)^2/2/sigam^2 )
#
#      S: Flux in the energy intervall 1GV - 100 GeV
#      S_0 = 1e-10 cm^2 s^-1
#
def _get_Cp_simple_cutoff_model_GammaDistr( Emin, Emax, E2min, E2max,  A, beta, mu, sigma, E_b, k=1.  ):
    Gamma_min = np.amax([1.0, mu-5*sigma])
    Gamma_max = np.amin([3.0, mu+5*sigma])
    return integrate.quad(lambda Gamma: _get_Cp_simple_cutoff_model(Emin, Emax, E2min, E2max,  A=A, Gamma=Gamma, beta=beta, E_b=E_b, k=k)*1./(sigma * np.sqrt(2 * np.pi)) * np.exp( - (Gamma - mu)**2 / (2 * sigma**2) ), Gamma_min, Gamma_max)[0]
get_Cp_simple_cutoff_model_GammaDistr = np.vectorize(_get_Cp_simple_cutoff_model_GammaDistr)


######### Functions for simple Cp model - in the limits exploit simpler calculations - errors at the percent level
#
#   dN_ph/dE      ~ E^-Gamma * exp(-E/E_b)
#   dN/(dSdGamma) = A (S/S_0)^-beta * 1/sqrt(2pi)/sigma * exp( - (Gamma-mu)^2/2/sigam^2 )
#
#      S: Flux in the energy intervall 1GV - 100 GeV
#      S_0 = 1e-10 cm^2 s^-1
#
#
#   IF E_b > 100 * max(Emax, E2max) -> ignore exponential cutoff
#   IF sigma < 0.01                 -> ignore Gamma Distr.
#
def _get_Cp_simple_model_optimized                         ( Emin, Emax, E2min, E2max, A, beta, mu,       sigma=0.0,      E_b=1e5, k=1.  ):
    #
    if sigma <= 0.01:
        if E_b/np.amax([Emax, E2max, 100.]) >= 100:
            return _get_Cp_simple_model                    ( Emin, Emax, E2min, E2max, A, beta, Gamma=mu,                          k=k )
        else:
            return _get_Cp_simple_cutoff_model             ( Emin, Emax, E2min, E2max, A, beta, Gamma=mu,                 E_b=E_b, k=k )
    else:
        if E_b/np.amax([Emax, E2max, 100.]) >= 100:
            return _get_Cp_simple_model_GammaDistr         ( Emin, Emax, E2min, E2max, A, beta, mu=mu   , sigma=sigma,             k=k )
        else:
            return _get_Cp_simple_cutoff_model_GammaDistr  ( Emin, Emax, E2min, E2max, A, beta, mu=mu   , sigma=sigma,    E_b=E_b, k=k )
get_Cp_simple_model_optimized = np.vectorize(_get_Cp_simple_model_optimized)



######### Functions for simple dN/dS model
#
#   dN_ph/dE ~ E^-Gamma
#   dN/dS    = A (S/S_0)^-beta
#
#      S: Flux in the [Emin, Emax];  [E]=GeV, [S]=cm^-2 s^-1
#      S_0 = 1e-10 cm^-2 s^-1
#
#      [A]     = cm^-2 s^-1 sr
#      [dN/dS] = cm^-2 s^-1 sr^-1
#
def _get_dNdS_simple_model( S, Emin, Emax, A, beta, Gamma ):
    S_0 = 1e-10
    integral_E     = powerlaw_integral(  Emin,  Emax,          -Gamma       )
    integral_1000  = powerlaw_integral(  1.,    100. ,         -Gamma       )
    f = integral_E/integral_1000
    return A * ( S / S_0 )**-beta * f**(beta-1)
get_dNdS_simple_model = np.vectorize(_get_dNdS_simple_model)


######### Functions for simple dN/dS model with exponential cutoff
#
#   dN_ph/dE ~ E^-Gamma * exp(-E/E_b)
#   dN/dS    = A (S/S_0)^-beta
#
#      S: Flux in the [Emin, Emax];  [E]=GeV, [S]=cm^-2 s^-1
#      S_0 = 1e-10 cm^-2 s^-1
#
#      [A]     = cm^-2 s^-1 sr
#      [dN/dS] = cm^-2 s^-1 sr^-1
#
def _get_dNdS_simple_cutoff_model( S, Emin, Emax, A, beta, Gamma, E_b ):
    S_0 = 1e-10
    integral_1000  = powerlaw_cutoff_integral(    1.,  100.,   -Gamma, E_b  )
    integral_E     = powerlaw_cutoff_integral(  Emin,  Emax,   -Gamma, E_b  )
    f = integral_E/integral_1000
    return A * ( S / S_0 )**-beta * f**(beta-1)
get_dNdS_simple_cutoff_model = np.vectorize(_get_dNdS_simple_cutoff_model)







######## Tests Cp predictions
#
#import matplotlib                  as mpl
#mpl.use('Agg')
#import matplotlib.pyplot           as plt
#import matplotlib.colors           as colors
#import MKastro.basic.PlotFunctions as pf
#
#A     = 1
#beta  = 1.8
#Gamma = 2.3
#sigma = 0.01
#E_b   = 80
#
### Plot Cp diagonla
#plot, fig = pf.new_plot(r'$E\;\;\mathrm{[GeV]}$', r'$\frac{E^4}{(\Delta E)^2} C_p\;\;\mathrm{[GeV^2 cm^{-4} s^{-2} sr^{-1}]}$', 'log', 'log', label_size=18)
#
#Emin     =  data__ebins[:,0]
#Emax     =  data__ebins[:,1]
#E        =  np.sqrt(Emin*Emax)
#DeltaE   =  Emax - Emin
#fac      =  E**4/DeltaE**2
#xerr     =  [E-Emin, Emax-E ]
#
#Cp_simple               = get_Cp_simple_model                   ( Emin, Emax, Emin, Emax,  A, beta, Gamma,                          k=1.  )
#Cp_simple_cutoff        = get_Cp_simple_cutoff_model            ( Emin, Emax, Emin, Emax,  A, beta, Gamma,                 E_b,     k=1.  )
#Cp_simple_Gamma         = get_Cp_simple_model_GammaDistr        ( Emin, Emax, Emin, Emax,  A, beta, mu=Gamma, sigma=sigma,          k=1.  )
#Cp_simple_cutoff_Gamma  = get_Cp_simple_cutoff_model_GammaDistr ( Emin, Emax, Emin, Emax,  A, beta, mu=Gamma, sigma=sigma, E_b=E_b, k=1.  )
#
#print( np.amax( np.fabs( 1 - Cp_simple/Cp_simple_Gamma ) )               )
#print( np.amax( np.fabs( 1 - Cp_simple_cutoff/Cp_simple_cutoff_Gamma ) ) )
#
#plot.errorbar( E, fac * Cp_simple             ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple'                     )
#plot.errorbar( E, fac * Cp_simple_cutoff      ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff'              )
#plot.errorbar( E, fac * Cp_simple_Gamma       ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple Gamma distr.'        )
#plot.errorbar( E, fac * Cp_simple_cutoff_Gamma,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff Gamma distr.' )
#
#plot.legend(frameon=False, loc='lower left', bbox_to_anchor=(0.05, 0.05), ncol=1, fontsize=15)
#plt.savefig('test_Cp_diagonal.pdf' )
#
#
### Plot Cp cross energies in one plot per energy bin
#for j in range(n_bins):
#    plot, fig = pf.new_plot(r'$E_j\;\;\mathrm{[GeV]}$', r'$\frac{E_i^2 E_j^2}{\Delta E_i\Delta E_j} C_p(E_i,E_j)\;\;\mathrm{[GeV^2 cm^{-4} s^{-2} sr^{-1}]}$', 'log', 'log', label_size=18)
#    Emin1     = data__ebins[j,0]
#    Emax1     = data__ebins[j,1]
#    E_1       = np.sqrt(Emin1*Emax1)
#    Emin2     = data__ebins[:,0]
#    Emax2     = data__ebins[:,1]
#    E_2       = np.sqrt(Emin2*Emax2)
#    DeltaE1   = Emax1 - Emin1
#    DeltaE2   = Emax2 - Emin2
#    fac       = E_1**2*E_2**2/DeltaE1/DeltaE2
#    xerr      = [E_2-Emin2, Emax2-E_2 ]
#
#    Cp_simple               = get_Cp_simple_model                   ( Emin1, Emax1, Emin2, Emax2,  A, beta, Gamma,                          k=1.  )
#    Cp_simple_cutoff        = get_Cp_simple_cutoff_model            ( Emin1, Emax1, Emin2, Emax2,  A, beta, Gamma,                 E_b,     k=1.  )
#    Cp_simple_Gamma         = get_Cp_simple_model_GammaDistr        ( Emin1, Emax1, Emin2, Emax2,  A, beta, mu=Gamma, sigma=sigma,          k=1.  )
#    Cp_simple_cutoff_Gamma  = get_Cp_simple_cutoff_model_GammaDistr ( Emin1, Emax1, Emin2, Emax2,  A, beta, mu=Gamma, sigma=sigma, E_b=E_b, k=1.  )
#
#
#    print( np.amax( np.fabs( 1 - Cp_simple/Cp_simple_Gamma ) )               )
#    print( np.amax( np.fabs( 1 - Cp_simple_cutoff/Cp_simple_cutoff_Gamma ) ) )
#
#    plot.errorbar( E_2 , fac * Cp_simple             ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple'                     )
#    plot.errorbar( E_2 , fac * Cp_simple_cutoff      ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff'              )
#    plot.errorbar( E_2 , fac * Cp_simple_Gamma       ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple Gamma distr.'        )
#    plot.errorbar( E_2 , fac * Cp_simple_cutoff_Gamma,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff Gamma distr.' )
#
#    #plot.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, fontsize=15)
#    plot.legend(frameon=False, loc='lower left', bbox_to_anchor=(0.05, 0.05), ncol=1, fontsize=15)
#    plt.savefig('test_Cp_crossE_data_%i.pdf' % j )
#
#
