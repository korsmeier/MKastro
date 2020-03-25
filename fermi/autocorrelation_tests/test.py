import numpy as np
import sys
sys.path.append('..')
import autocorrelation as AC


import matplotlib                  as mpl
mpl.use('Agg')
import matplotlib.pyplot           as plt
import matplotlib.colors           as colors
import MKastro.basic.PlotFunctions as pf


####### Tests Cp predictions

data__ebins, data__Cp, data__Cp_err = AC.get_Cp( )
n_bins                              = AC.n     ( )

A     = 1
beta  = 1.8
Gamma = 2.3
sigma = 0.2
E_b   = 80

## Plot Cp diagonla
plot, fig = pf.new_plot(r'$E\;\;\mathrm{[GeV]}$', r'$\frac{E^4}{(\Delta E)^2} C_p\;\;\mathrm{[GeV^2 cm^{-4} s^{-2} sr^{-1}]}$', 'log', 'log', label_size=18)

Emin     =  data__ebins[:,0]
Emax     =  data__ebins[:,1]
E        =  np.sqrt(Emin*Emax)
DeltaE   =  Emax - Emin
fac      =  E**4/DeltaE**2
xerr     =  [E-Emin, Emax-E ]

Cp_simple               = AC.get_Cp_simple_model                   ( Emin, Emax, Emin, Emax,  A, beta, Gamma,                          k=1.  )
Cp_simple_cutoff        = AC.get_Cp_simple_cutoff_model            ( Emin, Emax, Emin, Emax,  A, beta, Gamma,                 E_b,     k=1.  )
Cp_simple_Gamma         = AC.get_Cp_simple_model_GammaDistr        ( Emin, Emax, Emin, Emax,  A, beta, mu=Gamma, sigma=sigma,          k=1.  )
Cp_simple_cutoff_Gamma  = AC.get_Cp_simple_cutoff_model_GammaDistr ( Emin, Emax, Emin, Emax,  A, beta, mu=Gamma, sigma=sigma, E_b=E_b, k=1.  )


v2_Cp_simple               = AC.get_Cp_simple_model_optimized         ( Emin, Emax, Emin, Emax,  A, beta, Gamma,    sigma=0.0  , E_b=1e6, k=1.  )
v2_Cp_simple_cutoff        = AC.get_Cp_simple_model_optimized         ( Emin, Emax, Emin, Emax,  A, beta, Gamma,    sigma=0.0  , E_b=E_b, k=1.  )
v2_Cp_simple_Gamma         = AC.get_Cp_simple_model_optimized         ( Emin, Emax, Emin, Emax,  A, beta, mu=Gamma, sigma=sigma, E_b=1e6, k=1.  )
v2_Cp_simple_cutoff_Gamma  = AC.get_Cp_simple_model_optimized         ( Emin, Emax, Emin, Emax,  A, beta, mu=Gamma, sigma=sigma, E_b=E_b, k=1.  )



#print( np.amax( np.fabs( 1 - Cp_simple/Cp_simple_Gamma ) )               )
#print( np.amax( np.fabs( 1 - Cp_simple_cutoff/Cp_simple_cutoff_Gamma ) ) )

plot.errorbar( E, fac * Cp_simple             ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple'                     )
plot.errorbar( E, fac * Cp_simple_cutoff      ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff'              )
plot.errorbar( E, fac * Cp_simple_Gamma       ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple Gamma distr.'        )
plot.errorbar( E, fac * Cp_simple_cutoff_Gamma,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff Gamma distr.' )

plot.legend(frameon=False, loc='lower left', bbox_to_anchor=(0.05, 0.05), ncol=1, fontsize=15)
plt.savefig('test_Cp_diagonal.pdf' )


plot.errorbar( E, fac * v2_Cp_simple             ,   lw=2, yerr=0, xerr=0,    fmt='o', markersize=10, label='Cp simple'                     )
plot.errorbar( E, fac * v2_Cp_simple_cutoff      ,   lw=2, yerr=0, xerr=0,    fmt='o', markersize=10, label='Cp simple cutoff'              )
plot.errorbar( E, fac * v2_Cp_simple_Gamma       ,   lw=2, yerr=0, xerr=0,    fmt='o', markersize=10, label='Cp simple Gamma distr.'        )
plot.errorbar( E, fac * v2_Cp_simple_cutoff_Gamma,   lw=2, yerr=0, xerr=0,    fmt='o', markersize=10, label='Cp simple cutoff Gamma distr.' )


plot.legend(frameon=False, loc='lower left', bbox_to_anchor=(0.05, 0.05), ncol=1, fontsize=15)
plt.savefig('test_Cp_diagonal_optimized.pdf' )


## Plot Cp cross energies in one plot per energy bin
for j in range(n_bins):
    plot, fig = pf.new_plot(r'$E_j\;\;\mathrm{[GeV]}$', r'$\frac{E_i^2 E_j^2}{\Delta E_i\Delta E_j} C_p(E_i,E_j)\;\;\mathrm{[GeV^2 cm^{-4} s^{-2} sr^{-1}]}$', 'log', 'log', label_size=18)
    Emin1     = data__ebins[j,0]
    Emax1     = data__ebins[j,1]
    E_1       = np.sqrt(Emin1*Emax1)
    Emin2     = data__ebins[:,0]
    Emax2     = data__ebins[:,1]
    E_2       = np.sqrt(Emin2*Emax2)
    DeltaE1   = Emax1 - Emin1
    DeltaE2   = Emax2 - Emin2
    fac       = E_1**2*E_2**2/DeltaE1/DeltaE2
    xerr      = [E_2-Emin2, Emax2-E_2 ]

    Cp_simple               = AC.get_Cp_simple_model                   ( Emin1, Emax1, Emin2, Emax2,  A, beta, Gamma,                          k=1.  )
    Cp_simple_cutoff        = AC.get_Cp_simple_cutoff_model            ( Emin1, Emax1, Emin2, Emax2,  A, beta, Gamma,                 E_b,     k=1.  )
    Cp_simple_Gamma         = AC.get_Cp_simple_model_GammaDistr        ( Emin1, Emax1, Emin2, Emax2,  A, beta, mu=Gamma, sigma=sigma,          k=1.  )
    Cp_simple_cutoff_Gamma  = AC.get_Cp_simple_cutoff_model_GammaDistr ( Emin1, Emax1, Emin2, Emax2,  A, beta, mu=Gamma, sigma=sigma, E_b=E_b, k=1.  )


    #print( np.amax( np.fabs( 1 - Cp_simple/Cp_simple_Gamma ) )               )
    #print( np.amax( np.fabs( 1 - Cp_simple_cutoff/Cp_simple_cutoff_Gamma ) ) )

    plot.errorbar( E_2 , fac * Cp_simple             ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple'                     )
    plot.errorbar( E_2 , fac * Cp_simple_cutoff      ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff'              )
    plot.errorbar( E_2 , fac * Cp_simple_Gamma       ,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple Gamma distr.'        )
    plot.errorbar( E_2 , fac * Cp_simple_cutoff_Gamma,   lw=2, yerr=0, xerr=xerr, fmt='o', markersize=0, label='Cp simple cutoff Gamma distr.' )

    #plot.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, fontsize=15)
    plot.legend(frameon=False, loc='lower left', bbox_to_anchor=(0.05, 0.05), ncol=1, fontsize=15)
    plt.savefig('test_Cp_crossE_data_%i.pdf' % j )


