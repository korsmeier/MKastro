#!/usr/bin/env python

import numpy as np
import os
import sys
import math

import MKastro.fermi.cat4FGL as cat4FGL
import MKastro.fermi.cat4LAC as cat4LAC

_4pi_to_square_deg = 41252.9612
_sr_to_square_deg = 3282.8


n_4FGL = cat4FGL.n()


_4FGL_Source_name       = cat4FGL.get( 'Source_Name'        )
_4FGL_GLON              = cat4FGL.get( 'GLON'               )
_4FGL_GLAT              = cat4FGL.get( 'GLAT'               )
_4FGL_Flux1000          = cat4FGL.get( 'Flux1000'           )
_4FGL_Unc_Flux1000      = cat4FGL.get( 'Unc_Flux1000'       )

_4FGL_gamma             = cat4FGL.get( 'PL_Index'           )

_4FGL_z                 = -2. * np.ones( n_4FGL )

_4FGL_SED               = np.empty( n_4FGL, dtype='S5' )

_4LAC_z         = cat4LAC.get('Redshift'    )
_4LAC_SED       = cat4LAC.get('SED_class'   )




for i in range(n_4FGL):
    i4LAC = cat4LAC.find_source_from4FGL( _4FGL_Source_name[i] )
    if i4LAC>=0:
        _4FGL_z  [i] = _4LAC_z  [i4LAC]
        _4FGL_SED[i] = _4LAC_SED[i4LAC]
        if _4FGL_z[i]<0:
            _4FGL_z[i] = -1

'''
    Function to apply cuts on 4FGL (and 4LAC) sources
    
    param z      len(2) array [ z_lower, z_upper ], z=-1: no redshift measurement, z=-2: not in 4LAC
    param bLAT   minimal Latitude in deg
    param CLASS  string with source classes, use 'all' for any class. Example "fsrq bll bcu"
    param CLASS  string with source SED, use 'all' for any SED. Example "HSP LSP"
    
    return      list of 4FGL indices fullfilling the cuts
    '''
def get_list_cuts_4FGL_4LAC( z=[-3, 10], bLAT=30, CLASS='all', SED='all' ):
    index_list = []
    for i in range(n_4FGL):
        if not         _4FGL_z[i]                            >   z[0]       :
            continue
        if not         _4FGL_z[i]                            <   z[1]       :
            continue
        if not np.fabs(_4FGL_GLAT[i])                        >   bLAT       :
            continue
        if not 'all' in CLASS:
            if not     str(_4FGL_CLASS1[i]).lower().split(' ')[0]       in  CLASS.lower()               :  continue
        if not 'all' in SED:
            if not     (str(_4FGL_SED[i].decode('utf-8'))+' ').upper().split(' ')[0] in  SED.upper()    :  continue
        index_list.append(i)
    index_list = np.array(index_list)
    return index_list


'''
    Function to calculate flux in a specific energy bin.
    Assume powerlaw behaviour.
    '''
def _F_E__from__F_1000( F_1000, Gamma, Emin, Emax ):
    if np.fabs(-Gamma+1)>1e-5:
        N = np.power(Emax, -Gamma+1) - np.power(Emin, -Gamma+1)
        D = np.power(100., -Gamma+1) - np.power( 1.0, -Gamma+1)
    else:
        N = np.log(1000.)
        D = np.log(Emax/Emin)
    return F_1000 * N/D
F_E__from__F_1000 = np.vectorize(_F_E__from__F_1000)

'''
    Function to apply cuts on 4FGL (and 4LAC) sources
    
    param F      array, F bins
    param E      len(2) array [ E_lower, E_upper ], E in GeV
    param z      len(2) array [ z_lower, z_upper ], z=-1: no redshift measurement, z=-2: not in 4LAC
    param bLAT   minimal Latitude in deg
    param CLASS  string with source classes, use 'all' for any class. Example "fsrq bll bcu"
    param CLASS  string with source SED, use 'all' for any SED. Example "HSP LSP"
    
    return       F_average dN/dF   ( arrays of size len(F)-1 )
    '''
def get_dNdF( F,  E=[0.1,100], z=[-3, 10], bLAT=30, CLASS='all', SED='all'  ):
    global _4FGL_gamma, _4FGL_Flux1000
    #
    sources_list = get_list_cuts_4FGL_4LAC( z, bLAT, CLASS, SED )
    #
    F_list  = F_E__from__F_1000( _4FGL_Flux1000, _4FGL_gamma, E[0], E[1] )[sources_list]
    F_upper = F[1:  ]
    F_lower = F[ :-1]
    N       = np.histogram(F_list, F)[0]
    N_err   = np.sqrt(np.maximum(1., N))
    F_mean  = np.sqrt( F_upper * F_lower )
    Delta_F = F_upper - F_lower
    geom_fact = 1 + np.cos(math.pi * (0.5+bLAT/180.))
    dN_dF     = N     * 1./ (Delta_F*_4pi_to_square_deg*geom_fact)
    dN_dF_err = N_err * 1./ (Delta_F*_4pi_to_square_deg*geom_fact)
    #
    return F_mean, dN_dF, dN_dF_err



#import matplotlib                       as      mpl
#mpl.use('Agg')
#import matplotlib.pyplot                as      plt
#from matplotlib.lines                   import  Line2D
#
#import MKastro.basic.PlotFunctions      as      pf
#import MKastro.basic.Constants          as      cst
#import MKastro.basic.Cosmology          as      cosmo
#import MKastro.basic.Transformation     as      tr
# test:
#E_l, E_u = [1,2]
#plot, fig = pf.new_plot(r'$F$', '$dN/dF \;\; \mathrm{[GeV]}$', 'log', 'log')
#F_bins = np.power(10, np.linspace(-13, -7, 20))
#F, dNdF, dNdF_unc = get_dNdF( F_bins, E=[E_l, E_u], z=[-3, 10], bLAT=10, CLASS='all', SED='all' )
#plot.errorbar(F, dNdF*F**2, yerr=dNdF_unc*F**2,
#              xerr=0,
#              color=pf.colors[0],
#              fmt='o',
#              label='$E=%.1f$ to $%.1f GeV$, bcut = $30$' % (E_l, E_u)
#             )
#
#F, dNdF, dNdF_unc = get_dNdF( F_bins, E=[E_l, E_u], z=[-3, 10], bLAT=30, CLASS='all', SED='all' )
#plot.errorbar(F, dNdF*F**2, yerr=dNdF_unc*F**2,
#              xerr=0,
#              color=pf.colors[4],
#              fmt='o',
#              label='$E=%.1f$ to $%.1f GeV$, bcut = $10$' % (E_l, E_u)
#              )
#
#
##leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('dNdF_test.pdf')