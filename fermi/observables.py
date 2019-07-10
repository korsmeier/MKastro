#!/usr/bin/env python

import numpy  as np
import os
import sys
import math
from   scipy  import interpolate

import MKastro.fermi.cat4FGL as cat4FGL
import MKastro.fermi.cat4LAC as cat4LAC

CDIR = os.path.dirname(os.path.realpath(__file__))

_4pi_to_square_deg = 41252.9612
_sr_to_square_deg  =  3282.8


n_4FGL = cat4FGL.n()


_4FGL_Source_name       = cat4FGL.get( 'Source_Name'        )
_4FGL_GLON              = cat4FGL.get( 'GLON'               )
_4FGL_GLAT              = cat4FGL.get( 'GLAT'               )
_4FGL_Flux1000          = cat4FGL.get( 'Flux1000'           )
_4FGL_Unc_Flux1000      = cat4FGL.get( 'Unc_Flux1000'       )
_4FGL_gamma             = cat4FGL.get( 'PL_Index'           )
_4FGL_CLASS1            = cat4FGL.get( 'CLASS1'             )

_4FGL_z                 = -2. * np.ones( n_4FGL )
_4FGL_SED               = np.empty( n_4FGL, dtype='S5' )

_4LAC_z                 = cat4LAC.get('Redshift'    )
_4LAC_SED               = cat4LAC.get('SED_class'   )


for i in range(n_4FGL):
    i4LAC = cat4LAC.find_source_from4FGL( _4FGL_Source_name[i] )
    if i4LAC>=0:
        _4FGL_z  [i] = _4LAC_z  [i4LAC]
        _4FGL_SED[i] = _4LAC_SED[i4LAC]
        if _4FGL_z[i]<0:
            _4FGL_z[i] = -1

def get_4FGL_4LAC( observable, row=-1  ):
    if observable=='Redshift':
        if row<0:
            return np.copy(_4FGL_z[:])
        return _4FGL_z[row]
    if observable=='SED_class':
        if row<0:
            return np.copy(_4LAC_SED[:])
        return _4LAC_SED[row]
    return cat4FGL.get( observable, row )

'''
    Function to apply cuts on 4FGL (and 4LAC) sources
    
    param z      len(2) array [ z_lower, z_upper ], z=-1: no redshift measurement, z=-2: not in 4LAC
    param bLAT   minimal Latitude in deg
    param CLASS  string with source classes, use 'all' for any class. Example "fsrq bll bcu"
    param CLASS  string with source SED, use 'all' for any SED. Example "HSP LSP"
    
    return      list of 4FGL indices fullfilling the cuts
    '''
def get_list_cuts_4FGL_4LAC( z=[-3, 10], Gamma=[0,5], bLAT=30, CLASS='all', SED='all' ):
    index_list = []
    for i in range(n_4FGL):
        if not         _4FGL_z[i]                            >   z[0]       : continue
        if not         _4FGL_z[i]                            <   z[1]       : continue
        if not         _4FGL_gamma[i]                        >   Gamma[0]   : continue
        if not         _4FGL_gamma[i]                        <   Gamma[1]   : continue
        if not np.fabs(_4FGL_GLAT[i])                        >   bLAT       : continue
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
def _F_E__from__F_Efrom( F_Efrom, Gamma, Emin, Emax, Emin_from, Emax_from ):
    if np.fabs(-Gamma+1)>1e-5:
        N = np.power(Emax,      -Gamma+1) - np.power(Emin,      -Gamma+1)
        D = np.power(Emax_from, -Gamma+1) - np.power(Emin_from, -Gamma+1)
    else:
        N = np.log(1000.)
        D = np.log(Emax/Emin)
    return F_Efrom * N/D
F_E__from__F_Efrom = np.vectorize(_F_E__from__F_Efrom)

#'''
#    Function to calculate flux in a specific energy bin.
#    Assume powerlaw behaviour.
#    '''
#def _F_E__from__F_1000( F_1000, Gamma, Emin, Emax ):
#    if np.fabs(-Gamma+1)>1e-5:
#        N = np.power(Emax, -Gamma+1) - np.power(Emin, -Gamma+1)
#        D = np.power(100., -Gamma+1) - np.power( 1.0, -Gamma+1)
#    else:
#        N = np.log(1000.)
#        D = np.log(Emax/Emin)
#    return F_1000 * N/D
def F_E__from__F_1000( F_1000, Gamma, Emin, Emax ):
    return F_E__from__F_Efrom( F_1000, Gamma, Emin, Emax, 1.0, 100.0 )


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
def get_dNdF( F,  E=[0.1,100], z=[-3, 10], Gamma=[0,5], bLAT=30, CLASS='all', SED='all'  ):
    global _4FGL_gamma, _4FGL_Flux1000
    #
    sources_list = get_list_cuts_4FGL_4LAC( z, Gamma, bLAT, CLASS, SED )
    #
    F_list  = F_E__from__F_1000( _4FGL_Flux1000, _4FGL_gamma, E[0], E[1] )[sources_list]
    F_upper = F[1:  ]
    F_lower = F[ :-1]
    if len(F_upper)==1:
        F_upper = F_upper[0]
        F_lower = F_lower[0]
    N       = np.histogram(F_list, F)[0]
    N_err   = np.sqrt(np.maximum(1., N))
    F_mean  = np.sqrt( F_upper * F_lower )
    Delta_F = F_upper - F_lower
    geom_fact = 1 + np.cos(math.pi * (0.5+bLAT/180.))
    dN_dF     = N     * 1./ (Delta_F*_4pi_to_square_deg*geom_fact)
    dN_dF_err = N_err * 1./ (Delta_F*_4pi_to_square_deg*geom_fact)
    #
    return F_mean, dN_dF, dN_dF_err


index_4FGL_F_thresh = 11
d_F_thresh__Gamma = np.genfromtxt(CDIR+'/data/8yr_flux_threshold.dat')
d_F_thresh__Gamma_4FGL = d_F_thresh__Gamma[:index_4FGL_F_thresh,:]
d_F_thresh__Gamma_3FHL = d_F_thresh__Gamma[index_4FGL_F_thresh:,:]
i_F_thresh__Gamma_4FGL = interpolate.interp1d( d_F_thresh__Gamma_4FGL[:,1], np.log(d_F_thresh__Gamma_4FGL[:,0]), fill_value='extrapolate' )
def F_thresh__Gamma_4FGL( Gamma ):
    return np.exp(i_F_thresh__Gamma_4FGL(Gamma) )
i_F_thresh__Gamma_3FHL = interpolate.interp1d( d_F_thresh__Gamma_3FHL[:,1], np.log(d_F_thresh__Gamma_3FHL[:,0]), fill_value='extrapolate' )
def F_thresh__Gamma_3FHL( Gamma ):
    return np.exp(i_F_thresh__Gamma_3FHL(Gamma) )



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
#
#
#plot, fig = pf.new_plot( r'$\Gamma$', r'$S \;\mathrm{[cm^{-2}s^{-1}]}$','linear', 'log')
#vec_Gamma = np.arange(0,5.01,0.1)
#plot.plot(d_F_thresh__Gamma_4FGL[:,1],d_F_thresh__Gamma_4FGL[:,0],label='4FGL', color=pf.colors[0])
#plot.plot(vec_Gamma, F_thresh__Gamma_4FGL(vec_Gamma),label='4FGL', color=pf.colors[0], dashes=(5,5))
#plt.savefig('Scut_Gamma_4FGL.png')
#plot.plot(d_F_thresh__Gamma_3FHL[:,1],d_F_thresh__Gamma_3FHL[:,0],label='3FHL', color=pf.colors[4])
#plot.plot(vec_Gamma, F_thresh__Gamma_3FHL(vec_Gamma),label='3FHL', color=pf.colors[4], dashes=(5,5))
#plt.legend()
#plt.savefig('Scut_Gamma.png')
#
