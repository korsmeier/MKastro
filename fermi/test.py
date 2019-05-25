#!/usr/bin/env python

import numpy as np

import cat3LAC
import cat4FGL

cat3LAC.print_columns()
cat4FGL.print_columns()

n_4FGL = cat4FGL.n()

_4FGL_RA                = cat4FGL.get( 'RAJ2000'            )
_4FGL_DE                = cat4FGL.get( 'DEJ2000'            )
_4FGL_GLON              = cat4FGL.get( 'GLON'               )
_4FGL_GLAT              = cat4FGL.get( 'GLAT'               )
_4FGL_Flux_Band         = cat4FGL.get( 'Flux_Band'          )
_4FGL_Unc_Flux_Band     = cat4FGL.get( 'Unc_Flux_Band'      )
_4FGL_CLASS1            = cat4FGL.get( 'CLASS1'             )
_4FGL_AngSep            = cat4FGL.get( 'Conf_95_SemiMajor'  )

_4FGL_z                 = -1. * np.ones( len(_4FGL_RA) )
_4FGL_SED               = -1. * np.ones( len(_4FGL_RA) )

Ebins_4FGL = cat4FGL.get_energy_bounds()

print( (_4FGL_Flux_Band    [0]) )
print( (_4FGL_Unc_Flux_Band[0]) )
print( (_4FGL_Flux_Band    [1]) )
print( (_4FGL_Unc_Flux_Band[1]) )


#_4FGL_CLASS2    = cat4FGL.get('CLASS2'      )
#print(_4FGL_CLASS1[::100])
#print(_4FGL_AngSep[::100])

_3LAC_z         = cat3LAC.get('Redshift'    )
_3LAC_SED       = cat3LAC.get('SED Class'   )


c = 0
for i in range(n_4FGL):
    i3LAC = cat3LAC.find_source( _4FGL_RA[i], _4FGL_DE[i], 2*_4FGL_AngSep[i] )
    if i3LAC>=0:
        _4FGL_z  [i] = _3LAC_z[i3LAC]
        _4FGL_SED[i] = _3LAC_z[i3LAC]
        c += 1
print( c )


import matplotlib               as mpl
mpl.use('Agg')
import matplotlib.pyplot        as plt
from matplotlib.lines           import Line2D
import MKastro.PlotFunctions    as pf

import MKastro.Constants        as cst
import MKastro.Cosmology        as cosmo
import MKastro.Transformation   as tr
#import MKastro.DM_spectrum      as dm


## Plot some
plot, fig = pf.new_plot(r'$E\;\;\mathrm{[MeV]}$', r'$E^2\; F \;\;\mathrm{[MeV^2 cm^{-2} s^{-1}]}$', 'log', 'log')

for i in range(n_4FGL):
    if not _4FGL_z > 0:
        continue
    if not _4FGL_CLASS1.upper() == 'BBL':
        continue

    E    = np.sqrt(Ebins_4FGL[:,0]*Ebins_4FGL[:,1] )
    E_err= np.copy(Ebins_4FGL)

    E_err[:,0] =-Ebins_4FGL[:,0]+E
    E_err[:,1] = Ebins_4FGL[:,1]-E

    F    = _4FGL_Flux_Band    [i]
    F_err= _4FGL_Unc_Flux_Band[i]

    plot.errorbar(  E,  E**2*F, xerr=[E_sl,E_su], yerr=0, color=pf.colors[-1], fmt='' )

    break

#plot.set_ylim( (2e-16, 5e-10) )
#leg1 = plot.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.03, 0.97), ncol=2, fontsize=legend_fontsize)
plt.savefig('test.pddf')






