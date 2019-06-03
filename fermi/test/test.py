#!/usr/bin/env python

import numpy as np
import os
import sys

import matplotlib                       as      mpl
mpl.use('Agg')
import matplotlib.pyplot                as      plt
from matplotlib.lines                   import  Line2D

import MKastro.basic.PlotFunctions      as      pf
import MKastro.basic.Constants          as      cst
import MKastro.basic.Cosmology          as      cosmo
import MKastro.basic.Transformation     as      tr

import MKastro.gammaray.optic_depth     as      optic_depth


#plot, fig = pf.new_plot(r'$E\;\;[GeV]$', r'$z$', 'log', 'log')
#plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.15)
#_z = np.linspace( 0.001, 100,500 )
#_E = np.power(10, np.linspace( 0.0, 5.0, 100 ) )
#E, z = np.meshgrid(_E, _z, indexing='ij' )
#tau = optic_depth.tau_dominguez11( E, z, verbosity=0 )
#CS = plot.contourf( E, z, np.exp(-tau), 10, cmap='magma' )
#cbar = plt.colorbar(CS)
#y = (_E/100.)**-1.3
#plot.plot( _E, y, color='red' )
#rtg = mpl.patches.Rectangle( (30, 0.01) , width=30000-30, height=2-0.01, ec='green', fc=(0,0,0,0) )
#plot.add_artist(rtg)
#cbar.ax.set_ylabel(r'$\exp(-\tau)$')
#plt.savefig('tau.png')
#
#sys.exit(0)

import MKastro.fermi.cat4FGL as cat4FGL
import MKastro.fermi.cat4LAC as cat4LAC


#cat4FGL.print_columns()
#cat4LAC.print_columns()

n_4FGL = cat4FGL.n()

_4FGL_Source_name       = cat4FGL.get( 'Source_Name'        )
_4FGL_RA                = cat4FGL.get( 'RAJ2000'            )
_4FGL_DE                = cat4FGL.get( 'DEJ2000'            )
_4FGL_GLON              = cat4FGL.get( 'GLON'               )
_4FGL_GLAT              = cat4FGL.get( 'GLAT'               )
_4FGL_Flux_Band         = cat4FGL.get( 'Flux_Band'          )
_4FGL_Unc_Flux_Band     = cat4FGL.get( 'Unc_Flux_Band'      )
_4FGL_CLASS1            = cat4FGL.get( 'CLASS1'             )
_4FGL_AngSep            = cat4FGL.get( 'Conf_95_SemiMajor'  )

_4FGL_gamma             = cat4FGL.get( 'PL_Index'      )
_4FGL_alpha             = cat4FGL.get( 'LP_Index'      )
_4FGL_beta              = cat4FGL.get( 'LP_beta'       )
_4FGL_alpha_unc         = cat4FGL.get( 'Unc_LP_Index'  )
_4FGL_beta_unc          = cat4FGL.get( 'Unc_LP_beta'   )
_4FGL_E0                = cat4FGL.get( 'Pivot_Energy'   )



_4FGL_z                 = -1. * np.ones( len(_4FGL_RA) )
_4FGL_SED               = np.empty( len(_4FGL_RA), dtype='S5' )

Ebins_4FGL = cat4FGL.get_energy_bounds()


_4LAC_z         = cat4LAC.get('Redshift'    )
_4LAC_SED       = cat4LAC.get('SED_class'   )


c = 0
for i in range(n_4FGL):
    i4LAC = cat4LAC.find_source_from4FGL( _4FGL_Source_name[i] )
    if i4LAC>=0:
        _4FGL_z  [i] = _4LAC_z  [i4LAC]
        _4FGL_SED[i] = _4LAC_SED[i4LAC]
        c += 1
print( c )



def get_list_cuts( z=[0.00001, 5], bLAT=30, CLASS='all', SED='all' ):
    index_list = []
    for i in range(n_4FGL):
        if not         _4FGL_z[i]                            >   z[0]       :  continue
        if not         _4FGL_z[i]                            <   z[1]       :  continue
        if not np.fabs(_4FGL_GLAT[i])                        >   bLAT       :  continue
        if not 'all' in CLASS:
            #print(str(_4FGL_CLASS1[i]).lower().split(' ')[0])
            if not     str(_4FGL_CLASS1[i]).lower().split(' ')[0]       in  CLASS      :  continue
        if not 'all' in SED:
            #print( (str(_4FGL_SED   [i].decode('utf-8'))+' ').upper().split(' ')[0] )
            if not     (str(_4FGL_SED[i].decode('utf-8'))+' ').upper().split(' ')[0] in  SED      :  continue
        index_list.append(i)
    index_list = np.array(index_list)
    return index_list


#bll  = get_list_cuts(CLASS='bll' )
#bcu  = get_list_cuts(CLASS='bcu' )
#fsrq = get_list_cuts(CLASS='fsrq')
#
#
#plot, fig = pf.new_plot(r'$\alpha$', r'$\beta$', 'linear', 'linear')
#plot.scatter( _4FGL_alpha[bll ], _4FGL_beta[bll ], c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_alpha[bcu ], _4FGL_beta[bcu ], c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_alpha[fsrq], _4FGL_beta[fsrq], c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_alpha_beta.pdf')
#
#
#plot, fig = pf.new_plot(r'$\Gamma$', r'$\alpha$', 'linear', 'linear')
#plot.scatter( _4FGL_gamma[bll ], _4FGL_alpha[bll ], c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_gamma[bcu ], _4FGL_alpha[bcu ], c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_gamma[fsrq], _4FGL_alpha[fsrq], c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_gamma_alpha.pdf')
#
#
#plot, fig = pf.new_plot(r'$\Gamma$', r'$\beta$', 'linear', 'linear')
#plot.scatter( _4FGL_gamma[bll ], _4FGL_beta[bll ], c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_gamma[bcu ], _4FGL_beta[bcu ], c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_gamma[fsrq], _4FGL_beta[fsrq], c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_gamma_beta.pdf')
#
#
#plot, fig = pf.new_plot(r'$\Gamma$', '$E_0 \;\; \mathrm{[GeV]}$', 'linear', 'log')
#plot.scatter( _4FGL_gamma[bll ], _4FGL_E0[bll ]/1000, c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_gamma[bcu ], _4FGL_E0[bcu ]/1000, c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_gamma[fsrq], _4FGL_E0[fsrq]/1000, c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_gamma_E0.pdf')
#
#
#plot, fig = pf.new_plot(r'$\Gamma$', '$E_0 (1+z) \;\; \mathrm{[GeV]}$', 'linear', 'log')
#plot.scatter( _4FGL_gamma[bll ], _4FGL_E0[bll ]/1000*(1+_4FGL_z[bll ]), c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_gamma[bcu ], _4FGL_E0[bcu ]/1000*(1+_4FGL_z[bcu ]), c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_gamma[fsrq], _4FGL_E0[fsrq]/1000*(1+_4FGL_z[fsrq]), c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_gamma_E0z.pdf')
#
#plot, fig = pf.new_plot(r'$E_0 \;\; \mathrm{[GeV]}$', r'$1+z$', 'log', 'log')
#plot.scatter( _4FGL_E0[bll ], 1+_4FGL_z[bll ], c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_E0[bcu ], 1+_4FGL_z[bcu ], c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_E0[fsrq], 1+_4FGL_z[fsrq], c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_z_E0.pdf')
#
#plot, fig = pf.new_plot(r'$\alpha$', r'$1+z$', 'linear', 'log')
#plot.scatter( _4FGL_alpha[bll ], 1+_4FGL_z[bll ], c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_alpha[bcu ], 1+_4FGL_z[bcu ], c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_alpha[fsrq], 1+_4FGL_z[fsrq], c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_z_alpha.pdf')
#
#
#plot, fig = pf.new_plot(r'$\beta$', r'$1+z$', 'linear', 'log')
#plot.scatter( _4FGL_beta[bll ], 1+_4FGL_z[bll ], c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_beta[bcu ], 1+_4FGL_z[bcu ], c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_beta[fsrq], 1+_4FGL_z[fsrq], c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_z_beta.pdf')
#
#plot, fig = pf.new_plot(r'$\gamma$', r'$1+z$', 'linear', 'log')
#plot.scatter( _4FGL_gamma[bll ], 1+_4FGL_z[bll ], c=pf.colors[2* 0], label='bll ' )
#plot.scatter( _4FGL_gamma[bcu ], 1+_4FGL_z[bcu ], c=pf.colors[2* 1], label='bcu ' )
#plot.scatter( _4FGL_gamma[fsrq], 1+_4FGL_z[fsrq], c=pf.colors[2* 2], label='fsrq' )
#leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
#plt.savefig('corr_z_gamma.pdf')


hsp   = get_list_cuts(CLASS='bll', SED='HSP'     )

print('hsp')
print(len(hsp))

lisp  = get_list_cuts(CLASS='bll', SED='LSP ISP' )


def correct_F_unc(F, F_unc):
    for i in range(len(F_unc[:,0])):
        if F_unc[i,0]!=F_unc[i,0]:
            F_unc[i,0] = F[i]*(1.-1e-7)


data_PL_x     = np.ones(7)
data_PL_y     = np.ones(7)
data_PL_y_unc = np.ones(7)
def chiSq_PL( par ):
    log10A = par[0]
    Gamma  = par[1]
    model_PL_y = 10**log10A * np.power( data_PL_x, -Gamma )
    return (np.power( (model_PL_y-data_PL_y)/data_PL_y_unc, 2 )).sum()

# prepare minuit fit



## Loop over
_4FGL_GammaIntr_HSP = []


from iminuit import Minuit
my_minuit = Minuit.from_array_func(chiSq_PL, (-6, 2),
                                   error=(1, 0.5),
                                   fix=(False, False),
                                   limit=((-10,0), (0.5, 3.5)),
                                   name = ('log(A)', 'Gamma'),
                                   errordef=1,
                                   print_level=0)
my_minuit.get_param_states()


for j in range(len(hsp)):
    print( '%4i/%i' % (j, len(hsp)) )
    i = hsp[j]
    
    z  = _4FGL_z[i]
    dL = cosmo.d_Lum__z(z)
    
    E0    = np.sqrt(Ebins_4FGL[:,0]*Ebins_4FGL[:,1] )
    E0_unc= np.copy(Ebins_4FGL)
    E0_unc[:,0] =-Ebins_4FGL[:,0]+E0
    E0_unc[:,1] = Ebins_4FGL[:,1]-E0
    
    DE0 = E0_unc[:,1]-E0_unc[:,0]

    E = E0         * (1+z)
    E_unc = E0_unc * (1+z)


    tau = optic_depth.tau_dominguez11( E/1000., z, verbosity=0 )

    F0    = _4FGL_Flux_Band[i]     * np.exp(tau)   /  DE0
    F0_unc= _4FGL_Unc_Flux_Band[i]
    F0_unc[:,0] *=                  -np.exp(tau)   /  DE0
    F0_unc[:,1] *=                   np.exp(tau)   /  DE0

    F     = F0     * (1+z)
    F_unc = F0_unc * (1+z)

    correct_F_unc(F, F_unc)

    E2F_unc = np.copy(F_unc)
    E2F_unc[:,0] *= E**2
    E2F_unc[:,1] *= E**2
    
    
    #plot.errorbar(  E,  E**2*F, xerr=E_unc.transpose(), yerr=E2F_unc.transpose(), color=pf.cmap_rainbow(j/10.), fmt='.' )
    #plot.errorbar(  E,  E**2*F * np.exp(tau), xerr=(1+z)*E0_unc.transpose(), yerr=E2F_unc.transpose(), color=pf.cmap_rainbow(np.amin([2.*z, 1])), fmt='.' )
    plot, fig = pf.new_plot(r'$E\;\;\mathrm{[MeV]}$', r'$E^2\; F \;\;\mathrm{[MeV cm^{-2} s^{-1}]}$', 'log', 'log')
    plot.errorbar(  E,  E**2*F, xerr=E_unc.transpose(), yerr=E2F_unc.transpose(), color=pf.colors[2], fmt='o' )
    
    data_PL_x    [:] = E[:]
    data_PL_y    [:] = F[:]
    data_PL_y_unc[:] = 0.5*(F_unc[:,0] + F_unc[:,1])
    my_minuit.migrad()
    my_minuit.hesse()
#    print( my_minuit.values )
#    print( my_minuit.errors )
    x = 10**np.linspace(1, 6, 100)
    y = 10**my_minuit.values[0] * np.power(x,2.-my_minuit.values[1])
    plot.plot(  x, y, color=pf.colors[-1] )
    
    _4FGL_GammaIntr_HSP.append(my_minuit.values[1])
    
    plot.set_ylim( (np.amin(y), np.amax(y)) )
    plt.savefig( 'hsp_%05i.pdf' % j )

#    if j>5:
#        break

print(_4FGL_GammaIntr_HSP)

_4FGL_GammaIntr_HSP = np.array(_4FGL_GammaIntr_HSP)

plot, fig = pf.new_plot(r'$\Gamma$ 4FGL', r'my $\Gamma$', 'linear', 'linear')
plot.scatter( _4FGL_gamma[hsp ], _4FGL_GammaIntr_HSP, c=pf.colors[2* 0], label='HSP' )
leg1 = plot.legend(frameon=False, loc='upper right', ncol=1, fontsize=15)
plt.savefig('hsp_corr_gamma.pdf')

sys.exit(0)
#plot.set_ylim( 1e5, 1e8 )
#plot.set_ylim( 1e-10, 1e-4 )
#leg1 = plot.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.03, 0.97), ncol=2, fontsize=legend_fontsize)



plot, fig = pf.new_plot(r'$E\;\;\mathrm{[MeV]}$', r'$E\; F \;\;\mathrm{[MeV cm^{-2} s^{-1}]}$', 'log', 'log')
for j in range(len(lisp)):
    
    i = lisp[j]
    
    z  = _4FGL_z[i]
    dL = cosmo.d_Lum__z(z)
    
    if z<1e-4:
        continue

    #print(z)

    E0    = np.sqrt(Ebins_4FGL[:,0]*Ebins_4FGL[:,1] )

    E0_unc= np.copy(Ebins_4FGL)
    E0_unc[:,0] =-Ebins_4FGL[:,0]+E0
    E0_unc[:,1] = Ebins_4FGL[:,1]-E0

    DE0 = E0_unc[:,1]-E0_unc[:,0]
    
    F    = _4FGL_Flux_Band[i]          / DE0
    F_unc= _4FGL_Unc_Flux_Band[i]
    
    F_unc[:,0] /= -DE0
    F_unc[:,1] /=  DE0
    
    E = E0*(1+z)
    tau = optic_depth.tau_dominguez11( E/1000., z, verbosity=0 )
    
    fact = 1./F[3]/(1+z)**2
    
    F    *= fact
    F_unc*= fact
    correct_F_unc(F, F_unc)
    
    
    E2F_unc = np.copy(F_unc)
    E2F_unc[:,0] *= E**2 * np.exp(tau)
    E2F_unc[:,1] *= E**2 * np.exp(tau)
    
    
    #plot.errorbar(  E,  E**2*F, xerr=E_unc.transpose(), yerr=E2F_unc.transpose(), color=pf.colors[-1], fmt='.' )
    plot.errorbar(  E,  E**2*F * np.exp(tau), xerr=(1+z)*E0_unc.transpose(), yerr=E2F_unc.transpose(), color=pf.cmap_rainbow(np.amin([2.*z, 1])), fmt='.' )


plot.set_ylim( 1e5, 1e8 )
#leg1 = plot.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.03, 0.97), ncol=2, fontsize=legend_fontsize)
plt.savefig('lisp.pdf')


