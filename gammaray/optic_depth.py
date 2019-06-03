import numpy as np
import os

CDIR = os.path.dirname(os.path.realpath(__file__))

dominguez11     = np.genfromtxt( CDIR+'/data/tau_dominguez11.out.txt' )

dominguez11_tau = np.zeros( (len(dominguez11[:,0])-1,len(dominguez11[0,:]) )  )
dominguez11_E   = np.zeros( len( dominguez11[:,0] )-1 )
dominguez11_z   = np.zeros( len( dominguez11[0,:] )   )


dominguez11_tau[:,1:]  = dominguez11[1:,1:]
dominguez11_E[ :]      = dominguez11[1:,0]
dominguez11_z[1:]      = dominguez11[0,1:]

dominguez11_z[0] = 0.0

from scipy import interpolate

f_dominguez11_E_z = interpolate.interp2d( dominguez11_E, dominguez11_z, np.exp(-dominguez11_tau).transpose() )

'''
    Optic depth as provided by
    
    param E    Energy in GeV
    prarm z    Redshift distance
    param h    (optional) H_0/100km/s/Mpc
    
    return optic depth tau(
    '''
def _tau_dominguez11( E, z, h=0.7, verbosity=1 ):
    if verbosity>0:
        if h!=0.7:
            print('Warning: tau_dominguez11 is only implemented for 0.7 at this moment.')
        if z<0.00 or z>2.0 or E<0.1 or E>30000:
            print('Warning: tau_dominguez11 is out of tabulated boundary, E/GeV is in [%.2e,%.2e] and z is in [%.3f, %.3f], use extrapolation' %(30, 30000, 0.01, 2.))
    if z>2.0:
        E = E * (z/2.0)**(1./1.3)
        z = 2.0
    if E<30.0:
        z = z * (30./E)**(-1.3)
        E = 30.
    E = E/1e3
    return -np.log(f_dominguez11_E_z(E,z)[0])
tau_dominguez11 = np.vectorize(_tau_dominguez11)

#print( tau_dominguez11(30, 0.5)  )
#print( tau_dominguez11(100, 1.5) )
#print( tau_dominguez11(100, 2.2) )


#import matplotlib                       as      mpl
#mpl.use('Agg')
#import matplotlib.pyplot                as      plt
#from matplotlib.lines                   import  Line2D
#
#import MKastro.basic.PlotFunctions      as      pf
#
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
#plt.savefig('tau_extrapolation.png')
