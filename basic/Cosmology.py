import numpy as np
from scipy import integrate

import MKastro.basic.Constants as cst

c       = cst.c['m/s']                                                      # in m/s
H_0     = cst.H_0['km/s/Mpc'] * cst.length['km->m'] / cst.length['Mpc->m']  # in s
Omega_M = cst.Omega_M
Omega_R = cst.Omega_R

# Assume flat universe
Omega_L = 1. - Omega_M - Omega_R

'''
    Hubble constant as function of z in LCDM (flat universe)
    
    param  z
    return Hubble constant
    '''
def H__z(z): # in s
    global H_0, Omega_R, Omega_M, Omega_L
    return H_0 * np.sqrt( Omega_R*(1+z)**4 + Omega_M*(1+z)**3 + Omega_L )

'''
    Luminosity function as function of z in LCDM (flat universe)
    
    param  z
    return luminosity disatnce
    '''
def _d_Lum__z(z): # in m
    i = integrate.quad( lambda _z: 1./H__z(_z) , 0, z)
    if i[0] != 0:
        if i[1]/i[0] > 0.01:
            print('Warning. Error in d_Lum__z is larger that 1 percent for z = %f' % z )
    return i[0] * c * (1+z)
d_Lum__z = np.vectorize(_d_Lum__z)

#z = 0.001
#print( 'H_0      = %.3e s^-1' % H_0         )
#print( 'z        = %.3f     ' % z           )
#print( 'H__z     = %.3e s   ' % H__z(z)     )
#print( 'd_Lum__z = %.3e Mpc ' % (d_Lum__z(z)*cst.length['m->Mpc']) )

