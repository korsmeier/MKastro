
__author__='M. Korsmeier - (korsmeier@physik.rwth-aachen.de)'
__version__='1.0'

import os
import numpy as np

CDIR = os.path.dirname(os.path.realpath(__file__))

  ################
 #  Antiproton  #
################

data_cirelli_pbar          = np.genfromtxt( CDIR+'/AtProduction_antiprotons_with_ZZstar_WWstar.dat',      skip_header=1 )

def getSpectrum_Cirelli(mDM, ff = 13):
    
    for i in range(1000):
        if data_cirelli_pbar[i*179,0]>mDM:
            break
    i-=1
    mDM_l = data_cirelli_pbar[ i   *179,0]
    mDM_u = data_cirelli_pbar[(i+1)*179,0]
    
    
    x       = data_cirelli_pbar[ i   *179:(i+1)*179, 1]
    dNdx_l  = data_cirelli_pbar[ i   *179:(i+1)*179, ff] + np.ones(179)*1e-90
    dNdx_u  = data_cirelli_pbar[(i+1)*179:(i+2)*179, ff] + np.ones(179)*1e-90
    
    log_dNdx_m  = np.log(dNdx_l)    +    (  np.log(dNdx_u)-np.log(dNdx_l) ) / (np.log(mDM_u)-np.log(mDM_l)) * (np.log(mDM)-np.log(mDM_l))
    
    return np.power(10,x)*mDM, np.exp(log_dNdx_m)/(  np.power(10,x)*mDM  )/np.log(10)
