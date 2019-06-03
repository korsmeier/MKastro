import numpy as np
import os

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



