import numpy as np
import os
import sys

CDIR = os.path.dirname(os.path.realpath(__file__))

from astropy.io import fits
from astropy.table import Table

hdul        = fits.open(CDIR+'/data/4LAC.fit')


dictionary  = hdul['4LAC AGNs'].columns.names
format      = hdul['4LAC AGNs'].columns.formats
unit        = hdul['4LAC AGNs'].columns.units
data        = hdul['4LAC AGNs'].data
table       = Table(data)

def get( observable, row=-1  ):
    if not observable in dictionary:
        print('4LAC does not contain %s.' % observable )
        return None
    if row<0:
        return np.copy(data[observable][:])
    return data[observable][row]

def get_energy_bounds():
    return np.copy(energy_bounds)

def print_columns():
    print('Columns in 4LAC')
    for i, d in enumerate( dictionary ):
        print(' %-30s format: %-10s unit: %-10s' % (d, format[i], unit[i]) )
def columns():
    return dictionary

_ra  = data['RAJ2000']
_dec = data['DEJ2000']
def _find_source( ra, dec, min_separation=0.1 ):
    global _ra, _dec
    d2 = (ra - _ra)**2 + (dec - _dec)**2
    i = np.argmin(d2)
    if np.sqrt(d2[i])<min_separation:
        return i
    return -1
find_source = np.vectorize(_find_source)


_4LAC_names = data['Source_Name']
def _find_source_from4FGL( _4FGL_name ):
    global _4LAC_names
    ia = np.where(_4LAC_names==_4FGL_name)[0]
    if len(ia):
        return ia[0]
    return -1
find_source_from4FGL = np.vectorize(_find_source_from4FGL)

def n():
    return len(_ra)
