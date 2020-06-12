import numpy as np
import os

CDIR = os.path.dirname(os.path.realpath(__file__))

from astropy.io import fits
from astropy.table import Table

hdul        = fits.open(CDIR+'/data/3FGL.fit')
#print(hdul)
#print(hdul.info())
dictionary  = hdul['LAT_Point_Source_Catalog'].columns.names
#print(dictionary)
format      = hdul['LAT_Point_Source_Catalog'].columns.formats

unit        = hdul['LAT_Point_Source_Catalog'].columns.units
data        = hdul['LAT_Point_Source_Catalog'].data
table       = Table(data)


def get( observable, row=-1  ):
    if not observable in dictionary:
        print('3FGL does not contain %s.' % observable )
        return None
    if row<0:
        return np.copy(data[observable][:])
    return data[observable][row]

def print_columns():
    print('Columns in 4FGL')
    for i, d in enumerate( dictionary ):
        print(' %-30s format: %-10s unit: %-10s' % (d, format[i], unit[i]) )

#print_columns()

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

def n():
    return len(_ra)
