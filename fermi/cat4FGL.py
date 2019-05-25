import numpy as np
import os

CDIR = os.path.dirname(os.path.realpath(__file__))

from astropy.io import fits
from astropy.table import Table

hdul        = fits.open(CDIR+'/data/4FGL.fit')
dictionary  = hdul['LAT_Point_Source_Catalog'].columns.names
format      = hdul['LAT_Point_Source_Catalog'].columns.formats
unit        = hdul['LAT_Point_Source_Catalog'].columns.units
data        = hdul['LAT_Point_Source_Catalog'].data
table       = Table(data)

energy_bounds = hdul['EnergyBounds'].data
_energy_bound_lower = energy_bounds['LowerEnergy']
_energy_bound_upper = energy_bounds['UpperEnergy']
energy_bound_lower=[]
energy_bound_upper=[]
for e in _energy_bound_lower:
    if not e in energy_bound_lower:
        energy_bound_lower.append(e)
for e in _energy_bound_upper:
    if not e in energy_bound_upper:
        energy_bound_upper.append(e)
energy_bounds = (np.array( [ energy_bound_lower, energy_bound_upper ] )).transpose()


def get( observable, row=-1  ):
    if not observable in dictionary:
        print('4FGL does not contain %s.' % observable )
        return None
    if row<0:
        return np.copy(data[observable][:])
    return data[observable][row]

def get_energy_bounds():
    return np.copy(energy_bounds)

def print_columns():
    print('Columns in 4FGL')
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

def n():
    return len(_ra)
