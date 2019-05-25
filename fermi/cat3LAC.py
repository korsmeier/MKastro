import numpy as np
import os

CDIR = os.path.dirname(os.path.realpath(__file__))

data_3LAC_num    = np.genfromtxt( CDIR+'/data/3LAC.txt', delimiter=',', skip_header=1 )
col = len(data_3LAC_num[0,:])
row = len(data_3LAC_num[:,0])
data_3LAC_string = [x[:] for x in [[''] * col] * row]

data_3LAC_string = np.array(data_3LAC_string, dtype=object)

dictionary = []
f = open(CDIR+'/data/3LAC.txt')
for i, l in enumerate(f):
    if i==0:
        dictionary = (l.split('\n')[0]).split(',')
    for j, it in enumerate(  (l.split('\n')[0]).split(',')  ):
        data_3LAC_string[i-1,j] = it
f.close()

string_columns = ['Fermi name', 'Counterpart name', 'Bzcat5 name', 'Optical Class', 'Spectral Type', 'SED Class', 'Radio flag', 'CLEAN']

def n():
    return row

def get( observable, row=-1  ):
    if not observable in dictionary:
        print('3LAC does not contain %s.' % observable )
        return None
    col = dictionary.index(observable)
    is_string = False
    if observable in string_columns:
        is_string = True
    if row<0:
        if is_string:
            return data_3LAC_string[:, col]
        else:
            return np.copy(data_3LAC_num[:, col])
    if is_string:
        return data_3LAC_string[row, col]
    else:
        return np.copy(data_3LAC_num[row, col])

def print_columns():
    print('Columns in 3LAC')
    for d in dictionary:
        type = 'float'
        if d in string_columns: type = 'str'
        print(' %-30s%-10s' % (d, type) )
def columns():
    return dictionary

col_ra  = dictionary.index('RA (J2000.0)' )
col_dec = dictionary.index('Dec (J2000.0)')
def find_source( ra, dec, min_separation=0.1 ):
    global col_ra, col_dec
    _ra  = data_3LAC_num[:,col_ra ]
    _dec = data_3LAC_num[:,col_dec]
    d2 = (ra - _ra)**2 + (dec - _dec)**2
    i = np.argmin(d2)
    if np.sqrt(d2[i])<min_separation:
        return i
    return -1

