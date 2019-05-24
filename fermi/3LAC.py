#!/usr/bin/env python

import numpy as np

data_3LAC_num    = np.genfromtxt( '3LAC.txt', delimiter=',', skip_header=1 )
col = len(data_3LAC_num[0,:])
row = len(data_3LAC_num[:,0])
data_3LAC_string = [x[:] for x in [[''] * col] * row]

data_3LAC_string = np.array(data_3LAC_string, dtype=object)

dictionary = []
f = open('3LAC.txt')
for i, l in enumerate(f):
    if i==0:
        dictionary = (l.split('\n')[0]).split(',')
    for j, it in enumerate(  (l.split('\n')[0]).split(',')  ):
        data_3LAC_string[i-1,j] = it
f.close()

print(dictionary)

string_columns = ['Fermi name', 'Counterpart name', 'Bzcat5 name', 'Optical Class', 'Spectral Type', 'SED Class', 'Radio flag', 'CLEAN']

def get( observable, row=-1  ):
    col = dictionary.index(observable)
    is_string = False
    if observable in string_columns:
        is_string = True
    if col<0:
        print('3LAC does not contain %s.' % observable )
        return None
    if row<0:
        if is_string:
            return data_3LAC_string[:, col]
        else:
            return np.copy(data_3LAC_num[:, col])
    if is_string:
        return data_3LAC_string[row, col]
    else:
        return np.copy(data_3LAC_num[row, col])

