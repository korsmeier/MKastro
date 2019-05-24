#!/usr/bin/env python

import numpy as np



def get_strings_from_column( name, col=0 ):
    vec = []
    f = open(name)
    for i, l in enumerate(f):
        if l[0]=='#':
            continue
        str_list = l.split(' ')
        while '' in str_list:
            str_list.remove('')
        vec.append( str_list[col].lower()  )
    f.close()
    return vec


_4FGL_pos  = np.genfromtxt( '4FGL_positions.dat.txt' )[:,:2]
_4FGL_name = get_strings_from_column( '4FGL_positions.dat.txt', 2 )

#print( _4FGL_pos [:10] )
#print( _4FGL_name[:10] )

#import sys
#sys.exit(0)

_4FGL_silvia_name = get_strings_from_column( '4FGL_an_bcut30deg.dat.txt', 0 )
#print( _4FGL_silvia_name[:10] )


_3LAC      = np.genfromtxt( '3LAC.txt', delimiter=',' )

_3LAC_pos  = _3LAC[:,3:5]
_3LAC_z    = _3LAC[:, 5]
#print( _3LAC_pos[:10] )
#print( _3LAC_z[:10] )


def get_4FGL_pos( __4FGL_name ):
    for i,n in enumerate(_4FGL_name):
        if n==__4FGL_name:
            return _4FGL_pos[i,:]
    return (None,None)

pm = 0
def get_4FGL_name( position ):
    global pm
    name='empty'
    for i,p in enumerate(_4FGL_pos):
        d2 = ( np.power( position-p , 2 ) ).sum()
        if d2<0.2**2:
            pm+=1
            name=_4FGL_name[i]
            break
    return name

nz = 0
def get_z_from_pos( position ):
    global nz
    z=-1
    for i,p in enumerate(_3LAC_pos):
        d2 = ( np.power( position-p , 2 ) ).sum()
        if d2<0.2**2:
            z=_3LAC_z[i]
            nz += 1
            break
    return z



#print( get_4FGL_name( np.array([ 0.3151 , -7.7971 ]) ) )
#for i,p in enumerate( _3LAC_pos ):
#    print( get_4FGL_name(p) )
#
#print( pm )
#import sys
#sys.exit(0)




#for i,n in enumerate(_4FGL_silvia_name):
#    n = n[5:]
#    print(n)
#    ps = get_4FGL_pos(n)
#    print( ps )
#    if ps[0] != None:
#        z = get_z_from_pos( ps )
##    else:
##        print('>%s<' % n)
#        print(z)
#    print('-----')
##    if i>20:
##        break
#
#print('Identified number of redshifts:')
#print(nz)


s = ''
f = open('4FGL_an_bcut30deg.dat.txt')
for j, l in enumerate(f):
#    if j>3:
#        break
#    print(l)
    if j==0:
        s += l
        continue
    i = j-1
    n = _4FGL_silvia_name[i]
    n = n[5:]
    ps = get_4FGL_pos(n)
    z = -1
    if ps[0] != None:
        z = get_z_from_pos( ps )

    l = list(l)
    l[25:33] = list( '%2.5f' % z )[:]
    l = ''.join(l)
    s += l
#    print(l)


f.close()

print('Identified number of redshifts:')
print(nz)

f = open('4FGL_an_bcut30deg.dat.txt', 'w')
f.write(s)
f.close()

