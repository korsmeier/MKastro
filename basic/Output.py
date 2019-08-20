
__author__='M. Korsmeier - (korsmeier@physik.rwth-aachen.de)'
__version__='1.0'


  ###############
 #  Output     #
###############

# progress bar

'''
    Example:
    
    In python2:
    *  print( '\r'+out.bar(100.*i/n), end='' )
    *  sys.stdout.flush()
    
    '''
def bar(percent, length=50):
    """ Function returns a status bar.
        """
    hash =  int(percent*length/100.+0.01)
    empty = length-hash
    ret = '   |'
    for i in range(hash):
        ret += '#'
    for i in range(empty):
        ret += ' '
    return ret+('| %4.1f %%                   ' % percent)

