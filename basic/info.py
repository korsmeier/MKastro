import os
import subprocess
    
CDIR = os.path.dirname(os.path.realpath(__file__))

def write_git_SHA( fname ):
    os.system( 'echo "`date`" >> %s' % (fname) )
    if os.path.isfile(CDIR+'/../Git-SHA.txt'):
        os.system( 'echo "    `cat %s/../Git-SHA.txt`" >> %s' % (CDIR,fname) )
    else:
        os.system( 'echo "    `tail -n 1 %s/../.git/logs/HEAD | awk \'{print $2}\'`" >> %s' % (CDIR,fname) )

