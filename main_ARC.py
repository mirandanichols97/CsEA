

import sys, os
rootDir = '/var/www/html/'
sys.path.append(rootDir)
os.chdir(rootDir)

from arc import *
from arc.web_functionality import *

def transitTime (beam_waist,mean_speed):
    # in s
    beam_fwhm = beam_waist / 0.8493218
    return sqrt(pi)*beam_fwhm/(2.*mean_speed)


# = = = = = = = = = WRITE YOUR CODE BELOW = = = = = = = = =

if __name__ == '__main__':
    atom = Rubidium()
    
