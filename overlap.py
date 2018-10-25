#Global Imports

import numpy as np
from numpy import pi
from numpy import exp
from numpy import linalg
from numpy.linalg import *
from numpy import *

from scipy.special import *



from primgauss import *
from length_different import *

# overlap takes in two tuples of with information: (exponent, orbital centre)
# and returns overlap
# need to give centre as a python list



def overlap(tup1,tup2):
    
    a, Ra = primgauss(tup1[0],tup1[1])
    
    b, Rb = primgauss(tup2[0],tup2[1])
    
    p = a+b
    
    
    return ((((2*a/pi)**0.75)*((2*b/pi)**0.75)*(((pi/p)))**1.5)*exp(((-a*b)/p)*square(length_different(Ra,Rb))))
