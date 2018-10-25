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


from overlap import *









# kinetic takes in two tuples of information: (exponent, orbital centre)
# and returns kinetic energy element
# need to give centre as python list

def kinetic(tup1, tup2):
    
    a, Ra = primgauss(tup1[0],tup1[1])
    
    b, Rb = primgauss(tup2[0],tup2[1])
    
    p = a+b
    
    return (a*b/p)*(3+2*((-a*b)/p)*square(length_different(Ra,Rb)))*overlap(tup1,tup2)