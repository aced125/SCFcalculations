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

from HFvariables import *
from Fo import *







# potential takes in two tuples of information: (exponent, orbital centre)
# and returns potential energy element of jth atom
# need to give centre as python list
# also need to feed atom_coordinate as a number starting from 1, not 0



def potential(tup1, tup2, atom_number):
    
    a, Ra = primgauss(tup1[0],tup1[1])
    
    b, Rb = primgauss(tup2[0],tup2[1])
    
    p = a+b
    
    Rp = (1/(a+b))*(a*Ra+b*Rb)
    
    Rc = np.asarray(atomcentres[atom_number-1])
    
    Zc = Z[atom_number-1]
    
    
    return -2*Zc*((p/pi)**0.5)*overlap(tup1,tup2)*Fo((a+b)*square(length_different(Rp,Rc)))