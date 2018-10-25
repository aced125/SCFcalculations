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






# multi takes in two tuples of information: (exponent, orbital centre)
# and returns the tensor (ab|cd)
# need to give centre as python list
# also need to feed atom_coordinate as python list, starting from 1

def multi(tup1, tup2, tup3, tup4):
    
    a, Ra = primgauss(tup1[0],tup1[1])
    
    b, Rb = primgauss(tup2[0],tup2[1])
    
    c, Rc = primgauss(tup3[0],tup3[1])
    
    d, Rd = primgauss(tup4[0],tup4[1])
    
    p = a+b
    
    q = c+d
    
    Rp = (1/(a+b))*(a*Ra+b*Rb)
    
    Rq = (1/(c+d))*(c*Rc+d*Rd)
    
    
    return (2*pi**2.5*((a+b)*(c+d)*(p+q)**0.5)**-1)*overlap(tup1,tup2)*overlap(tup3,tup4)*(p/pi)**1.5*(q/pi)**1.5*Fo(square(length_different(Rp,Rq))*(a+b)*(c+d)/(p+q))