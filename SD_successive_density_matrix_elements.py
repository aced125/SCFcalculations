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

from HFvariables import *



# Takes in two matrices and returns their standard deviation

def SD_successive_density_matrix_elements(Ptilde, P):
    
    x = 0
    
    for i in range(basis_set_size):
        for j in range(basis_set_size):
            x += (basis_set_size**-2)*(Ptilde[i,j]-P[i,j])**2
            
    return x**0.5