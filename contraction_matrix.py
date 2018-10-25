# matrix alpha's ijth element corresponds to the ith sum in the contraction
# and the jth atom

#Global Imports

import numpy as np
from numpy import pi
from numpy import exp
from numpy import linalg
from numpy.linalg import *
from numpy import *

from scipy.special import *




import HFvariables
from HFvariables import *


alphavec = [0.109818, 0.405771, 2.22766]

alpha = []

for i in range(Natoms):
    alpha += [x*zeta[i]**2 for x in alphavec]
    

    
alpha = np.asarray(alpha).reshape(2,3).T