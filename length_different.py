#Global Imports
from numpy.linalg import *
from numpy import *
from scipy.special import *



# length_different takes in two vectors a and b and returns the length of a-b

def length_different(a,b):
    
    return norm(a-b)