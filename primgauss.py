# primgauss takes in the two parameters of the primitive gaussian, its exponent
# and which atom it is centred on, and returns a tuple of the exponent and 
# vector of the atom it is centred on

from atomvector import *


def primgauss(exponent, atom):
    
    return exponent, atomvector(atom)
