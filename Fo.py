from numpy import *

# Fo function as described in Szabo and Ostlund pp415

def Fo(t):
    
    if t==0:
        return 1
    
    else:
        return (0.5*(pi/t)**0.5)*erf(t**0.5)
    