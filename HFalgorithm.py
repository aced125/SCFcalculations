#Global Imports
from numpy import *
from scipy.special import *



#Personal Imports
from HFvariables import *
from SD_successive_density_matrix_elements import *



#Closed form expression for integrals
 
from overlap import *
from kinetic import *
from potential import *
from multi import *

import contraction_matrix
from contraction_matrix import *



# First evaluate all the required integrals that don't change in the iterative process

#Overlap matrix in the contracted gaussian basis

S = np.zeros((basis_set_size, basis_set_size))

for mew in range(basis_set_size):
    for v in range(basis_set_size):
        for p in range(STOLG):
            for q in range(STOLG):
                
                S[mew, v] += d[p]*d[q]*overlap((alpha[p,mew],atomcentres[mew]),(alpha[q,v],atomcentres[v]))



#Kinetic energy matrix in the contracted gaussian basis

T = np.zeros((basis_set_size, basis_set_size))

for mew in range(basis_set_size):
    for v in range(basis_set_size):
        for p in range(STOLG):
            for q in range(STOLG):
                
                T[mew, v] += d[p]*d[q]*kinetic((alpha[p,mew],atomcentres[mew]),(alpha[q,v],atomcentres[v]))               




#Potential energy matrix in the contracted gaussian basis               
                
V = np.zeros((basis_set_size, basis_set_size))

for i in range(1,Natoms+1):
    for mew in range(basis_set_size):
        for v in range(basis_set_size):
            for p in range(STOLG):
                for q in range(STOLG):
                    V[mew, v] += d[p]*d[q]*potential((alpha[p,mew],atomcentres[mew]),(alpha[q,v],atomcentres[v]),i)







# The multi electron tensor Mijkl = (ij|kl) where chemists' notation has been used

multi_electron_tensor = np.zeros((basis_set_size,basis_set_size,basis_set_size,basis_set_size))

for mew in range(basis_set_size):
    for v in range(basis_set_size):
        for x in range(basis_set_size):
            for y in range(basis_set_size):
                for p in range(3):
                    for q in range(3):
                        for pprime in range(3):
                            for qprime in range(3):
                                multi_electron_tensor[mew, v, x, y] += d[p]*d[q]*d[pprime]*d[qprime]*multi((alpha[p,mew],atomcentres[mew]),(alpha[pprime,v],atomcentres[v]),(alpha[q,x],atomcentres[x]),(alpha[qprime,y],atomcentres[y]))




# Symmetric orthogonalisation of basis pp144 Atilla and Szabo

evalS, U = eig(S)

diagS = dot(U.T,dot(S,U))

diagS_minushalf = diag(diagonal(diagS)**-0.5)

X = dot(U,dot(diagS_minushalf,U.T))





Hcore = T + V



#Initial guess at P

P = np.zeros((basis_set_size,basis_set_size))
Pprevious = np.zeros((basis_set_size,basis_set_size))

Plist = []


#algorithm

Decider = 5

while Decider > 10**-4:
    # Calculate Fock matrix with this guesss
    Fock = Hcore

    Hcore = T + V

    for i in range(basis_set_size):
        for j in range(basis_set_size):
            for x in range(basis_set_size):
                for y in range(basis_set_size):
                    Fock[i,j] += P[x,y]*(multi_electron_tensor[i,j,y,x]-0.5*multi_electron_tensor[i,x,y,j])


    # Calculate Fock matrix in orthogonalised base

    Fockprime = dot(X.T,dot(Fock, X))

    evalFockprime, Cprime = eig(Fockprime)


    #Correct ordering of eigenvalues and eigenvectors (starting from ground MO as first column of C, else we get the wrong P)

    idx = evalFockprime.argsort()
    evalFockprime = evalFockprime[idx]
    Cprime = Cprime[:,idx]

    C = dot(X,Cprime) 

    #Form new P

    for i in range(basis_set_size):
        for j in range(basis_set_size):
            for a in range(int(n/2)):
                P[i,j] = 2*C[i,a]*C[j,a]
                
                
    Plist.append(P)

    Decider = SD_successive_density_matrix_elements(Pprevious,P)
    Pprevious = P.copy()
    
print('\n')
print('STO3G Restricted Closed Shell HF algorithm took {} iterations to converge'.format(len(Plist)))
print('\n')
print('The orbital energies are {} and {} Hartrees'.format(evalFockprime[0],evalFockprime[1]))
