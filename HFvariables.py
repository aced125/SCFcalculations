

# Variables and constants

#R = Bond length in Bohr
#zeta = list of Slater coefficients for the different atoms (starting with A)
#Za = atom A atomic number (Helium)
#Zb = atom B atomic number (hydrogen)
#n = number of electrons
#STOLG = number of primitive gaussians used in the expansion of the contracted gaussian



IOP = 2
R = 1.4632
zeta = [2.0925, 1.2400]
Natoms = 2
Z = [2,1]
n = 2
STOLG = 3

d = [0.444635, 0.535328, 0.154329]

basis_set_size = 2


atomA = (zeta[0], [0,0,0])	
atomB = (zeta[0], [0,0,R])

atomcentres = [[0,0,0],[0,0,R]]