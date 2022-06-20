import dendro
from sympy import *

###############################################################
#  initialize
###############################################################

PI = symbols('PI')

# declare functions
rho_e = dendro.scalar("rho_e","[pp]")
Gamma = dendro.scalar("Gamma","[pp]")

E = dendro.vec3("E", "[pp]") 
A = dendro.vec3("A", "[pp]") 

d = dendro.set_first_derivative('grad')    # first argument is direction


###############################################################
#  constraint equation 
###############################################################

divA = sum( d(i,A[i]) for i in dendro.e_i ) - Gamma 

divE = sum( d(i,E[i]) for i in dendro.e_i ) - 4.0 * PI * rho_e 

###############################################################
#  evolution equations
###############################################################

outs = [divA, divE]
vnames = ['divA', 'divE']
dendro.generate(outs, vnames, '[pp]')
