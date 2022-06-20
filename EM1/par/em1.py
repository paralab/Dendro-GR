import dendro
from sympy import *

###############################################################
#  initialize
###############################################################

PI = symbols('PI')

# declare functions
psi   = dendro.scalar("psi",  "[pp]")

A = dendro.vec3("A", "[pp]") 
E = dendro.vec3("E", "[pp]") 
J = dendro.vec3("J", "[pp]") 

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
d2 = dendro.d2


###############################################################
#  evolution equations
###############################################################

A_rhs = [ - E[i] - d(i,psi) for i in dendro.e_i ]

E_rhs = [ - sum( d2(j,j,A[i]) for j in dendro.e_i ) + sum( d2(i,j,A[j]) for j in dendro.e_i ) - 4.0 * PI * J[i] for i in dendro.e_i ] 

psi_rhs = - sum( d(i,A[i]) for i in dendro.e_i) 

###############################################################
#  evolution equations
###############################################################

outs = [A_rhs, E_rhs, psi_rhs]
vnames = ['A_rhs', 'E_rhs', 'psi_rhs']
dendro.generate(outs, vnames, '[pp]')
