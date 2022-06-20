import dendro
from sympy import *

###############################################################
#  initialize
###############################################################

from sympy import LeviCivita

PI = symbols('PI')
kappa_1 = symbols('kappa_1') 
kappa_2 = symbols('kappa_2') 

# declare functions
rho_e = dendro.scalar("rho_e","[pp]") 
Phi = dendro.scalar("Phi","[pp]") 
Psi = dendro.scalar("Psi","[pp]") 

E = dendro.vec3("E", "[pp]") 
B = dendro.vec3("B", "[pp]") 
J = dendro.vec3("J", "[pp]") 

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
d2 = dendro.d2


###############################################################
#  evolution equations
###############################################################

B_rhs = [ - sum(sum(LeviCivita(i,j,k)*d(j,E[k]) for k in dendro.e_i) for j in dendro.e_i) + d(i,Phi)  for i in dendro.e_i ] 

E_rhs = [  sum(sum(LeviCivita(i,j,k)*d(j,B[k]) for k in dendro.e_i) for j in dendro.e_i) - 4.0 * PI * J[i] - d(i,Psi)  for i in dendro.e_i ] 

Phi_rhs = sum(d(i,B[i]) for i in dendro.e_i) - kappa_2 * Phi 

Psi_rhs = 4 * PI * rho_e - sum(d(i,E[i]) for i in dendro.e_i) - kappa_1 * Psi 

###############################################################
#  evolution equations
###############################################################

outs = [B_rhs, E_rhs, Phi_rhs, Psi_rhs]	 
vnames = ['B_rhs', 'E_rhs', 'Phi_rhs', 'Psi_rhs']
dendro.generate(outs, vnames, '[pp]')
