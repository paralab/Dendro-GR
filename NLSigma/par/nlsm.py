import dendro
from sympy import *

###############################################################
#  initialize
###############################################################

r = symbols('r')

# declare functions
chi = dendro.scalar("chi","[pp]")
phi = dendro.scalar("phi","[pp]")

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
d2 = dendro.d2



###############################################################
#  evolution equations
###############################################################

phi_rhs = sum( d2(i,i,chi)  for i in dendro.e_i ) - sin(2*chi)/r**2

chi_rhs = phi

###############################################################
#  evolution equations
###############################################################

outs = [phi_rhs, chi_rhs]
vnames = ['phi_rhs', 'chi_rhs']

dendro.generate(outs, vnames, '[pp]')
