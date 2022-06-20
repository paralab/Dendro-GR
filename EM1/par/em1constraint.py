import dendro
from sympy import *

###############################################################
#  initialize
###############################################################

PI = symbols('PI')

# declare functions
rho_e = dendro.scalar("rho_e","[pp]")

E = dendro.vec3("E", "[pp]") 

d = dendro.set_first_derivative('grad')    # first argument is direction


###############################################################
#  constraint equation 
###############################################################

divE = sum( d(i,E[i]) for i in dendro.e_i ) - 4.0 * PI * rho_e 

###############################################################
#  evolution equations
###############################################################

outs = [divE]
vnames = ['divE']
dendro.generate(outs, vnames, '[pp]')
