import dendro
from sympy import *

###############################################################
#  initialize
###############################################################

PI = symbols('PI')

# declare functions
rho_e = dendro.scalar("rho_e","[pp]")

E = dendro.vec3("E", "[pp]") 
B = dendro.vec3("B", "[pp]") 

d = dendro.set_first_derivative('grad')    # first argument is direction


###############################################################
#  constraint equation 
###############################################################

divE = sum( d(i,E[i]) for i in dendro.e_i ) - 4.0 * PI * rho_e 
divB = sum( d(i,B[i]) for i in dendro.e_i )			   

###############################################################
#  evolution equations
###############################################################

outs = [divE, divB]
vnames = ['divE', 'divB']
dendro.generate(outs, vnames, '[pp]')
