####################################################################
# May.2020
# Quad grav rhs generator new version
#####################################################################

import dendro
from sympy import *

import numpy as np
###################################################################
# initialize
###################################################################

# declare variables (BSSN vars)
a   = dendro.scalar("alpha", "[pp]")
chi = dendro.scalar("chi", "[pp]")
K   = dendro.scalar("K", "[pp]")

Gt  = dendro.vec3("Gt", "[pp]")
b   = dendro.vec3("beta", "[pp]")
B   = dendro.vec3("B", "[pp]")

gt  = dendro.sym_3x3("gt", "[pp]")
At  = dendro.sym_3x3("At", "[pp]")

Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")

Atr = dendro.scalar("Atr", "[pp]")
Aij  = dendro.sym_3x3("Aij", "[pp]")
# From V_ab 
Btr = dendro.scalar("Btr", "[pp]")
Bij  = dendro.sym_3x3("Bij", "[pp]")


# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

# specify the functions for computing first and second derivatives
# d = dendro.set_first_derivative('grad')    # first argument is direction
# d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
# ad = dendro.set_advective_derivative('agrad')  # first argument is direction
# kod = dendro.set_kreiss_oliger_dissipation('kograd')
d   = lambda i,x : symbols("grad_%d_%s"%(i,x))
ad  = lambda i,x : symbols("agrad_%d_%s"%(i,x))
kod = lambda i,x : symbols("kograd_%d_%s"%(i,x))
d2  = lambda i,j,x : symbols("grad2_%d_%d_%s"%(min(i,j),max(i,j),x))
dendro.d    = d
dendro.ad   = ad
dendro.kod  = kod
dendro.d2   = d2

#f = Function('f')

# generate metric related quantities
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)
###################################################################
# evolution equations
###################################################################
# normal vector n^a
n_vec = Matrix([[1/a, -b[0]/a, -b[1]/a, -b[2]/a]])
# Precomputataion of covariant derviative
Djni= Matrix([d(n_vec[i],j) + sum([dendro.C3[k,j,i]*n_vec[k] for k in dendro.e_i]) for i,j in dendro.e_ij]).reshape(3,3)
Aij_rhs = Djni

###################################################################
# generate code
###################################################################

outs = [Aij_rhs]
vnames = ['Aij_rhs']
dendro.generate(outs, vnames, '[pp]')
