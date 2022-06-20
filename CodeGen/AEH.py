
## Symbolic code for computing the Apperent Event Horizon. 

import sys as sys
import dendro
from sympy import *

a   = dendro.scalar("alpha", "[pp]")
chi = dendro.scalar("chi", "[pp]")
K   = dendro.scalar("K", "[pp]")
Gt  = dendro.vec3("Gt", "[pp]")
b   = dendro.vec3("beta", "[pp]")
B   = dendro.vec3("B", "[pp]")
s   = dendro.vec3("s", "[pp]")

gt  = dendro.sym_3x3("gt", "[pp]")
At  = dendro.sym_3x3("At", "[pp]")

Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')

d2 = dendro.d2


dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

KK  = dendro.sym_3x3("KK", "[pp]")
theta_rhs = - (sum([d(i,s[i]) for i in dendro.e_i]) + sum( [KK[i,j]*s[i]*s[j] for i,j in dendro.e_ij ] ) -K ) 


outs = [theta_rhs]
vnames = ['theta_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

