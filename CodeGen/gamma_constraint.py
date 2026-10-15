####################################################################
#
# Python script that generates the conformal connection (Gamma)
# constraint C^i = Gt^i - gt^{jk} Gt^i_{jk} and its magnitude
# C = sqrt(gt_{ij} C^i C^j). Output goes to BSSN_GR/src/physconeqs_gt.cpp.
#
####################################################################

#!/usr/bin/env/ python3

import dendro
from sympy import *

chi = dendro.scalar("chi", "[pp]")
Gt = dendro.vec3("Gt", "[pp]")
gt = dendro.sym_3x3("gt", "[pp]")

d = dendro.set_first_derivative("grad")
d2 = dendro.set_second_derivative("grad2")

dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)

cgt = [Gt[i] - CalGt[i] for i in dendro.e_i]
cgt_mag = sqrt(sum(gt[i, j] * cgt[i] * cgt[j] for i, j in dendro.e_ij))

outs = [cgt, cgt_mag]
vnames = ["cgt", "cgt_mag"]
dendro.generate_cpu(outs, vnames, "[pp]")
