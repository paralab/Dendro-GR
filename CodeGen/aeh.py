
## Symbolic code for computing the Apparent Event Horizon (AEH) finder 
import sys as sys
import dendro
import sympy

gt  = dendro.sym_3x3("gt" , "[pp]")
At  = dendro.sym_3x3("At" , "[pp]")

K   = dendro.scalar("K"   , "[pp]")
chi = dendro.scalar("chi" , "[pp]")

d   = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad  = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')
d2  = dendro.d2

dendro.set_metric(gt)
# inverse of the spatial metric g_{ij} = (1/chi) * gt_{ij}
ig = chi * dendro.get_inverse_metric()
#igt = dendro.sym_3x3("igt", "")
#dendro.inv_metric = igt
C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()

# full Christoffel symbols w.r.t g_{ij}
C3 = dendro.get_complete_christoffel(chi)

# r,theta, phi = sympy.symbols("r, theta, phi")

# def sph_real(l,m, theta, phi):
#     """
#     computes the real spherical harmonics
#     """
#     assert abs(m)<=l

#     if m==0:
#         sph_norm_fac = sympy.sqrt(sympy.Rational((2 * l + 1), 4))/sympy.sqrt(sympy.pi) 
#     else:
#         sph_norm_fac = (-1)**m * sympy.sqrt(sympy.Rational((2 * l + 1) * sympy.factorial(l-sympy.Abs(m)) , 4 * sympy.factorial(l+sympy.Abs(m)))) /sympy.sqrt(sympy.pi)

#     if m < 0:
#         return sph_norm_fac * sympy.functions.special.polynomials.assoc_legendre(l, m, sympy.cos(theta)) * sympy.sin(sympy.Abs(m) * phi)
#     elif m==0:
#         return sph_norm_fac * sympy.functions.special.polynomials.assoc_legendre(l, m, sympy.cos(theta)) 
#     else:
#         return sph_norm_fac * sympy.functions.special.polynomials.assoc_legendre(l, m, sympy.cos(theta)) * sympy.cos( m * phi)

# l_max = 5
# lm_modes = [(l, m) for l in range(0, l_max + 1) for m in range(-l, -l+1)]

F      = dendro.scalar("F","[pp]")
s_dk   = [d(i,F) for i in dendro.e_i]
s_uk   = [sympy.simplify(sum([ig[i,j] * s_dk[j] for j in dendro.e_i])) for i in dendro.e_i]
s_norm = sympy.sqrt(sympy.simplify(sum([s_dk[i] * s_uk[i] for i in dendro.e_i])))
n_uk   = [sympy.simplify(s_uk[i]/s_norm) for i in dendro.e_i]

Kij    = (1/chi) * (At + sympy.Rational(1,3) * gt * K)

H      = sum([(dendro.DiDj(F)[i,j]/s_norm - Kij[i,j]) * (ig[i,j] - n_uk[i] * n_uk[j]) for i in dendro.e_i for j in dendro.e_i])

outs   = [H]
vnames = ['H']
dendro.generate_cpu(outs, vnames, '[pp]')


#mu_ij = sympy.Matrix([[sympy.simplify(ig[i,j] - n_uk[i] * n_uk[j]) for j in dendro.e_i] for i in dendro.e_i ])
