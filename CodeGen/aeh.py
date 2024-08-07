## Symbolic code for computing the Apparent Event Horizon (AEH) finder 
import sys as sys
import dendro
import sympy
import contextlib 
import io

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

F      = dendro.scalar("F","[pp]")
s_dk   = [d(i,F) for i in dendro.e_i]
s_uk   = [sympy.simplify(sum([ig[i,j] * s_dk[j] for j in dendro.e_i])) for i in dendro.e_i]
s_norm = sympy.sqrt(sympy.simplify(sum([s_dk[i] * s_uk[i] for i in dendro.e_i])))
n_uk   = [sympy.simplify(s_uk[i]/s_norm) for i in dendro.e_i]

Kij    = (1/chi) * (At + sympy.Rational(1,3) * gt * K)

H      = sum([(dendro.DiDj(F)[i,j]/s_norm - Kij[i,j]) * (ig[i,j] - n_uk[i] * n_uk[j]) for i in dendro.e_i for j in dendro.e_i])

outs   = [H]
vnames = ['H']

with contextlib.redirect_stdout(io.StringIO()) as f:
    dendro.generate_cpu(outs, vnames, '[pp]')

with open("../BSSN_GR/src/expansion_aeh.cpp", "w") as f_out:
    f_out.writelines(f.getvalue())

dendro.d   = lambda i,x : sympy.Symbol("grad_%d_%s"%(i,str(x).split('[')[0]))
dendro.d2  = lambda i,j,x : sympy.Symbol("grad2_%d_%d_%s"%(min(i,j),max(i,j),str(x).split('[')[0]))

dendro.ad  = dendro.d
dendro.kod = dendro.undef

d=dendro.d
d2=dendro.d2
ad=dendro.ad
    
g      = gt/chi
A      = dendro.vec3("A","[pp]")


# Define s
s = sympy.symbols('xx yy zz')

# Define the Killing vectors for rotational symmetry
Killing_x = [0, -(s[2]), (s[1])]
Killing_y = [(s[2]), 0, -(s[0])]
Killing_z = [-(s[1]), (s[0]), 0]

# Define the integrand for each angular momentum component
def integrand(phi, Kij):
    # phi - is a Cartesian killing vector
    # Kij - is the ADM extrinsic curvature
    # n - is the unit normal to the surface
    
    # redefine without [pp] index, [pp] index is needed for the other parts of the computation. 
    F      = dendro.scalar("F","")
    s_dk   = [d(i,F) for i in dendro.e_i]
    s_uk   = [sympy.simplify(sum([ig[i,j] * s_dk[j] for j in dendro.e_i])) for i in dendro.e_i]
    s_norm = sympy.sqrt(sympy.simplify(sum([s_dk[i] * s_uk[i] for i in dendro.e_i])))
    n_uk   = [sympy.simplify(s_uk[i]/s_norm) for i in dendro.e_i]
    
    return sum(sum( phi[i] * Kij[i, j] * n_uk[j] for i in range(3)) for j in range(3))

md_ij        = sympy.Matrix([[sympy.simplify(sum([g[a,b] * d(i, A[a]) * d(j, A[b]) for a in dendro.e_i for b in dendro.e_i])) for j in range(2)] for i in range(2)])
det_mij      = sympy.simplify(sympy.det(md_ij))
sqrt_det_mij = sympy.sqrt(sympy.Abs(det_mij))

# Define integral for angular momentum components
J_x = sqrt_det_mij * sympy.simplify(integrand(Killing_x, Kij))
J_y = sqrt_det_mij * sympy.simplify(integrand(Killing_y, Kij))
J_z = sqrt_det_mij * sympy.simplify(integrand(Killing_z, Kij))

# Add the integrals to the list of outputs
outs = [sqrt_det_mij, J_x, J_y, J_z]
vnames = ['sqrt_det_m_ab', 'J_x', 'J_y', 'J_z']

with contextlib.redirect_stdout(io.StringIO()) as f:
    dendro.generate_cpu(outs, vnames, '[pp]')

with open("../BSSN_GR/src/det_metric_aeh.cpp", "w") as f_out:
    f_out.writelines(f.getvalue())
