import dendro
import networkx as nx
import nxgraph 
from sympy import *
import sympy
import queue
from dataclasses import dataclass, field
from typing import Any
import contextlib 
import io

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

# Additional parameters for damping term
R0 = symbols('BSSN_ETA_R0')
ep1, ep2 = symbols('BSSN_ETA_POWER[0] BSSN_ETA_POWER[1]')

xi1, xi2, xi3 = symbols('BSSN_XI[0] BSSN_XI[1] BSSN_XI[2] ')
eta_damp=symbols("eta")

dendro.d   = lambda i,x : Symbol("grad_%d_%s"%(i,str(x).split('[')[0]))
dendro.d2  = lambda i,j,x : Symbol("grad2_%d_%d_%s"%(min(i,j),max(i,j),str(x).split('[')[0]))

dendro.ad  = dendro.d
dendro.kod = dendro.undef

d=dendro.d
d2=dendro.d2
ad=dendro.ad

# declare variables
a   = dendro.scalar("alpha", "[pp]")
chi = dendro.scalar("chi", "[pp]")
K   = dendro.scalar("K", "[pp]")

Gt  = dendro.vec3("Gt", "[pp]")
b   = dendro.vec3("beta", "[pp]")
B   = dendro.vec3("B", "[pp]")

gt  = dendro.sym_3x3("gt", "[pp]")
At  = dendro.sym_3x3("At", "[pp]")

Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

dendro.set_metric(gt)
igt = dendro.get_inverse_metric()
C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
C2_spatial = dendro.get_complete_christoffel(chi)
[R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

a_rhs = l1*dendro.lie(b, a) - 2*a*K 
        
b_rhs = [(Rational(3,4) * (lf0 + lf1*a) * B[i] + l2 * dendro.vec_j_ad_j(b, b[i])) for i in dendro.e_i ] 
    
gt_rhs = dendro.lie(b, gt, weight) - 2*a*At 

chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K)

AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) 

K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) 

At_UU = dendro.up_up(At)

Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
        Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
        Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
        Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
        Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
        Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
        Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])
    
Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]
_Gt_rhs = dendro.vec3("Gt_rhs","[pp]")
B_rhs = [ (_Gt_rhs[i] - eta_damp * B[i] +
        l3 * dendro.vec_j_ad_j(b, B[i]) -
        l4 * dendro.vec_j_ad_j(b, Gt[i]))
        for i in dendro.e_i ]

outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']


g=nxgraph.ExpressionGraph()
g.add_expressions(outs,vnames,'[pp]')
