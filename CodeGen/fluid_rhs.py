'''
 The RHS of the fluid equations -- including geometric pieces. 
'''
import sys as sys
import dendro
from sympy import *
import math 

###################################################################
# initialize
###################################################################

# declare geometric variables
a    = dendro.scalar("alpha", "[pp]")
chi  = dendro.scalar("chi", "[pp]")
K    = dendro.scalar("K", "[pp]")
beta = dendro.vec3("beta", "[pp]")
gt   = dendro.sym_3x3("gt", "[pp]")
At   = dendro.sym_3x3("At", "[pp]")

#Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")

# declare hydrodynamic quantities (without B field)  
#vd = dendro.vec3 ("vd" , "[pp]")  
Sd = dendro.vec3 ("Sd" , "[pp]")  
D = dendro.scalar ("D", "[pp]") 
tau = dendro.scalar ("tau", "[pp]") 
Tdd_stress_dens = dendro.sym_3x3("Tdd_stress_dens", "[pp]")  

# declare RHSs
Sd_rhs = dendro.vec3 ("Sd_rhs" , "[pp]") 
D_rhs = dendro.scalar ("D_rhs" , "[pp]") 
tau_rhs = dendro.scalar ("tau_rhs" , "[pp]") 

# specify the functions for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction

#d2 = dendro.d2

dendro.set_metric(gt)
igt = dendro.get_inverse_metric()


#  This should be provided before coming into these RHS routines.  If not, it  
#  is defined here but note that we will need the velocity (vd[i]) somehow.  
#  Note, the two versions here are the same up to densitization.  The first  
#  is not densitized while the second is.  

#Tdd_stress_dens = Matrix([ 0.5*( vd[i]*Sd[j] + vd[j]*Sd[i] ) + (P/chi)*gt[i,j] for i,j in dendro.e_ij ])  

#Tdd_stress_dens = chi^(-1.5) * Matrix([ 0.5*( vd[i]*Sd[j] + vd[j]*Sd[i] ) + (P/chi)*gt[i,j] for i,j in dendro.e_ij ])  


#  As part of the RHS of the momentum equation, we need the usual Christoffels
#  built from the full (non-conformal) metric. We call it C3 in dendro.py:  
C3 = dendro.get_complete_christoffel(chi)


#  This version of rhs of the Sd (momentum) equation assumes that the 
#  spatial components of the stress tensor (\perp T or Tdd_stress_dens) 
#  as well the components of Sd, tau and D come in as *densitized* 
#  variables, i.e. carrying a factor of the determinant of the (full or 
#  non-conformal) 3-metric (usually referred to as gamma_ij or h_ij).  
Tud_stress_dens = chi * dendro.up_down(Tdd_stress_dens) 

Sd_rhs = Matrix([ a * sum([ sum([ C3[i,j,k] * Tud_stress_dens[j,i] for j in dendro.e_i ])  for i in dendro.e_i ]) + sum([ Sd[i] * d(k,beta[i]) for i in dendro.e_i ]) - (tau+D) * d(k,a) for k in dendro.e_i ]) 

Sd_rhs = [item for sublist in Sd_rhs.tolist() for item in sublist]


# The tau rhs.  But first the stress tensor with up indices.  
Tuu_stress_dens = chi * chi * dendro.up_up(Tdd_stress_dens) 

tau_rhs = (a/chi) * ( sum([ Tuu_stress_dens[i,i] * (At[i,i] + gt[i,i]*(K/3)) for i in dendro.e_i ]) + 2*sum([ (Tuu_stress_dens[i,j] * (At[i,j] + gt[i,j]*(K/3))) for i,j in dendro.e_ij_offdiag ]) ) - sum([ (Sd[i] * d(j,a) * igt[i,j]) for i,j in dendro.e_ij ])   


D_rhs = sympify(0.0) 

# the magnetic field if and when we use MHD 
#Bu_rhs = simplify(0.0) 
        

###################################################################
# generate code
###################################################################

#outs = [Sd_rhs, tau_rhs, D_rhs, Bu_rhs]
#vnames = ['Sd_rhs', 'tau_rhs', 'D_rhs', 'Bu_rhs']
outs = [Sd_rhs, tau_rhs, D_rhs]
vnames = ['Sd_rhs', 'tau_rhs', 'D_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

#numVars=len(outs)
#for i in range(0,numVars):
#    dendro.generate_separate([outs[i]],[vnames[i]],'[pp]')

