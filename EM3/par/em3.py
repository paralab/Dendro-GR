import dendro
from sympy import *

###############################################################
#  initialize
###############################################################

from sympy import LeviCivita

PI = symbols('PI')

# declare functions


B = dendro.vec3("B", "[pp]") 
E = dendro.vec3("E", "[pp]") 
J = dendro.vec3("J", "[pp]") 

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
d2 = dendro.d2


###############################################################
#  evolution equations
###############################################################

#B_rhs = - [ if i = 0: d(i+1,B[i+2]) - d(i+2,B[i+1])
#		  elif = 1: d(i+1,B[i-1]) - d(i-1,B[i+1])
#			elif i = 2: d(i-2,B[i-1]) - d(i-1,B[i-2])		
#				for i in dendro.e_i ]

#B_rhs = - [ sum(sum(epsilon[i][j][k]*d(j,E[k]),for k in dendro.e_i)for j in dendro.e_i)] 
B_rhs = [ - sum(sum(LeviCivita(i,j,k)*d(j,E[k]) for k in dendro.e_i) for j in dendro.e_i) for i in dendro.e_i ] 
E_rhs = [  sum(sum(LeviCivita(i,j,k)*d(j,B[k]) for k in dendro.e_i) for j in dendro.e_i)
           - 4.0 * PI * J[i] for i in dendro.e_i ] 

#E_rhs =  [ sum(sum(epsilon[i][j][k]*d(j,B[k]),for k in dendro.e_i)for j in dendro.e_i) - 4.0 * PI * J[i] for i in dendro.e_i ] 


#E_rhs = [ if i = 0: d(i+1,E[i+2]) - d(i+2,E[i+1])
#		 elif i = 1: d(i+1,E[i-1]) - d(i-1,E[i+1])
#			elif i = 2: d(i-2,E[i-1]) - d(i-1,E[i-2])		
#				- 4.0 * PI * J[i] for i in dendro.e_i ] 


###############################################################
#  evolution equations
###############################################################

outs = [B_rhs, E_rhs]	 
vnames = ['B_rhs', 'E_rhs']
dendro.generate(outs, vnames, '[pp]')
