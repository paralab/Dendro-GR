#!/usr/bin/env python3
##Author: Milinda Fernando
##School of Computing, University of Utah 
##03/1/2018
##Computes finite difference coefficients for given sample points given derivative.  
##

from __future__ import print_function
from sympy import *
x, x0, h = symbols('x, x_i, h')
# specify the order of the fd approximation
p=2
sample_loc=[-1,0,1]
F=symbols('F:3')
c=symbols('c:3')
# define a polynomial of degree n
def P(h,n):
	return sum( ((1/factorial(i))*(h)**i for i in range(n)) )

# now we make a matrix consisting of the coefficients
# of the c_i in the nth degree polynomial P
# coefficients of c_i evaluated at x_i -h 
m11 = P(x0-h , x0, c, n).diff(c[0])
m12 = P(x0-h , x0, c, n).diff(c[1])
m13 = P(x0-h , x0, c, n).diff(c[2])
# coefficients of c_i evaluated at x_i 
m21 = P(x0, x0, c, n).diff(c[0])
m22 = P(x0, x0, c, n).diff(c[1])
m23 = P(x0, x0, c, n).diff(c[2])
# coefficients of c_i evaluated at x_i + h
m31 = P(x0+h, x0, c, n).diff(c[0])
m32 = P(x0+h, x0, c, n).diff(c[1])
m33 = P(x0+h, x0, c, n).diff(c[2])
# matrix of the coeffcients is 3x3 in this case
M = Matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
print(M)
#Now that we have the matrix of coefficients we next form the right-hand-side and solve by inverting `M`:
# matrix of the function values...actually a vector of right hand sides
R = Matrix([[Fi], [Fim1], [Fip1]])
# matrix form of the three equations for the c_i is M*X = R
# solution directly inverting the 3x3 matrix M:
X =  M.inv() * R
# note that all three coefficients make up the solution
# the first derivative is coefficient c_1 which is X[1].
print("The second-order accurate approximation for the first derivative is: ")
print(together(X[1]))

