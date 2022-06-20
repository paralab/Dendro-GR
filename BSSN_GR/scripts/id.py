'''

Mon Jan  1 11:50:27 MST 2018  
David Neilsen

This script estimates the black hole puncture parameters for initial data.
The script implements the method outlined in Section IV of

  W. Tichy and P. Marronetti, "Simple method to set up low eccentricity 
      initial data for moving puncture simulations," Phys. Rev. D 83, 
      024012 (2011).

The Post-Newtonian equations are from 

  L. Kidder, "Coalescing binary systems of compact objects to 
      (post)^{5/2}-Newtonian order. V. Spin effects," 
      Phys. Rev. D 52 821 (1995).

'''


import argparse
import numpy as np
from math import *
import sys




# Set default values for the Gravitational constant and the speed of light.
G = 1.0
c = 1.0

def normal_vector(v):
    '''
    create a normal vector. If v = [ 0, 0, 0], then return the unit vector zhat
    '''
    vmag = np.linalg.norm(v)
    if vmag > 1.0e-9:
        n = v / vmag
    else:
        n = np.array([0,0,1])
    return n


#----------------------------------------------------------------------
#  Set up commandline arguments
#----------------------------------------------------------------------
zd = sqrt(2.0) * 1.0e-5

parser = argparse.ArgumentParser(description='Estimate parameters for puncture initial data.')

#parser.add_argument('-o','--order', type=str, default="2", 
#     help='Post-Newtonian order of the equations', 
#     choices=["0","1","2","3"])
#parser.add_argument('--m1', type=float, default=1.0, 
#       help='Mass of star 1 (Default 1.0)', metavar='NUM')
#parser.add_argument('--m2', type=float, default=1.0, 
#       help='Mass of star 2 (Default 1.0)', metavar='NUM')
parser.add_argument('separation', type=float, 
   help='seperation of the binary star system', metavar = 'SEPARATION')

parser.add_argument('-M', '--total_mass', type=float, default=1.0, 
       help='Total ADM Mass (Default 1.0)', metavar='NUM')
parser.add_argument('-r', '--mass_ratio', type=float, default=1.0, 
       help='Mass ratio r = m1/m2, r >= 1 (Default 1.0)', metavar='NUM')
parser.add_argument('--s1', type=str, default="0, 0, 0", 
     help='BH1 spin "spin parameter, theta, phi". (Default "0,0,0")', 
     metavar='"NUM, NUM, NUM"')
parser.add_argument('--s2', type=str, default="0, 0, 0", 
     help='BH2 spin "spin parameter, theta, phi". (Default "0,0,0")', 
     metavar='"NUM, NUM, NUM"')
parser.add_argument('-z', '--zoffset', type=float, default=zd, 
       help='z coordinate of bh plane (Default sqrt(2)*1e-5)', metavar='NUM')

args = parser.parse_args()

d = args.separation
M = args.total_mass
mr = args.mass_ratio
zoffset = args.zoffset


if d < 6.0:
    print('separation between black holes too small. d >= 6.0')
    sys.exit()

# Calculate the ADM masses of each black hole, m1 and m2.
m2 = M/(mr + 1.0)
m1 = M - m2

mu = m1 * m2 / M
nu = mu / M
GMD = G * M / d

# Calculate the bare masses. 
# These are used for puncture initial data, but the ADM masses are used 
# in the PN expressions. The expression for the bare masses is derived
# from Eq. (22) of Tichy and Marronetti.
mb1 = 0.5*(m1 - m2 - 2.0*d + sqrt( 8.0*m1*d + (m2 - m1 + 2.0*d)**2))
mb2 = 0.5*(m2 - m1 - 2.0*d + sqrt( 8.0*m1*d + (m2 - m1 + 2.0*d)**2))


#
# The spin parameters are spin magnitude, theta, phi. Theta, phi are the
# standard spherical coordinate angles. Use this to construct the spin vectors
# S1 and S2 in Cartesian coordinates.
#
s1pars = np.fromstring(args.s1,count=3,sep=',')
s2pars = np.fromstring(args.s2,count=3,sep=',')

xi1 = s1pars[0]
theta1 = s1pars[1]
phi1 = s1pars[2]
xi2 = s2pars[0]
theta2 = s2pars[1]
phi2 = s2pars[2]
S1norm = xi1 * m1**2
S2norm = xi2 * m2**2

S1 = [ S1norm * sin(theta1) * cos(phi1), S1norm * sin(theta1) * sin(phi1), S1norm * cos(theta1) ]
S2 = [ S2norm * sin(theta2) * cos(phi2), S2norm * sin(theta2) * sin(phi2), S2norm * cos(theta2) ]



'''
Calculate the tangential momentum for the zero-spin case using Eq. (45) from 

   B. Walther, B Bruegmann, and D. Mueller, "Numerical black hole initial 
       data with low eccentricity based on post-Newtonian orbital parameters,"
       arXiv:0901.0993v3 [gr-qc] 2009.

This is now deprecated.

if args.order == "0":
    ptns = mu * sqrt(GMD)
elif args.order == "1":
    ptns = mu * (sqrt(GMD) + 1.0/c**2 * GMD**1.5)
elif args.order == "2":
    ptns = mu * (sqrt(GMD) + 1.0/c**2 * GMD**1.5 + 
               1.0/(16.0 * c**4) * (42 - 43 * nu) * GMD**2.5)
elif args.order == "3":
    ptns = mu * (sqrt(GMD) + 1.0/c**2 * GMD**1.5 + 
               1.0/(16.0 * c**4) * (42 - 43 * nu) * GMD**2.5 + 
               1.0/(128.0*c**6) * (480.0 + (163*pi**2 - 4556)*nu +
                   104*nu**2)*GMD**3.5)
else:
    print('Unknown PN order = ' + args.order)
'''


#
# Calculate the tangential momentum following Tichy & Marronetti.
# Set the Newtonian angular momentum (LNewt) to be along the z-axis.
#
LNewt = np.array([0,0,1])

#
# S1N and S2N are \hat{\bf s}_A in Kidder.
#
S1N = normal_vector(S1)
S2N = normal_vector(S2)

#
# Calculate the PN2 angular momentum L from Eqs. (4.7) and (2.8) in Kidder.
#
T0 = mu*sqrt(M*d)*(1.0 + 2.0*GMD - 0.25*(xi1*np.dot(LNewt,S1N)*(8.0*(m1/M)**2 + 7.0*nu) + xi2 *np.dot(LNewt,S2N)*(8.0*(m2/M)**2 + 7.0*nu))*GMD**(1.5) + (0.5*(5.0 - 9.0*nu) - 0.75*nu*xi1*xi2*(np.dot(S1N,S2N) - 3.0*np.dot(LNewt, S1N)*np.dot(LNewt, S2N)))*GMD**2)
T1 = -0.25*mu*sqrt(M*d)*xi1*(4.0*(m1/M)**2 + nu)*(GMD)**1.5
T2 = -0.25*mu*sqrt(M*d)*xi2*(4.0*(m2/M)**2 + nu)*(GMD)**1.5

L = np.zeros(3)
L = T0*LNewt + T1*S1N + T2*S2N

#
# Rotate L such that it is aligned along the z-axis. The calculation of the
# rotation matrix R is from
# https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
#

LN = normal_vector(L)

zN = np.array([0,0,1.0])   # the normal vector along z.

v = np.cross(LN, zN)
v1 = v[0]
v2 = v[1]
v3 = v[2]
vc = np.dot(LN,zN)

vx = np.array([[0, -v3, v2],[v3, 0, -v1], [-v2, v1, 0]])
vx2 = vx @ vx
T4 = 1.0/(1.0 + vc)
R = np.identity(3) + vx + T4 * vx2   # Rotation matrix

# Rotate L, S1, and S2.
Lrot = R @ L
bh1_spin = R @ S1
bh2_spin = R @ S2

# The tangential momentum pt for the punctures. Eq. (20) in Tichy and 
# Marronetti.
pt = np.linalg.norm(Lrot) / d

# Tichy and Marronetti note that the PN2 parameters for BBH punctures
# will not lead to quasi-circular orbits in full GR. Experience has shown
# that The tangential the largest errors are in the tangential momentum (pt),
# rather than the radial momentum. As an guide to the size of variations
# in pt that may be needed to obtain orbits with a smaller eccentricity, 
# we calculate dpt using Eq. (23)
#
dpt = mu*sqrt(GMD)*(11.29*GMD**3 - 92.37*GMD**4)

#
# Calculate the radial momentum pr. See Eq. (21) in Tichy and Marronetti.
# The equation for rdot is Eq. (4.12) in Kidder.
#
rdot = - 64.0/5.0 * nu * (GMD)**3 * (1.0 - 1.0/336.0*(1751.0 + 588.0*nu)*GMD - (7.0/12.0*(xi1*np.dot(LNewt,S1N)*(19.0*(m1/M)**2 + 15.0*nu) + xi2*np.dot(LNewt,S2N)*(19.0*(m2/M)**2 + 15.0*nu)) - 4.0*pi)*GMD**1.5 - 5.0/48.0*nu*xi1*xi2*(59.0*np.dot(S1N,S2N) - 173.0*np.dot(LNewt,S1N)*np.dot(LNewt,S2N))*GMD**2)
pr = mu*np.absolute(rdot)

#
# Calculate the black hole coordinates and momenta in Cartesian coordinates.
# (See Tichy and Marronetti.)
#
bh1_x = np.array([m2*d/M, 0.0, zoffset])
bh2_x = np.array([-m1*d/M, 0.0, zoffset])
bh1_p = np.array([-pr, pt, 0.0])
bh2_p = np.array([pr, -pt, 0.0])



'''
#  ...Debugging stuff

print('m1  = ' + str(m1))
print('m2  = ' + str(m2))
print('mb1 = ' + str(mb1))
print('mb2 = ' + str(mb2))
print('S1N = ' + str(S1N))
print('S2N = ' + str(S2N))
print('|S1| = ' + str(np.linalg.norm(S1)))
print('|S2| = ' + str(np.linalg.norm(S2)))
print('|bh1_spin| = ' + str(np.linalg.norm(bh1_spin)))
print('|bh2_spin| = ' + str(np.linalg.norm(bh2_spin)))
print('|S1N| = ' + str(np.linalg.norm(S1N)))
print('|S2N| = ' + str(np.linalg.norm(S2N)))
print('L = ' + str(L))
print('LN = ' + str(LN))
print('Lrot = ' + str(Lrot))
print('|L| = ' + str(np.linalg.norm(L)))
print('|Lrot| = ' + str(np.linalg.norm(Lrot)))
'''

print('----------------------------------------------------------------')
print(' Initial data parameters')
print('----------------------------------------------------------------')
print('Total Mass:        ' + str(M))
print('Mass ratio:        ' + str(mr))
print('m1 (ADM):          ' + str(m1))
print('m2 (ADM):          ' + str(m2))
print('separation:        ' + str(d))
print('Spin 1:            ' + str(S1))
print('Spin 2:            ' + str(S2))
print('z offset:          ' + str(zoffset))
print(' ')

print('----------------------------------------------------------------')
print(' PUNCTURE PARAMETERS')
print('----------------------------------------------------------------')
print('BH1')
print('    bare mass: ' + str(mb1))
print('    coords:    ' + str(bh1_x))
print('    momentum:  ' + str(bh1_p))
print('    spin:      ' + str(bh1_spin))

print(' ')

print('BH2')
print('    bare mass: ' + str(mb2))
print('    coords:    ' + str(bh2_x))
print('    momentum:  ' + str(bh2_p))
print('    spin:      ' + str(bh2_spin))
print(' ')

print('----------------------------------------------------------------')
print(' PUNCTURE PARAMETERS (par file foramt)')
print('----------------------------------------------------------------')

print("\"BSSN_BH1\": {")
print("    \"MASS\":%f," %(mb1))
print("    \"X\":%f,"%(bh1_x[0]))
print("    \"Y\":%f,"%(bh1_x[1]))
print("    \"Z\": %f,"%(bh1_x[2]))
print("    \"V_X\": %f,"%(bh1_p[0]))
print("    \"V_Y\": %f,"%(bh1_p[1]))
print("    \"V_Z\": %f,"%(bh1_p[2]))
print("    \"SPIN\": %f,"%(bh1_spin[0]))
print("    \"SPIN_THETA\":%f,"%(bh1_spin[1]))
print("    \"SPIN_PHI\": %f"%(bh1_spin[2]))
print("},")
print(" \"BSSN_BH2\": {")
print("      \"MASS\":%f,"%(mb2))
print("      \"X\":%f," %(bh2_x[0]))
print("      \"Y\":%f," %(bh2_x[1]))
print("      \"Z\":%f," %(bh2_x[2]))
print("      \"V_X\": %f," %(bh2_p[0]))
print("      \"V_Y\": %f," %(bh2_p[1]))
print("      \"V_Z\": %f," %(bh2_p[2]))
print("      \"SPIN\": %f," %(bh2_spin[0]))
print("      \"SPIN_THETA\":%f,"%(bh2_spin[1]))
print("      \"SPIN_PHI\": %f" %(bh2_spin[2]))
print("}")

print('The tangential momentum is just an estimate, and the value for a')
print('for a circular orbit is likely between (' + str(pt-dpt) + ', ' + str(pt+dpt) + ')')



