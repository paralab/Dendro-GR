import ctypes
import numpy as np
import sympy

sph = ctypes.CDLL("./libsph.so")
sph.real_spherical_harmonic.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double]
sph.real_spherical_harmonic.restype  = ctypes.c_double

sph.real_spherical_harmonic_d_theta.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
sph.real_spherical_harmonic_d_theta.restype  = ctypes.c_double

sph.real_spherical_harmonic_d_phi.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
sph.real_spherical_harmonic_d_phi.restype  = ctypes.c_double

# ll, mm, tt, pp = sympy.symbols("l m t p")
# ylm            = sympy.functions.special.spherical_harmonics.Ynm(ll, mm, tt, pp)

# def Ylm(l,m,theta, phi):
#     return ylm.eval(l,m, theta, phi)

def sph_real(l,m, theta, phi):
    """
    computes the real spherical harmonics
    """
    assert abs(m)<=l
    if m==0:
        sph_norm_fac = sympy.sqrt(sympy.Rational((2 * l + 1), 4))/sympy.sqrt(sympy.pi) 
    else:
        sph_norm_fac = (-1)**m * sympy.sqrt(2) * sympy.sqrt(sympy.Rational((2 * l + 1) * sympy.factorial(l-sympy.Abs(m)) , 4 * sympy.factorial(l+sympy.Abs(m)))) /sympy.sqrt(sympy.pi)

    if m < 0:
        return sph_norm_fac * sympy.functions.special.polynomials.assoc_legendre(l, sympy.Abs(m), sympy.cos(theta)) * sympy.sin(sympy.Abs(m) * phi)
    elif m==0:
        return sph_norm_fac * sympy.functions.special.polynomials.assoc_legendre(l, m, sympy.cos(theta)) 
    else:
        return sph_norm_fac * sympy.functions.special.polynomials.assoc_legendre(l, m, sympy.cos(theta)) * sympy.cos( m * phi)

def Ylm(l, m, theta, phi):
    x, y = sympy.symbols('x y')
    return sph_real(l,m, x, y).subs({x:theta, y:phi}).evalf()

def Ylm_d_theta(l,m, theta, phi):
    x, y = sympy.symbols('x y')
    ylm  = sph_real(l,m, x, y)
    return sympy.diff(ylm, x).subs({x:theta, y:phi}).evalf()

def Ylm_d_phi(l,m, theta, phi):
    x, y = sympy.symbols('x y')
    ylm  = sph_real(l,m, x, y)
    return sympy.diff(ylm, y).subs({x:theta, y:phi}).evalf()

    

l=2; mp1=1; mm1=-1;

a1 = Ylm(l, mm1, 0.2 , 0.1)
a2 = sph.real_spherical_harmonic(l, mm1, 0.2, 0.1)
print("a1 = %.12E a2 = %.12E rel = %.12E"%(a1, a2, abs(1-a1/a2)))
assert abs(1-a1/a2) < 1e-14

a1 = Ylm(l, mp1, 0.2 , 0.1)
a2 = sph.real_spherical_harmonic(l, mp1, 0.2, 0.1)
print("a1 = %.12E a2 = %.12E rel = %.12E"%(a1, a2, abs(1-a1/a2)))
assert abs(1-a1/a2) < 1e-14

print("\n\n")
a1 = Ylm_d_phi(l, mm1, 0.2 , 0.1)
a2 = sph.real_spherical_harmonic_d_phi(l, mm1, 0.2, 0.1, 1)
print("a1 = %.12E a2 = %.12E rel = %.12E"%(a1, a2, abs(1-a1/a2)))
assert abs(1-a1/a2) < 1e-14

a1 = Ylm_d_phi(l, mp1, 0.2 , 0.1)
a2 = sph.real_spherical_harmonic_d_phi(l, mp1, 0.2, 0.1, 1)
print("a1 = %.12E a2 = %.12E rel = %.12E"%(a1, a2, abs(1-a1/a2)))
assert abs(1-a1/a2) < 1e-14

print("\n\n")
a1 = Ylm_d_theta(l, mm1, 0.2 , 0.1)
a2 = sph.real_spherical_harmonic_d_theta(l, mm1, 0.2, 0.1, 1)
print("a1 = %.12E a2 = %.12E rel = %.12E"%(a1, a2, abs(1-a1/a2)))
assert abs(1-a1/a2) < 1e-14

a1 = Ylm_d_theta(l, mp1, 0.2 , 0.1)
a2 = sph.real_spherical_harmonic_d_theta(l, mp1, 0.2, 0.1, 1)
print("a1 = %.12E a2 = %.12E rel = %.12E"%(a1, a2, abs(1-a1/a2)))
assert abs(1-a1/a2) < 1e-14