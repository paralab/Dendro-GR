from __future__ import print_function
import numpy as np


from sympy.physics.units import G, kg, s, m, au, year
import scipy.constants 
from numpy import pi
import numpy as np

parsec = 648000/pi * au
megaparsec = 10**6 * parsec
msun = 4*pi**2*au**3/G/year**2

c=scipy.constants.speed_of_light*m/s


def GeometricTime_To_MKS_Time(t_over_m, mass):
  """
  t_over_m = numpy array with t/M in geometric units
  mass = total mass. Should have units recognized by sympy

  Return value is in units of seconds but the units are removed.
  """

  factor = np.float64(((G*mass/c**3/s)).simplify())
  
  return t_over_m * factor

def GeometricStrain_TO_Observer_Strain(r_h_over_m, mass, distance):
  """
  r_h_over_m = numpy array with strain in geometric units
  mass = total mass. Should have units recognized by sympy
  distance = distance to source. Should have units recognized by sympy

  Return value is strain at observer location
  """
  factor = np.float64((G*mass/c**2 / distance).n().simplify())
  return r_h_over_m  * factor
