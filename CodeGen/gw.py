"""
@brief: This file constains the symbolic code to perform gravitational wave
extraction for Dendro-GR framework.

The extraction method is based on "Extraction of Gravitational Waves in Numerical Relativity"
https://arxiv.org/abs/1606.02532 : Chapter 2 , 3, 4.

For Spin Weighted Spherical Harmonics (SWSH) we use the SphericalFunction module based on
https://arxiv.org/pdf/1604.08140.pdf

"""

import dtypes
import dendro

import numpy as np
import spherical_functions as sf
import quaternion

# import quadpy
import math
from sympy.printing.dot import dotprint
from sympy.printing import print_ccode
from sympy import re, im, pprint, linsolve, symbols

from sympy.functions import arg

import sys

import cog

##############################################################################################################
# Note: Traditional spherical coord notation, theta angle betwen point an z axis, phi angle on the xy plane.
################################################################################################################


class GWExtract:
    def __init__(self, N, qprec, lmodes, spin, varNames):
        # number of different radii for psi4 polynomial approx
        self.N = N

        # precision for the lebedev quadrature
        self.qprec = qprec

        # l modes for SWSH
        self.lmodes = lmodes

        self.numLModes = len(self.lmodes)

        # spin values for SWSH (2)
        self.spin = spin

        # computing the SWSH decomposition
        # self.r=dtypes.scalar(varNames[0])
        # self.psi4_s2=dtypes.scalar(varNames[1])
        # self.psi4_lm=dtypes.scalar(varNames[2])

        # part of generating symbolic expression for the polynomial fit.
        # self.V=dtypes.mat(varNames[3],N,N)
        self.vR = dtypes.vec("r", N)  # radius of the extraction
        self.psi4 = dtypes.vec(
            "psi4_", N
        )  # different lm values at the extraction radius
        self.argpsi4 = dtypes.vec(
            "apsi4_", N
        )  # different lm values at the extraction radius
        self.A = dtypes.vec("a", N)  # coefficients for the polynomial fit.
        self.B = dtypes.vec("b", N)  # coefficennts for the phase polynomical fit.

        # fname=sys.path[0]+"/Lebedev/lebedev_" +"{:03d}".format(self.qprec)+".txt"
        fname = "Lebedev/lebedev_" + "{:03d}".format(self.qprec) + ".txt"
        self.lebedev_theata = []
        self.lebedev_phi = []
        self.lebedev_w = []

        deg_radians = math.pi / 180.0
        # read lebedev points
        with open(fname, "r") as f:
            for row in f:
                row = row.strip().split()
                self.lebedev_theata.append(deg_radians * float(row[1]))
                self.lebedev_phi.append(deg_radians * (float(row[0]) + 180.0))
                self.lebedev_w.append(float(row[2]))

    """
    @brief: Evaluate a given SWSH at a given theta, phi angle. 
    Note that we compute the SWSH based on the Wigner matrix defined on the quaternions
    for specified theta and phi. 

    _sY_{lm}(R(\theta,\phi))= (-1)^s \sqrt{\frac{2l+1}{4\pi}} D^{l}_{m,-s}(R(\theta,\phi))

    @param s: spin 
    @param l: l of the SWSH 
    @param m: m of the SWSH
    @param theta: Spherical coord \theta
    @param phi: Spherical coord \phi

    """

    def swsh(self, s, l, m, theta, phi):
        # get the rotor related for Wigner matrix
        # Note: We need to specify the (-theta,-phi direction to match the convenction used in HAD)
        # 10/12/20: David has pointed out this should be theta, phi to match the SWSH used in the NR community. so changed back to theta phi.
        # this is the reason that we have a sign difference between the LazEv code.
        R_tp = quaternion.from_spherical_coords(theta, phi)
        W = sf.Wigner_D_element(R_tp, l, m, -s)
        # print(W)
        return ((-1) ** (s)) * math.sqrt((2 * l + 1) / (4 * math.pi)) * W

    def initVars(self, varnames):
        cog.outl("#define %s %d \n" % (varnames[0], len(self.lebedev_theata)))

        # allocate Lebedev quadrature points
        cog.outl(
            "static const double %s [] = { %s };"
            % (
                varnames[1],
                ",".join(["{:.20f}".format(x) for x in self.lebedev_theata]),
            )
        )
        cog.outl(
            "static const double %s [] = { %s };"
            % (varnames[2], ",".join(["{:.20f}".format(x) for x in self.lebedev_phi]))
        )
        cog.outl(
            "static const double %s [] = { %s };"
            % (varnames[3], ",".join(["{:.20f}".format(x) for x in self.lebedev_w]))
        )

        for l in self.lmodes:
            for m in range(-l, l + 1):
                lm_value = []
                for q in range(0, len(self.lebedev_theata)):
                    lm_value.append(
                        self.swsh(-2, l, m, self.lebedev_theata[q], self.lebedev_phi[q])
                    )

                if m < 0:
                    cog.outl(
                        "static const DendroComplex %s_sp_%d_l_%d_m_minus_%d [] = { %s };"
                        % (
                            varnames[4],
                            2,
                            l,
                            abs(m),
                            ",".join(
                                [
                                    "DendroComplex({:.20f},{:.20f})".format(
                                        x.real, x.imag
                                    )
                                    for x in lm_value
                                ]
                            ),
                        )
                    )
                else:
                    cog.outl(
                        "static const DendroComplex %s_sp_%d_l_%d_m_plus_%d [] = { %s };"
                        % (
                            varnames[4],
                            2,
                            l,
                            abs(m),
                            ",".join(
                                [
                                    "DendroComplex({:.20f},{:.20f})".format(
                                        x.real, x.imag
                                    )
                                    for x in lm_value
                                ]
                            ),
                        )
                    )

        ptr_str = []
        for i in range(0, len(self.lmodes)):
            for m in range(-self.lmodes[i], self.lmodes[i] + 1):
                if m < 0:
                    ptr_str.append(
                        "%s_sp_%d_l_%d_m_minus_%d"
                        % (varnames[4], 2, self.lmodes[i], abs(m))
                    )
                else:
                    ptr_str.append(
                        "%s_sp_%d_l_%d_m_plus_%d"
                        % (varnames[4], 2, self.lmodes[i], abs(m))
                    )

        cog.out(
            "static const DendroComplex* %s[]={%s};" % (varnames[4], ",".join(ptr_str))
        )

    """
    @brief: Computes the polynomial fit for the psi4 extraction at the far radius.     
    """

    def psi4PolyFit(self):
        for r in range(0, self.N):
            for j in range(0, self.N):
                self.V[r, j] = 1 / (self.vR[r] ** (j + 1))

        self.Vinv = self.V ** (-1)
        self.A = self.Vinv * self.psi4
        self.B = self.Vinv * self.argpsi4

        pprint(self.A)
        pprint(self.B)

        # pprint(Vinv*self.psi4)
        # linsolve((self.V,self.psi4),self.A)
        # print(self.A[0])

    """
    Computes the t* propergation time computation. 
    """

    def propTimeCompute(self):
        print("not implemented yet")


# ======================================================================================

# gw=GWExtract(N=4,qprec=9,lmodes=[2,3,4],spin=2,varNames=["r[rIndex]","psi4_s2[sIndex]","psi4_lm[lmIndex]","Van[w]"])
