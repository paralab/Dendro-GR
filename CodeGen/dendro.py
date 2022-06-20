##########################################################################
# module: dendro
# author: Hari Sundar
# email:  hari@cs.utah.edu
#
# python module to generate efficient code for General Relativity.
#
# (c) 2016 University of Utah, All rights reserved.
##########################################################################

from sympy import *
from sympy.tensor.array import *
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.utilities import numbered_symbols
from sympy.printing import print_ccode
from sympy.printing.dot import dotprint
import networkx as nx
import nxgraph
import sympy
import re as regex

import string
import random

# internal variables
undef = symbols('undefined')

metric = undef
inv_metric = undef
C1 = undef
C2 = undef
# C2_spatial
C3 = undef

# first derivative
d = undef
# second derivative
d2s = undef
# advective derivative
ad = undef

# Kreiss-Oliger dissipation operator
kod = undef

one = symbols('one_')
negone = symbols('negone_')

e_i = [0, 1, 2]
e_ij = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
e_ij_offdiag = [(0, 1), (0, 2), (1, 2)] 

Ricci = undef


def d2(i, j, a):
    global d2s
    if (i>j):
        return d2s(j,i,a)
    else:
        return d2s(i,j,a)


##########################################################################
# variable initialization functions
##########################################################################


def scalar(name, idx):
    """
    Create a scalar variable with the corresponding name. The 'name' will be during code generation, so should match the
    variable name used in the C++ code.
    """
    tname = name + idx
    return symbols(tname)


def vec3(name, idx):
    """
    Create a 3D vector variable with the corresponding name. The 'name' will be during code generation, so should match
    the variable name used in the C++ code. The returned variable can be indexed(0,1,2), i.e.,

    b = dendro.vec3("beta")
    b[1] = x^2
    """
    vname = ' '.join([name + repr(i) + idx for i in [0, 1, 2]])
    return symbols(vname)


def sym_3x3(name, idx):
    """
    Create a symmetric 3x3 matrix variables with the corresponding name. The 'name' will be during code generation, so
    should match the variable name used in the C++ code. The returned variable can be indexed(0,1,2)^2, i.e.,

    gt = dendro.sym_3x3("gt")
    gt[0,2] = x^2
    """

    vname = ' '.join([name + repr(i) + idx for i in range(6)])
    m1, m2, m3, m4, m5, m6 = symbols(vname)

    return Matrix([[m1, m2, m3], [m2, m4, m5], [m3, m5, m6]])


def mat_3x3(name, idx):
    """
    Create a 3x3 matrix variables with the corresponding name. The 'name' will be during code generation, so
    should match the variable name used in the C++ code. The returned variable can be indexed(0,1,2)^2, i.e.,

    gt = dendro.sym_3x3("gt")
    gt[0,2] = x^2
    """
    vname = ' '.join([name + repr(i) + idx for i in range(9)])
    m1, m2, m3, m4, m5, m6, m7, m8 , m9 = symbols(vname)
    return Matrix([[m1, m2, m3], [m4, m5, m6], [m7, m8, m9]])



##########################################################################
# derivative related functions
##########################################################################


def set_first_derivative(g):
    """
    Set how the stencil for the first derivative will be called. Here g is a string

    Typically,

    d_i u =  g(i, u)
    """
    global d
    d = Function(g)
    return d


def set_second_derivative(g):
    """
    Set how the stencil for the second derivative will be called. Here g is a string

    Typically,

    d_ij u =  g(i, j, u)
    """
    global d2s
    d2s = Function(g)
    return d2s

def set_advective_derivative(g):
    """
    Set how the stencil for the second derivative will be called. Here g is a string

    Typically,

    ad_i u =  g(i, u)
    """
    global ad
    ad = Function(g)
    return ad

def set_kreiss_oliger_dissipation(g):
    """
    Set how the stencil for Kreiss-Oliger dissipation will be called. Here g is a string.

    Typically,

    kod_i u = g(i, u)
    """
    global kod
    kod = Function(g)
    return kod

# Covariant Derivatives
def DiDj(a):
    """
    Actually this defines two covariant derivatives acting on a scalar.  The 
    derivative in this case is built from the full (non-conformal) metric as 
    C3 is built from the full (non-conformal) metric.  This object is 
    symmetric in both indices.
    """
    global d, C3

    m = Matrix([d2(i, j, a) - sum([C3[l, i, j] * d(l, a) for l in e_i]) for i, j in e_ij])
    return m.reshape(3, 3)


def _Di_Dj(a):
    """
    This defines two covariant derivatives acting on a scalar.
    The use of C2 below, however, suggests that this derivative is built
    from the conformal metric.  Such an operator and term shows up in the
    definition of the Ricci scalar which, in turn shows up in the trace-free
    term in the At evolution equation.  As with DiDj, this object is symmetric
    in both indices when acting on a scalar.
    """
    global d, C2

    m = Matrix([d2(i, j, a) - sum([C2[l, i, j] * d(l, a) for l in e_i]) for i, j in e_ij])
    return m.reshape(3, 3)


# Index Raising
def up_up(A):
    """
    raises both the indices of A, i.e., A_{ij} --> A^{ij}
    """
    global inv_metric

    m = Matrix([sum([inv_metric[i, k]*inv_metric[j, l]*A[k, l] for k, l in e_ij]) for i, j in e_ij])
    return m.reshape(3, 3)

# One index rasing
def up_down(A):
    """
    raises the first index of A, i.e., A_{ij} --> A^i_j
    """
    global inv_metric

    m = Matrix([sum([inv_metric[i, k]*A[k, j] for k in e_i]) for i, j in e_ij])
    return m.reshape(3, 3)


def lie(b, a, weight=0):
    """
    Computes the Lie derivative of a field, a, along the vector b.  Assumes
    the metric has been set.  An optional weight for the field can be
    specified.

    b must be of type dendro.vec3
    a can be scalar, vec3 or sym_3x3

    Computes L_b(v)
    """
    global d, ad

    # e_ij = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]

    if type(b) != tuple:
        raise ValueError('Dendro: The field wrt which the Lie derivative is calculated needs to be vec3.')

    if type(a) == Symbol:
        return sum([b[i] * ad(i, a) for i in e_i]) + weight*a*sum([d(i, b[i]) for i in e_i])
    elif type(a) == tuple:
        return [sum([b[j] * ad(j, a[i]) - a[j] * d(j, b[i]) + weight*a[i]*d(j, b[j]) for j in e_i]) for i in e_i]
    elif type(a) == Matrix:
        m = Matrix([sum([b[k]*ad(k, a[i, j]) + a[i, k]*d(j, b[k]) + a[k, j]*d(i, b[k]) + weight*a[i, j]*d(k, b[k]) for k in e_i]) for i, j in e_ij])
        return m.reshape(3, 3)
    else:
        raise ValueError('Dendro: Unknown type for input field to compute Lie derivative for.')

def kodiss(a):
    """
    Kreiss-Oliger dissipation operator
    """
    global kod

    if type(a) == Symbol:
        return sum( [ kod(i, a) for i in e_i ] )
    elif type(a) == tuple:
        return [ sum ( [ kod(i, a[j]) for i in e_i ] ) for j in e_i ]
    elif type(a) == Matrix:
        return Matrix( [ sum( [ kod(k, a[i, j]) for k in e_i ] ) for i, j in e_ij ]).reshape(3, 3)
    else:
        raise ValueError('Dendro: Unknown type for input to computer kodiss.')


def laplacian(a, chi):
    """
    Computes the laplacian of a scalar function with respect to the 3D metric
    gamma_ij.  Assumes that the conformally rescaled metric (called gt in
    various places) and the conformal factor (chi) is set.  Note that C3 is
    built from the same 3D metric.  The only place that this laplacian is
    used in the bssn equations is in the evolution equation for K and is
    the laplacian of alpha (the lapse).
    """
    global d, metric, C3

    full_metric = metric/chi
    inv_full_metric = simplify(full_metric.inv('ADJ'))

    # this could be optimized a tad as a symmetric x symmetric quantity is being summed over 
    #return ( sum([ inv_full_metric[i, i] * ( d2(i, i, a) - sum([C3[l, i, i] * d(l, a) for l in e_i]) ) for i in e_i ]) + 2*sum([ inv_full_metric[i,j] * ( d2(i,j,a) - sum([C3[l,i,j] * d(l,a) for l in e_i ]) ) for i,j in e_ij_offdiag ]) )
    
    return sum([ inv_full_metric[i, j] * ( d2(i, j, a) - sum([C3[l, i, j] * d(l, a) for l in e_i]) ) for i, j in e_ij])


def laplacian_conformal(a):
    """
    Computes the (conformal) laplacian of a scalar function with respect
    to the tilded or conformally rescaled metric (called gt in various
    places).  We assume the rescaled metric is set as well the conformal
    factor, chi.  Note that C2 is built from the conformally rescaled
    metric.  This (conformal) laplacian is only used in the definition of
    Ricci that shows up in the evolution equation for At (under the trace
    free operation), and even then only in the part that multiplies the
    metric and which will drop out on taking the trace free part.  So, in
    fact, the code could be written to completely ignore this operation
    in the evolution equations themselves.  However, if the constraints
    are included or the full Ricci is needed for another reason, this
    would be needed.
    """
    global d, inv_metric, C2

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    return sum([ inv_metric[i, j] * (d2(i, j, a) - sum([C2[l, i, j] * d(l, a) for l in e_i])) for i, j in e_ij])


def sqr(a):
    """
    Computes the square of the matrix. Assumes metric is set.
    """
    global inv_metric

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    return sum([a[i, j]*sum([inv_metric[i, k] * inv_metric[j, l] * a[k, l] for k in e_i for l in e_i]) for i, j in e_ij])
    
    #return ( sum([a[i,i]*sum([inv_metric[i, k] * inv_metric[i, l] * a[k, l] for k in e_i for l in e_i]) for i in e_i]) + 2*sum([a[i,j]*sum([inv_metric[i, k] * inv_metric[j, l] * a[k, l] for k in e_i for l in e_i]) for i,j in e_ij_offdiag]) ) 


def trace_free(x):
    """
    makes the operator trace-free
    """
    global metric, inv_metric

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    trace = sum([ inv_metric[i, j] * x[i, j] for i, j in e_ij])

    # X_{ab} - 1/3 gt_{ab} X.
    # tf = Matrix([x[i, j] - 1/3*metric[i,j]*trace for i, j in e_ij])
    tf = Matrix([x[i, j] - metric[i,j]*trace/3 for i, j in e_ij])

    return tf.reshape(3, 3)

def vec_j_del_j(b, a):
    """
    expands to  $\beta^i\partial_i \alpha$
    """
    return sum([b[i]*d(i, a) for i in e_i])


#[ewh] Adding this as this term needs to be in the beta equation as an
#      advective derivative ... and not as a regular (partial) derivative.
def vec_j_ad_j(b, f):
    """
    expands to  $\beta^i\partial_i f$
    """
    return sum([b[i]*ad(i, f) for i in e_i])

    #vec_k_del_k = vec_j_del_j

##########################################################################
# metric related functions
##########################################################################


def set_metric(g):
    """
    sets the metric variable, so that dendro knows how to compute the derived variables. This should be done fairly
    early on. e.g.,

    gt = dendro.sym_3x3("gt")
    dendro.set_metric(gt)
    """
    global metric

    metric = g


def get_inverse_metric():
    """
    Computes and returns the inverse metric. The variables need for be defined in advance. e.g.,

    gt = dendro.sym_3x3("gt")
    dendro.set_metric(gt)
    igt = dendro.get_inverse_metric()
    """
    global metric, inv_metric, undef

    if metric == undef:
        raise ValueError('Dendro: Metric not defined.')

    if inv_metric == undef:
        # method : ('GE', 'LU', or 'ADJ')
        inv_metric = simplify(metric.inv('ADJ'))

    return inv_metric


def get_first_christoffel():
    """
    Computes and returns the first Christoffel symbols (the quantity with three
    lower indices). It assumes the metric has been set. e.g.,

    dendro.set_metric(gt);

    Note, this is used to calculate the second Christoffel symbols as well 
    as the 3D (conformal) Ricci tensor.  This should be built from the 
    conformal metric, called "gamma tilde" or gt.  If the usual (non-conformal)
    metric gets used, this will be wrong.  

    C1 = dendro.get_first_christoffel();
    """
    global metric, inv_metric, undef, C1, d

    if inv_metric == undef:
        get_inverse_metric()

    if C1 == undef:
        C1 = MutableDenseNDimArray(range(27), (3, 3, 3))

        for k in e_i:
            for j in e_i:
                for i in e_i:
                    C1[k, i, j] = 0.5 * (d(j, metric[k, i]) + d(i, metric[k, j]) - d(k, metric[i, j]))

    return C1


def get_second_christoffel():
    """
    Computes and returns the second Christoffel symbols. Assumes the metric 
    has been set. Will compute the first Christoffel if not already 
    computed. e.g.,

    dendro.set_metric(gt);

    Note that this assumes that the metric that comes in is the conformal 
    metric, gt.  This calculates the "usual" Christoffel symbols that we 
    usually refer to.  This includes the factor with the inverse metric.  
    The final quantity has one index up and two indices down.  

    C2 = dendro.get_second_christoffel();
    """
    global C2, C1, inv_metric

    if C2 == undef:
        if C1 == undef:
            get_first_christoffel()

        igt_t = Array(inv_metric, (3, 3))
        C2 = tensorcontraction(tensorproduct(igt_t, C1), (1, 2))

    return C2


def get_complete_christoffel(chi):
    """
    Computes and returns the second Christoffel symbols (despite the 
    returned name of C3 -- sort of a bad choice).  These (in comparison 
    to C2) are built from the non-conformal metric using the conformal 
    metric ("gamma tilde") and the conformal factor, chi. Assumes the 
    metric (i.e. the conformal metric) has been set and pulls in the 
    conformal factor, chi. Will compute the first/second Christoffels if 
    not already computed. e.g.,

    dendro.set_metric(gt);

    C2_spatial = dendro.get_complete_christoffel(); 
    """
    global metric, inv_metric, undef, C1, C2, C3, d

    if C3 == undef:
        C3 = MutableDenseNDimArray(range(27), (3, 3, 3))

        if C2 == undef:
            get_second_christoffel()

        for k in e_i:
            for j in e_i:
                for i in e_i:
                    C3[i, j, k] = C2[i, j, k] - 0.5/(chi)*(KroneckerDelta(i, j) * d(k, chi) +
                                                           KroneckerDelta(i, k) * d(j, chi) -
                                                           metric[j, k]*sum([inv_metric[i, m]*d(m, chi) for m in e_i])
                                                           )

    return C3


def compute_ricci(Gt, chi):
    """
    Computes the (3D) Ricci tensor. e.g.,

    dendro.set_metric(gt)

    R = dendro.compute_ricci(Gt, chi)

    or

    dendro.compute_ricci(Gt, chi)

    and use

    dendro.ricci

    The conformal connection coefficient and the conformal variable needs 
    to be supplied.
    """
    global metric, inv_metric, C1, C2

    Lchi = laplacian_conformal(chi)

    #print(type(Lchi))

    #print('Done with Lphi') #simplify(Lchi))


    #ewh4 DKchiDkchi = Matrix([4*metric[i, j]*sum([sum([inv_metric[k, l]*d(l, chi) for l in e_i])*d(k, chi) for k in e_i]) for i, j in e_ij])
    DKchiDkchi = Matrix([0.25/chi/chi*metric[i, j]*sum([sum([inv_metric[k, l]*d(l, chi) for l in e_i])*d(k, chi) for k in e_i]) for i, j in e_ij])

    #print('done with DKchi') # simplify(DKchiDkchi))

    CalGt = [sum(inv_metric[k,l]*C2[i,k,l] for k, l in e_ij) for i in e_i]

    Rt = Matrix([-0.5*sum([inv_metric[l, m]*d2(l, m, metric[i, j]) for l, m in e_ij]) +
              0.5*sum([metric[k,i]*d(j, Gt[k]) + metric[k,j]*d(i, Gt[k]) for k in e_i]) +
              0.5*sum([CalGt[k]*(C1[i,j,k] + C1[j,i,k]) for k in e_i]) +
              sum([inv_metric[l,m]*(C2[k,l,i]*C1[j,k,m] + C2[k,l,j]*C1[i,k,m] + C2[k,i,m]*C1[k,l,j])
                   for k in e_i for l,m in e_ij]) for i,j in e_ij])

    #print('done with Rt') #simplify(Rt))

    #ewh5    Rphi_tmp = Matrix([2*metric[i, j]*Lchi - 4*d(i, chi)*d(j, chi) for i, j in e_ij])
    #dwn    Rphi_tmp = Matrix([ 0.5*metric[i, j]*Lchi/chi - 0.25*d(i, chi)*d(j, chi)/chi/chi for i, j in e_ij])

    #print(simplify(Rphi_tmp))

    #ewh6    Rphi = -2*_Di_Dj(chi) - Rphi_tmp.reshape(3, 3) - DKchiDkchi.reshape(3, 3)
    #dwn    Rphi = -0.5*_Di_Dj(chi)/chi - Rphi_tmp.reshape(3, 3) - DKchiDkchi.reshape(3, 3)
    xRphi = Matrix( [ 1/(2*chi)*(d2(i,j,chi) -
          sum(C2[k,j,i]*d(k,chi) for k in e_i)) -
          1/(4*chi*chi)*d(i,chi)*d(j,chi) for i, j in e_ij]).reshape(3,3)

    Rphi = xRphi + Matrix( [
           1/(2*chi)*metric[i,j] * ( sum(inv_metric[k,l]*(d2(k,l,chi) -
           3/(2*chi)*d(k,chi)*d(l,chi))  for k, l in e_ij) -
           sum(CalGt[m]*d(m,chi) for m in e_i))
           for i, j in e_ij ] ).reshape(3,3)

    return [Rt.reshape(3, 3) + Rphi, Rt.reshape(3,3), Rphi, CalGt]

##########################################################################
# code generation function
##########################################################################

# construct the common sub-expression ellimination tree
def construct_cse(ex,vnames,idx):
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i]+repr(j)+idx)
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i]+midx[j]+idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i]+idx)
    
    ee_name = 'DENDRO_' #''.join(random.choice(string.ascii_uppercase) for _ in range(5))
    ee_syms = numbered_symbols(prefix=ee_name)
    _v = cse(lexp, symbols=ee_syms, optimizations="basic")
    
    return [_v,count_ops(lexp)]



def generate_cpu(ex, vnames, idx):
    """
    Generate the C++ code by simplifying the expressions.
    """
    # generate_code_nx(ex,vnames,idx)
    # return
    # print(ex)
    custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i]+repr(j)+idx)
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i]+midx[j]+idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i]+idx)

    cse=construct_cse(ex,vnames,idx)
    _v=cse[0]
    
    print("// Dendro: {{{ ")
    print("// Dendro: original ops: %d " %(cse[1]))

    ee_name = 'DENDRO_' 
    ee_syms = numbered_symbols(prefix=ee_name)
    

    rops=0
    print('// Dendro: printing temp variables')
    for (v1, v2) in _v[0]:
        print('const double ', end='')
        print(change_deriv_names(ccode(v2, assign_to=v1, user_functions=custom_functions)))
        rops = rops + count_ops(v2)

    print()
    print('// Dendro: printing variables')
    for i, e in enumerate(_v[1]):
        print("//--")
        print(change_deriv_names(ccode(e, assign_to=lname[i], user_functions=custom_functions)))
        rops = rops + count_ops(e)

    print('// Dendro: reduced ops: %d' %(rops))
    print('// Dendro: }}} ')
    

def generate_cpu_no_cse(ex, vnames, idx):
    """
    Generate the C++ code by simplifying the expressions.
    """
    # print(ex)
    custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i]+repr(j)+idx)
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i]+midx[j]+idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i]+idx)

    for i, e in enumerate(lexp):
        print(change_deriv_names(ccode(expand(e), assign_to=lname[i], user_functions=custom_functions)))


def generate_fpcore(ex, vnames, idx):

    """
    Generate the FPCore code, 
    """
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i]+repr(j)+idx)
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i]+midx[j]+idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i]+idx)

    cse=construct_cse(ex,vnames,idx)
    _v=cse[0]
    
    #print("// Dendro: {{{ ")
    #print("// Dendro: original ops: %d " %(cse[1]))

    ee_name = 'DENDRO_' 
    ee_syms = numbered_symbols(prefix=ee_name)
    
    custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}
    rops=0
    
    #re_symbol=regex.compile(r"Symbol\('[a-z,A-Z,_]+[0-9,\[pp\],\[0-9\]]*'\)")
    re_symbol=regex.compile(r"Symbol\('([a-z,A-Z,0-9,_,\[\]]*)'\)")
    re_integer=regex.compile(r"Integer\(([\-,0-9]+)\)")
    re_float=regex.compile(r"Float\('([\-,0-9]*\.[0-9]*)'\s prec=([0-9]+)\)")                
    re_grad=regex.compile(r"Function\('([a-z]+[0-9]*)'\)\(Integer\(([0-9]+)\),\s*Symbol\('([a-z,A-Z]+[0-9]*\[pp\])'\)\)")

    subs_functions={"Add(":"(+ ", "Integer(-1)": "-1 ", "Mul(":"(* ", "Div(":"(/ ", "Pow(":"(pow ","Rational(":"(/ "}

    #print('// Dendro: printing temp variables')
    tmp_vars=list()
    for (v1, v2) in _v[0]:
        tmp_vars.append(str(v1))
        sym_sub=dict()
        srep=srepr(v2)
        #print(srep)

        res=re_grad.findall(srep)
        for g in res:
            s="Function('%s')(Integer(%s), Symbol('%s'))" %(g[0],g[1],g[2])
            #print(s)
            ss="Symbol('%s')"%(g[0]+"_"+g[1]+"_"+g[2])
            srep=srep.replace(s,ss)

        srep=srep.replace(","," ")            
        #print(srep)

        res=re_symbol.findall(srep)
        inp_params=list()
        #print(res)
        for s in res:
            ss=s.replace("[pp]","")
            for index in range(0,6):
                ss=ss.replace("["+str(index)+"]",str(index))
            inp_params.append(ss)
            tmp_vars.append(ss)
            sym_sub["Symbol(\'%s\')"%(s)]=ss

        int_sub=dict()
        res=re_integer.findall(srep)
        for s in res:
            int_sub["Integer(%s)"%(s)]=s

        float_sub=dict()
        res=re_float.findall(srep)
        
        for s in res:
            float_sub["Float('%s'  prec=%s)"%(s[0],s[1])]=s[0]

        for key,val in sym_sub.items():
            #print("{%s: %s} "%(key,val))
            srep=srep.replace(key,val)

        for key,val in int_sub.items():
            #print("{%s: %s} "%(key,val))
            srep=srep.replace(key,val)

        for key,val in float_sub.items():
            #print("{%s: %s} "%(key,val))
            srep=srep.replace(key,val)            

        for key, val in subs_functions.items():
            srep=srep.replace(key,val)

        
        print("(FPCore (%s)" %(" ".join(inp_params)))
        print("\t%s" %(srep))
        print(")\n")
    
    #print(tmp_vars)
    tmp_vars.clear();
    tmp_vars=list()
    for i, e in enumerate(_v[1]):
        srep=srepr(e)

        res=re_grad.findall(srep)
        for g in res:
            s="Function('%s')(Integer(%s), Symbol('%s'))" %(g[0],g[1],g[2])
            #print(s)
            ss="Symbol('%s')"%(g[0]+"_"+g[1]+"_"+g[2])
            srep=srep.replace(s,ss)

        srep=srep.replace(","," ")            
        
        res=re_symbol.findall(srep)
        inp_params=list()
        #print(res)
        for s in res:
            ss=s.replace("[pp]","")
            for index in range(0,6):
                ss=ss.replace("["+str(index)+"]",str(index))
            inp_params.append(ss)
            tmp_vars.append(ss)
            sym_sub["Symbol(\'%s\')"%(s)]=ss

        int_sub=dict()
        res=re_integer.findall(srep)
        for s in res:
            int_sub["Integer(%s)"%(s)]=s
        
        float_sub=dict()
        res=re_float.findall(srep)

        for s in res:
            float_sub["Float('%s'  prec=%s)"%(s[0],s[1])]=s[0]

        for key,val in sym_sub.items():
            #print("{%s: %s} "%(key,val))
            srep=srep.replace(key,val)

        for key,val in int_sub.items():
            #print("{%s: %s} "%(key,val))
            srep=srep.replace(key,val)

        for key,val in float_sub.items():
            #print("{%s: %s} "%(key,val))
            srep=srep.replace(key,val)

        for key, val in subs_functions.items():
            srep=srep.replace(key,val)


        tmp_vars=list(set(tmp_vars))
        print("(FPCore (%s)" %(" ".join(tmp_vars)))
        print("\t%s" %(srep))
        print(")")
        #print(")")
        #print(")")
        #print(change_deriv_names(ccode(e, assign_to=lname[i], user_functions=custom_functions)))
    


def generate_avx(ex, vnames, idx):
    
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i]+repr(j)+idx)
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i]+midx[j]+idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i]+idx)
    
    cse=construct_cse(ex,vnames,idx)
    _v=cse[0]
    
    print('// Dendro: {{{ ')
    print("// Dendro: original ops: %d " %(cse[1]))

    ee_name = 'DENDRO_' 
    ee_syms = numbered_symbols(prefix=ee_name)
    
    print('// Dendro vectorized code: {{{')
    oper = {'mul': 'dmul', 'add': 'dadd', 'load': '*'}
    prevdefvars = set()
    for (v1, v2) in _v[0]:
        vv = numbered_symbols('v')
        vlist = []
        gen_vector_code(v2, vv, vlist, oper, prevdefvars, idx)
        print('  double ' + repr(v1) + ' = ' + repr(vlist[0]) + ';')
    for i, e in enumerate(_v[1]):
        print("//--")
        vv = numbered_symbols('v')
        vlist = []
        gen_vector_code(e, vv, vlist, oper, prevdefvars, idx)
        #st = '  ' + repr(lname[i]) + '[idx] = ' + repr(vlist[0]) + ';'
        st = '  ' + repr(lname[i]) + " = " + repr(vlist[0]) + ';'
        print(st.replace("'",""))

    print('// Dendro vectorized code: }}} ')



def change_deriv_names(str):
    c_str=str
    derivs=['agrad','grad','kograd']
    for deriv in derivs:
        key=deriv+'\(\d, \w+\[pp\]\)'
        slist=regex.findall(key,c_str)
        for s in slist:
            #print(s)
            w1=s.split('(')
            w2=w1[1].split(')')[0].split(',')
            #print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep=w1[0]
            for v in w2:
                rep=rep+'_'+v.strip()
            #rep=rep+';'
            c_str=c_str.replace(s,rep)

    derivs2=['grad2']
    for deriv in derivs2:
        key=deriv+'\(\d, \d, \w+\[pp\]\)'
        slist=regex.findall(key,c_str)
        for s in slist:
            #print(s)
            w1=s.split('(')
            w2=w1[1].split(')')[0].split(',')
            #print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep=w1[0]
            for v in w2:
                rep=rep+'_'+v.strip()
            #rep=rep+';'
            c_str=c_str.replace(s,rep)
    
    func_list=['pow']
    for func in func_list:
        key=func+'\(\w+, \d\)'
        slist=regex.findall(key,c_str)
        for s in slist:
            #print(s)
            w1=s.split('(')
            w2=w1[1].split(')')[0].split(',')
            #print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            if(int(w2[1].strip())==2):
                rep = "("+w2[0].strip()+" * "+w2[0].strip()+")"
                c_str=c_str.replace(s,rep)
            # rep=w1[0]
            # for v in w2:
            #     rep=rep+'_'+v.strip()
            # #rep=rep+';'
            # c_str=c_str.replace(s,rep)
    return c_str



def generate_separate(ex, vnames, idx,prefix=""):
    """
    Generate the C++ code by simplifying the expressions.
    """
    # print(ex)
    if len(ex)!=1 :
        print ('pass each variable separately ',end='\n')
        return

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i]+repr(j)+idx)
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i]+midx[j]+idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i]+idx)

    # print(num_e)
    # print(len(lname))
    c_file=open(prefix + vnames[0]+'.cpp','w')
    print('generating code for '+vnames[0])
    print('    bssn::timer::t_rhs.start();',file=c_file)
    print('for (unsigned int k = 3; k < nz-3; k++) { ',file=c_file)
    print('    z = pmin[2] + k*hz;',file=c_file)

    print('for (unsigned int j = 3; j < ny-3; j++) { ',file=c_file)
    print('    y = pmin[1] + j*hy; ',file=c_file)

    print('for (unsigned int i = 3; i < nx-3; i++) {',file=c_file)
    print('    x = pmin[0] + i*hx;',file=c_file)
    print('    pp = i + nx*(j + ny*k);',file=c_file)
    print('    r_coord = sqrt(x*x + y*y + z*z);',file=c_file)
    print('    eta=ETA_CONST;',file=c_file)
    print('    if (r_coord >= ETA_R0) {',file=c_file)
    print('    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);',file=c_file)
    print('    }',file=c_file)




    print('// Dendro: {{{ ',file=c_file)
    print('// Dendro: original ops: ', count_ops(lexp),file=c_file)

    # print("--------------------------------------------------------")
    # print("Now trying Common Subexpression Detection and Collection")
    # print("--------------------------------------------------------")

    # Common Subexpression Detection and Collection
    # for i in range(len(ex)):
    #     # print("--------------------------------------------------------")
    #     # print(ex[i])
    #     # print("--------------------------------------------------------")
    #     ee_name = ''.join(random.choice(string.ascii_uppercase) for _ in range(5))
    #     ee_syms = numbered_symbols(prefix=ee_name)
    #     _v = cse(ex[i],symbols=ee_syms)
    #     # print(type(_v))
    #     for (v1,v2) in _v[0]:
    #         print("double %s = %s;" % (v1, v2))
    #     print("%s = %s" % (vnames[i], _v[1][0]))

    #mex = Matrix(ex)
    ee_name = 'DENDRO_' #''.join(random.choice(string.ascii_uppercase) for _ in range(5))
    ee_syms = numbered_symbols(prefix=ee_name)
    _v = cse(lexp, symbols=ee_syms, optimizations='basic')

    custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}

    rops=0
    print('// Dendro: printing temp variables',file=c_file)
    for (v1, v2) in _v[0]:
        # print("double %s = %s;" % (v1, v2)) # replace_pow(v2)))
        print('double ', end='', file=c_file)
        print(change_deriv_names(ccode(v2, assign_to=v1, user_functions=custom_functions)),file=c_file)
        rops = rops + count_ops(v2)


    print('// Dendro: printing variables',file=c_file)
    for i, e in enumerate(_v[1]):
        print("//--",file=c_file)
        # print("%s = %s;" % (lname[i], e)) # replace_pow(e)))
        f = open(str(vnames[0])+'.gv','w')
        print(dotprint(e), file=f)
        f.close()
        print(change_deriv_names(ccode(e, assign_to=lname[i], user_functions=custom_functions)),file=c_file)
        #c_file.write('\n')
        rops = rops + count_ops(e)

    print('// Dendro: reduced ops: ', rops,file=c_file)
    print('// Dendro: }}} ',file=c_file)


    print('  }',file=c_file)
    print(' }',file=c_file)
    print('}',file=c_file)
    print('     bssn::timer::t_rhs.stop();',file=c_file)
    c_file.close()
    print('generating code for '+vnames[0]+' completed')


def replace_pow(exp_in):
    """
    Convert integer powers in an expression to Muls, like a**2 => a*a
    :param exp_in: the input expression,
    :return: the output expression with only Muls
    """
    pows = list(exp_in.atoms(Pow))
    if any(not e.is_Integer for b, e in (i.as_base_exp() for i in pows)):
         raise ValueError("Dendro: Non integer power encountered.")
    repl = zip(pows, (Mul(*[b]*e, evaluate=False) for b, e in (i.as_base_exp() for i in pows)))
    return exp_in.xreplace(dict(repl))


def generate_debug (ex, vnames):
    """
    Generate the C++ code by simplifying the expressions.
    """
    # print(ex)

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    print('// Dendro: {{{ ')
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                # lexp.append(ev)
                print(vnames[i] + repr(j), end='')
                print(' = ', end='')
                print(replace_pow(ev), ';')
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                # lexp.append(e[k])
                print(vnames[i] + midx[j], end='')
                print(' = ', end='')
                print(replace_pow(e[k]), ';')
        else:
            num_e = num_e + 1
            #lexp.append(e)
            print(vnames[i], end='')
            print(' = ', end='')
            print(replace_pow(e), ';')

    print('// Dendro: }}} ')

def vec_print_str(tv, pdvars):
    """
    This returns a string that will be used to print a line of code. If the
    variable tv has not yet been used before, then the declaration of this
    variable must be included in the string. pdvars is the list of variables
    that have been previously defined.

        tv:          new temporary variable
        pdvars:      list of previously declared variables.
    """
    st = '  '
    if tv not in pdvars:
        st += 'double '
        pdvars.add(tv)
    return st

def gen_vector_code(ex, vsym, vlist, oper, prevdefvars, idx):
    """
    create vectorized code from an expression.
    options:
        ex:               expression
        vsym:             numbered symbols
        vlist:            an empty list that is used to process the tree. on return
                          this list contains the name of the variable with the final
                          result
        oper:             dictionary for '+' and '*' operators
        prevdefvars:      an empty set used to identify previously defined temporary variables.
        idx:              name of index for accessing arrays, i.e., alpha[idx].
    """
    one = symbols('one')
    negone = symbols('negone')
    #print (vlist)
    if isinstance(ex, Function):
        # check to see if we are processing a derivative
        if isinstance(ex, ad) or isinstance(ex, d) or isinstance(ex, kod) or isinstance(ex,d2s):
            #print('...ex and args: ',ex,ex.func,ex.args)
            tv = next(vsym)
            vlist.append(tv)
            st = vec_print_str(tv, prevdefvars)
            str_args = [repr(a) for a in ex.args]
            o1 = oper['load']
            o1s = repr(o1).replace("'","")
            idxn = idx.replace("[","")
            idxn = idxn.replace("]","")
            st += repr(tv) + ' = ' + o1s + '(' + repr(ex.func) + '_' + '_'.join(str_args) + '+' + idxn + ' );'
            # st += repr(tv) + ' = ' + repr(ex) + ';'
            print(st.replace(idx,""))
            return

    if isinstance(ex, Pow):
        # check to see if we are processing a simple pow
        a1, a2 = ex.args
        #print('processing pow...',ex,a1,a2)
        if isinstance(a1, Symbol) and isinstance(a2, Number):
            # This is a simple Pow function. Process it here and return
            tv = next(vsym)
            vlist.append(tv)
            st = vec_print_str(tv, prevdefvars)
            if (a2 == -1):
                st += repr(tv) + ' = 1.0 / ' + repr(a1) + ';'
            elif (a2 == 2):
                st += repr(tv) + ' = ' + repr(a1) + ' * ' + repr(a1) + ';'
            else:
                st += repr(tv) + ' = pow( ' + repr(a1) + ', ' + repr(a2) + ');'
            print(st)
            return

    # recursively process the arguments of the function or operator
    for arg in ex.args:
        gen_vector_code(arg, vsym, vlist, oper, prevdefvars, idx)

    if isinstance(ex, Number):
        if isinstance(ex, Integer) and ex == 1:
            vlist.append(one)
        elif isinstance(ex, Number) and ex == -1:
            vlist.append(negone)
        else:
            tv = next(vsym)
            vlist.append(tv)
            st = vec_print_str(tv, prevdefvars)
            if isinstance(ex, Rational):
                st += repr(tv) + ' = ' + repr(float(ex)) + ';'
            else:
                st += repr(tv) + ' = ' + repr(ex) + ';'
            print(st)
    elif isinstance(ex, Symbol):
        tv = next(vsym)
        vlist.append(tv)
        st = vec_print_str(tv, prevdefvars)
        st += repr(tv) +  ' = ' + repr(ex) + ';'
        print(st)
    elif isinstance(ex, Mul):
        nargs = len(ex.args)
        #print('mul..',len(vlist))
        for i in range(nargs-1):
            tv = next(vsym)
            st = vec_print_str(tv, prevdefvars)
            st += repr(tv) + ' = '
            v1 = vlist.pop()
            v2 = vlist.pop()
            #st += repr(v1) + ' * ' + repr(v2) + ';'
            o1 = oper['mul']
            st += repr(o1) + '(' + repr(v1) + ', ' + repr(v2) + ');'
            print(st.replace("'", ""))
            vlist.append(tv)
    elif isinstance(ex, Add):
        nargs = len(ex.args)
        #print('add..',len(vlist))
        for i in range(nargs-1):
            tv = next(vsym)
            st = vec_print_str(tv, prevdefvars)
            st += repr(tv) + ' = '
            v1 = vlist.pop()
            v2 = vlist.pop()
            o1 = oper['add']
            st += repr(o1) + '(' + repr(v1) + ', ' + repr(v2) + ');'
            print(st.replace("'",""))
            vlist.append(tv)
    elif isinstance(ex, Pow):
        tv = next(vsym)
        qexp = vlist.pop()
        qman = vlist.pop()
        a1, a2 = ex.args
        o1 = oper['mul']
        if isinstance(a2,Integer):
            if (a2 == -1):
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' =  1.0 / ' + repr(qman) + ';'
            elif (a2 == 2):
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = ' + repr(o1) + '(' + repr(qman) + ', ' + repr(qman) + ');'
            elif (a2 == -2):
                v1 = next(vsym)
                st = vec_print_str(v1, prevdefvars)
                st += repr(v1) + ' = ' + repr(o1) + '(' + repr(qman) + ', ' + repr(qman) + ');'
                print(st.replace("'",""))
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = 1.0 / ' + repr(v1) + ';'
            elif (a2 > 2 and a2 < 8):
                v1 = next(vsym)
                st = vec_print_str(v1, prevdefvars)
                st += repr(v1) + ' = ' + repr(o1) + '(' + repr(qman) + ', ' + repr(qman) + ');'
                print(st.replace("'",""))
                for i in range(a2-3):
                    v2 = next(vsym)
                    st = vec_print_str(v2, prevdefvars)
                    st += repr(v2) + ' = ' + repr(o1) + '(' + repr(v1) + ', ' + repr(qman) + ');'
                    print(st.replace("'",""))
                    v1 = v2
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = ' + repr(o1) + '(' + repr(v1) + ', ' + repr(qman) + ');'
            else:
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = pow(' + repr(qman) + ',' + repr(qexp) + ');'
        else:
            st = vec_print_str(tv, prevdefvars)
            st = repr(tv) + ' = pow(' + repr(qman) + ',' + repr(qexp) + ');'
        print(st.replace("'",""))
        vlist.append(tv)


def store_node(v, at_idx, local_mem):
    for key,val in local_mem.items():
        if val is None:
            local_mem[key]=v
            #print("storing node %s at %d"%(v,key))
            break
    
    at_idx[v]=key

def evict_node(v,at_idx,local_mem):
    local_mem[at_idx[v]]=None
    #print("release %s"%v)
    at_idx[v]=-1

def visit_node(G: nx.DiGraph, v, work_queue, local_mem, is_allocated):
    at_eval    = nx.get_node_attributes(G,"eval")
    at_func    = nx.get_node_attributes(G,"func") 
    at_args    = nx.get_node_attributes(G,"args") 
    at_idx     = nx.get_node_attributes(G, "idx")
    at_eval[v] = True

    descendents_list = list(G.successors(v))

    if( at_func[v] != sympy.core.add.Add  and at_func[v]!=sympy.core.mul.Mul):
        """
        Direct evaluation for non-reduction type functions such as pow
        """
        #print(at_func[v],v)
        if(at_func[v] == sympy.core.power.Pow):
            a1=at_args[v][0]
            a2=at_args[v][1]

            # print(a1,a2)
            # print("a1 is at DENDRO_%d\n"%at_idx[a1])
            # print(type(a2))
            
            if(at_idx[a1]==-1):
                assert False, "invalid traversal order"

            if(type(a2) == sympy.core.numbers.Integer or type(a2) == sympy.core.numbers.One or type(a2) == sympy.core.numbers.NegativeOne):
                if(at_idx[v]==-1):
                    store_node(v,at_idx,local_mem)

                if(is_allocated[at_idx[v]] == True):
                    print("DENDRO_%d = DENDRO_%d;\n"%(at_idx[v],at_idx[a1])) 
                else:
                    print("double DENDRO_%d = DENDRO_%d;\n"%(at_idx[v],at_idx[a1])) 
                    is_allocated[at_idx[v]]=True
                
                
                for i in range(abs(int(a2))-1):
                    print("DENDRO_%d *= DENDRO_%d;"%(at_idx[v],at_idx[a1])) 
                    
                if(int(a2)<0):
                    print("DENDRO_%d = 1/DENDRO_%d;"%(at_idx[v],at_idx[v])) 


            G.remove_edge(v,a1)
            G.remove_edge(v,a2)

            if(G.in_degree(a1)==0):
                evict_node(a1,at_idx,local_mem)

            if(G.in_degree(a2)==0):
                evict_node(a2,at_idx,local_mem)

        
        else:
            c_code = ccode(v)
            if(at_idx[v]==-1):
                store_node(v,at_idx,local_mem)
                if(is_allocated[at_idx[v]] == True):
                    print("DENDRO_%d = %s;"%(at_idx[v],c_code)) 
                else:
                    print("double DENDRO_%d = %s;"%(at_idx[v],c_code)) 
                    is_allocated[at_idx[v]]=True
    else:
        for u in descendents_list:
            if(at_eval[u] == True):
                if(at_idx[v]==-1):
                    store_node(v,at_idx,local_mem)
                    #print("\n// initialize reduce for %s"%v)
                    print("\n// initialize reduction for ")
                    if(at_func[v] == sympy.core.add.Add):
                        if(is_allocated[at_idx[v]] == True):
                            print("DENDRO_%d = 0;\n"%(at_idx[v])) 
                        else:
                            print("double DENDRO_%d = 0;\n"%(at_idx[v])) 
                            is_allocated[at_idx[v]]=True

                    elif(at_func[v] == sympy.core.mul.Mul):
                        if(is_allocated[at_idx[v]] == True):
                            print("DENDRO_%d = 1;\n"%(at_idx[v])) 
                        else:
                            print("double DENDRO_%d = 1;\n"%(at_idx[v]))
                            is_allocated[at_idx[v]]=True
                    
                if(at_func[v] == sympy.core.add.Add):
                    print("DENDRO_%d += DENDRO_%d;"%(at_idx[v],at_idx[u])) 
                elif(at_func[v] == sympy.core.mul.Mul):
                    print("DENDRO_%d *= DENDRO_%d;"%(at_idx[v],at_idx[u])) 
                
                G.remove_edge(v,u)
                if G.in_degree(u)==0:
                    evict_node(u,at_idx,local_mem)

                    
            else:
                work_queue.append(u)
                at_eval[v]=False
    
    nx.set_node_attributes(G, at_eval, "eval")
    nx.set_node_attributes(G, at_idx,  "idx")
    return G,local_mem


def generate_code_nx(ex, vnames, idx):
    """
    Generate the C++ code by simplifying the expressions.
    """
    # print(ex)
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i]+repr(j)+idx)
        elif type(e) == Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i]+midx[j]+idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i]+idx)
    


    g=nxgraph.ExpressionGraph()
    g.add_expressions(lexp,lname)
    G=g.composed_graph(verbose=False)
    # Gr=G.reverse()
    # tgr = [g for g in nx.topological_generations(Gr)]

    all_tsorts = [nx.topological_sort(G)]
    
    for i,ts in enumerate(all_tsorts):
        Gp = nx.DiGraph(G)
        at_eval       = nx.get_node_attributes(Gp,"eval")
        at_idx        = dict()
        local_mem     = dict()
        is_allocated  = dict()
        for node_id,v in enumerate(Gp.nodes):
            at_idx[v]=-1
            local_mem[node_id]=None
            is_allocated[node_id]=False
            # if Gp.out_degree(v) == 0:
            #     at_eval[v]=True

        nx.set_node_attributes(Gp,at_eval,"eval")
        nx.set_node_attributes(Gp,at_idx, "idx")
        W=list(ts)
        #print("ts = %d"%(i))
        while len(W)>0:
            v = W.pop()
            Gp, local_mem = visit_node(Gp, v, W, local_mem, is_allocated)
            at_eval    = nx.get_node_attributes(Gp,"eval")
            at_idx     = nx.get_node_attributes(Gp,"idx")
            if(not at_eval[v]):
                W.append(v)
            else:
                c_code  = ccode(v)
                if(at_idx[v]==-1):
                    store_node(v, at_idx, local_mem)
                    if(is_allocated[at_idx[v]] == True):
                        print("DENDRO_%d=%s;"%(at_idx[v],c_code))
                    else:
                        print("double DENDRO_%d=%s;"%(at_idx[v],c_code))
                        is_allocated[at_idx[v]]=True
                    
                    #print("if ( fabs((DENDRO_%d) - (%s))>1e-6) {printf(\"reduction error %s at DENDRO_%d=%%.8E expected=%%.8E \\n\",DENDRO_%d,%s);}"%(at_idx[v],c_code,v,at_idx[v],at_idx[v],c_code))

                else:
                    #print("if ( fabs((DENDRO_%d) - (%s))>1e-6) {printf(\"reduction error %s at DENDRO_%d=%%.8E expected=%%.8E \\n\",DENDRO_%d,%s);}"%(at_idx[v],c_code,v,at_idx[v],at_idx[v],c_code))
                    if v in lexp:
                        print("%s=DENDRO_%d;"%(lname[lexp.index(v)],at_idx[v]))
                        if(Gp.in_degree(v)==0):
                            evict_node(v,at_idx,local_mem)
                        
                
                
                        

            nx.set_node_attributes(Gp,at_eval,"eval")
            nx.set_node_attributes(Gp,at_idx, "idx")
        



