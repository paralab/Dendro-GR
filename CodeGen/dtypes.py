##########################################################################
# module: symbolic
# author: Hari Sundar
# email:  hari@cs.utah.edu
# 
# author: Milinda fernando
# email: milinda@cs.utah.edu
# (c) 2019 University of Utah, All rights reserved.
'''
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense,and/or sell copies
of the Software,and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software. 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.
'''
##########################################################################

from sympy import *
from sympy.tensor.array import *
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.utilities import numbered_symbols
from sympy.printing import print_ccode
from sympy.printing.dot import dotprint

import re as regex

import string
import random

# internal variables
undef = symbols('undefined')


##########################################################################
# variable initialization functions
##########################################################################

def scalar(name):
    """
    Create a scalar variable with the corresponding name. The 'name' will be during code generation, so should match the
    variable name used in the C++ code.
    """
    tname = name
    return symbols(tname)


def vec(name,n):
    """
    Create a nD vector variable with the corresponding name. The 'name' will be during code generation, so should match
    the variable name used in the C++ code. The returned variable can be indexed(0,1,2), i.e.,

    b = dendro.vec("beta")
    b[1] = x^2
    """

    vname=list()
    
    for i in range(0,n):
        nameR = name + repr(i) 
        vname.append(nameR)
    
    return Matrix([symbols(vname[i]) for i in range(0,n)])


def mat(name,n,m):
    
    """
    Creates a symbolic matrix of size nxm
    """
    vname=list()
    
    for i in range(0,n):
        nameR = name + repr(i) 
        nameC = ' '.join([nameR + repr(j)  for j in range(0,m)])
        vname.append(nameC)
    
    return Matrix([symbols(vname[i]) for i in range(0,n)])
    

def matI(n):
    
    """
    Creates a Idendity matrix of size nxn
    """
    
    return eye(n)
    



