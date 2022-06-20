"""
bssn code generation with network x
"""
import random
from numpy import Infinity
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
import queue
from collections import OrderedDict
import pickle
import argparse

DX  = "device::__blk1_deriv644_x"
DY  = "device::__blk1_deriv644_y"
DZ  = "device::__blk1_deriv644_z"

DXX = "device::__blk1_deriv644_xx"
DYY = "device::__blk1_deriv644_yy"
DZZ = "device::__blk1_deriv644_zz"

KDX = "device::__blk1_ko_deriv42_x"
KDY = "device::__blk1_ko_deriv42_y"
KDZ = "device::__blk1_ko_deriv42_z"

LD   = "device::__ld_blk_var1__"
SD   = "device::__st_blk_var1__"
SV   = "su"
DSV  = "Du"
DDSV = "DDu"
ST  = "__syncthreads();"

D = ["alpha", "chi", "K", "Gt0", "Gt1", "Gt2", "beta0", "beta1", "beta2",
    "B0", "B1", "B2", "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
    "At0", "At1", "At2", "At3", "At4", "At5" ]

# second derivs required for RHS
DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi",
    "alpha", "beta0", "beta1", "beta2" ]

# Kries-Oliger derivatives. 
# KO = [ "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
#        "At0", "At1", "At2", "At3", "At4", "At5",
#        "alpha", "beta0", "beta1", "beta2",
#        "chi", "Gt0", "Gt1", "Gt2", "K",
#        "B0", "B1", "B2"]

# first derivs required for constraints--no gauge variables
CONSTRAINT_D = [ "chi", "Gt0", "Gt1", "Gt2", "K",
                "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
                "At0", "At1", "At2", "At3", "At4", "At5" ]

# second derivs required for constraints--no gauge variables
CONSTRAINT_DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi"]

# custom functions for code generation in cse.
custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}


## NOTE: Prefix ordering matters!!! during the deriv allocation code generation. 
PREFIX_D   = ["grad_0_",    "grad_1_",      "grad_2_"]
PREFIX_DD  = ["grad2_0_0_", "grad2_0_1_",   "grad2_0_2_", "grad2_1_1_", "grad2_1_2_", "grad2_2_2_"]
#PREFIX_AD  = ["agrad_0_",   "agrad_1_",     "agrad_2_"]
PREFIX_KOD = ["kograd_0_",  "kograd_1_",    "kograd_2_"]


# first derivative in i direction
FUNC_D_I=[]
for f in D:
    for p in PREFIX_D:
        FUNC_D_I.append(p+f)

# second derivative in ij direction
FUNC_D_IJ=[]
for f in DD:
    for p in PREFIX_DD:
        FUNC_D_IJ.append(p+f)

FUNC_KOD_I=[]
for f in D:
    for p in PREFIX_KOD:
        FUNC_KOD_I.append(p+f)

FUNC_CONS=[]
for f in CONSTRAINT_D:
    for p in PREFIX_D:
        FUNC_CONS.append(p+f)
        
for f in CONSTRAINT_DD:
    for p in PREFIX_DD:
        FUNC_CONS.append(p+f)



@dataclass(order=True)
class PrioritizedItem:
    priority: int
    item: Any=field(compare=False)

def node_priority(node, G, at_eval):
    is_ready=True
    for u in G.successors(node):
        if at_eval[u]==False:
            is_ready=False
            break
    
    if is_ready:
        priority=0
    else:
        priority=Infinity
    
    return priority


def first_derivs_only(var, load_to_shared=True):
    c_code=""
    if(load_to_shared):
        c_code+=("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV, var))
    
    # x
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[0] +  var, DSV))
    c_code+=("%s\n"%ST)

    # y
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[1] +  var, DSV))
    c_code+=("%s\n"%ST)

    # z
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[2] +  var, DSV))
    c_code+=("%s\n"%ST)

    return c_code
        
def ko_derivs_only(var, load_to_shared=True):
    c_code=""
    if(load_to_shared):
        c_code+=("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV, var))
    
    # kox
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDX, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[0] +  var, DSV))
    c_code+=("%s\n"%ST)

    # koy
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDY, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[1] +  var, DSV))
    c_code+=("%s\n"%ST)

    # koz
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDZ, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[2] +  var, DSV))
    c_code+=("%s\n"%ST)

    return c_code
    
def first_and_second_derivs(var, load_to_shared=True):
    c_code=""
    # load
    if(load_to_shared):
        c_code+=("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV, var))
    #xx
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DXX, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[0] +  var, DSV))
    c_code+=("%s\n"%ST)

    
    #x
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[0] +  var, DSV))
    c_code+=("%s\n"%ST)

    #xy
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, DDSV, DSV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[1] +  var, DDSV))
    c_code+=("%s\n"%ST)

    #xz
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DDSV, DSV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[2] +  var, DDSV))
    c_code+=("%s\n"%ST)

    #yy
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DYY, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[3] +  var, DSV))
    c_code+=("%s\n"%ST)

    # y 
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[1] +  var, DSV))
    c_code+=("%s\n"%ST)

    # #yz
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DDSV, DSV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[4] +  var, DDSV))
    c_code+=("%s\n"%ST)

    #zz
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZZ, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[5] +  var, DSV))
    c_code+=("%s\n"%ST)

    #z
    c_code+=("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DSV, SV))
    c_code+=("%s\n"%ST)
    c_code+=("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[2] +  var, DSV))
    c_code+=("%s\n"%ST)

    return c_code

def apply_bc(var):
    bc_code=""
    if(var=="a_rhs"):
        bc_code = "radiative_bc_pt<pw , nx>(&a_rhs [gidx] , alpha[gidx], grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk);\n"
    elif (var == "b_rhs"):
        bc_code+= "radiative_bc_pt<pw , nx>(&b_rhs0[gidx] , beta0[gidx], grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk);\n"
        bc_code+= "radiative_bc_pt<pw , nx>(&b_rhs1[gidx] , beta1[gidx], grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk);\n"
        bc_code+= "radiative_bc_pt<pw , nx>(&b_rhs2[gidx] , beta2[gidx], grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk);\n"
    elif (var == "chi_rhs"):
        bc_code = "radiative_bc_pt<pw , nx>(&chi_rhs[gidx], chi[gidx]  , grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk);\n"
    elif (var == "K_rhs"):
        bc_code = "radiative_bc_pt<pw , nx>(&K_rhs[gidx],     K[gidx]  ,   grad_0_K,   grad_1_K,   grad_2_K, 1.0, 0.0, blk);\n"
    elif (var == "Gt_rhs"):
        bc_code+="radiative_bc_pt<pw , nx>(&Gt_rhs0[gidx], Gt0[gidx], grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&Gt_rhs1[gidx], Gt1[gidx], grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&Gt_rhs2[gidx], Gt2[gidx], grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk);\n"
    elif (var == "B_rhs"):
        bc_code+="radiative_bc_pt<pw , nx>(&B_rhs0[gidx], B0[gidx], grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&B_rhs1[gidx], B1[gidx], grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&B_rhs2[gidx], B2[gidx], grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk);\n"
    elif (var == "At_rhs"):
        bc_code+="radiative_bc_pt<pw , nx>(&At_rhs00[gidx], At0[gidx], grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&At_rhs01[gidx], At1[gidx], grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&At_rhs02[gidx], At2[gidx], grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&At_rhs11[gidx], At3[gidx], grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&At_rhs12[gidx], At4[gidx], grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&At_rhs22[gidx], At5[gidx], grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk);\n"
    elif (var == "gt_rhs"):
        bc_code+="radiative_bc_pt<pw , nx>(&gt_rhs00[gidx], gt0[gidx], grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&gt_rhs01[gidx], gt1[gidx], grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&gt_rhs02[gidx], gt2[gidx], grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&gt_rhs11[gidx], gt3[gidx], grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&gt_rhs12[gidx], gt4[gidx], grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk);\n"
        bc_code+="radiative_bc_pt<pw , nx>(&gt_rhs22[gidx], gt5[gidx], grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk);\n"

    return "if(blk->m_bflag != 0){\n%s}\n"%(bc_code)
    
def apply_ko(var):
    ko_code=""
    if(var == "a_rhs"):
        ko_code =  "a_rhs[pp]    += ko_sigma * (kograd_0_alpha + kograd_1_alpha + kograd_2_alpha);\n"
    elif(var == "b_rhs"):
        ko_code += "b_rhs0[pp]   += ko_sigma * (kograd_0_beta0 + kograd_1_beta0 + kograd_2_beta0);\n"
        ko_code += "b_rhs1[pp]   += ko_sigma * (kograd_0_beta1 + kograd_1_beta1 + kograd_2_beta1);\n"
        ko_code += "b_rhs2[pp]   += ko_sigma * (kograd_0_beta2 + kograd_1_beta2 + kograd_2_beta2);\n"
    elif(var == "gt_rhs"):
        ko_code += "gt_rhs00[pp] += ko_sigma * (kograd_0_gt0 + kograd_1_gt0 + kograd_2_gt0);\n"
        ko_code += "gt_rhs01[pp] += ko_sigma * (kograd_0_gt1 + kograd_1_gt1 + kograd_2_gt1);\n"
        ko_code += "gt_rhs02[pp] += ko_sigma * (kograd_0_gt2 + kograd_1_gt2 + kograd_2_gt2);\n"
        ko_code += "gt_rhs11[pp] += ko_sigma * (kograd_0_gt3 + kograd_1_gt3 + kograd_2_gt3);\n"
        ko_code += "gt_rhs12[pp] += ko_sigma * (kograd_0_gt4 + kograd_1_gt4 + kograd_2_gt4);\n"
        ko_code += "gt_rhs22[pp] += ko_sigma * (kograd_0_gt5 + kograd_1_gt5 + kograd_2_gt5);\n"
    elif(var == "chi_rhs"):
        ko_code =  "chi_rhs[pp]  += ko_sigma * (kograd_0_chi + kograd_1_chi + kograd_2_chi);\n"
    elif(var == "At_rhs"):
        ko_code +="At_rhs00[pp] += ko_sigma * (kograd_0_At0 + kograd_1_At0 + kograd_2_At0);\n"
        ko_code +="At_rhs01[pp] += ko_sigma * (kograd_0_At1 + kograd_1_At1 + kograd_2_At1);\n"
        ko_code +="At_rhs02[pp] += ko_sigma * (kograd_0_At2 + kograd_1_At2 + kograd_2_At2);\n"
        ko_code +="At_rhs11[pp] += ko_sigma * (kograd_0_At3 + kograd_1_At3 + kograd_2_At3);\n"
        ko_code +="At_rhs12[pp] += ko_sigma * (kograd_0_At4 + kograd_1_At4 + kograd_2_At4);\n"
        ko_code +="At_rhs22[pp] += ko_sigma * (kograd_0_At5 + kograd_1_At5 + kograd_2_At5);\n"
    elif(var =="K_rhs"):
        ko_code +="K_rhs[pp]    += ko_sigma * (kograd_0_K + kograd_1_K + kograd_2_K);\n"
    elif(var =="Gt_rhs"):
        ko_code +="Gt_rhs0[pp]  += ko_sigma * (kograd_0_Gt0 + kograd_1_Gt0 + kograd_2_Gt0);\n"
        ko_code +="Gt_rhs1[pp]  += ko_sigma * (kograd_0_Gt1 + kograd_1_Gt1 + kograd_2_Gt1);\n"
        ko_code +="Gt_rhs2[pp]  += ko_sigma * (kograd_0_Gt2 + kograd_1_Gt2 + kograd_2_Gt2);\n"
    elif(var == "B_rhs"):
        ko_code +="B_rhs0[pp]   += ko_sigma * (kograd_0_B0 + kograd_1_B0 + kograd_2_B0);\n"
        ko_code +="B_rhs1[pp]   += ko_sigma * (kograd_0_B1 + kograd_1_B1 + kograd_2_B1);\n"
        ko_code +="B_rhs2[pp]   += ko_sigma * (kograd_0_B2 + kograd_1_B2 + kograd_2_B2);\n"

    return ko_code + ST + "\n"

def store_node(v, at_idx, local_mem):
    # for key,val in local_mem.items():
    #     if val is None:
    #         local_mem[key]=v
    #         print("storing node %s at %d"%(v,key))
    #         break

    # key=list(local_mem.keys())[list(local_mem.values()).index(None)]
    # local_mem[key]=v
    if local_mem["MAX_TEMP_VARS"] <= at_idx:
        local_mem["MAX_TEMP_VARS"]=at_idx

    local_mem[v] = at_idx
    at_idx+=1
    return at_idx
    
def evict_node(v,at_idx,local_mem):
    local_mem[v]=-1
    at_idx+=-1
    return at_idx

def visit_node(G: nx.DiGraph, v, work_queue, local_mem, at_idx):
    at_eval    = nx.get_node_attributes(G,"eval")
    at_func    = nx.get_node_attributes(G,"func") 
    at_args    = nx.get_node_attributes(G,"args") 

    descendents_list = list(G.successors(v))
    if(at_func[v] != sympy.core.add.Add  and at_func[v]!=sympy.core.mul.Mul):
        """
        Direct evaluation for non-reduction type functions such as pow
        """
        #print(at_func[v],v)
        if(at_eval[v]==False):
            if(at_func[v] == sympy.core.power.Pow):
                a1=hash(at_args[v][0])
                a2=hash(at_args[v][1])

                # print(a1,a2)
                # print("a1 is at DENDRO_%d\n"%at_idx[a1])
                # print(type(a2))
                if(local_mem[a1]==-1):
                    return at_idx #assert False, "invalid traversal order"

                if(type(a2) == sympy.core.numbers.Integer or type(a2) == sympy.core.numbers.One or type(a2) == sympy.core.numbers.NegativeOne):
                    if(local_mem[v]==-1):
                        at_idx=store_node(v,at_idx,local_mem)
                        print("DENDRO_%d = DENDRO_%d;\n"%(local_mem[v],local_mem[a1])) 

                    for i in range(abs(int(a2))-1):
                        print("DENDRO_%d *= DENDRO_%d;"%(local_mem[v],local_mem[a1])) 
                        
                    if(int(a2)<0):
                        print("DENDRO_%d = 1/DENDRO_%d;"%(local_mem[v],local_mem[v]))
                 
                else:
                    c_code = ccode(v)
                    if(local_mem[v]==-1):
                        at_idx=store_node(v,at_idx,local_mem)
                        print("DENDRO_%d = %s;"%(local_mem[v],c_code)) 
                            
                G.remove_edge(v,a1)
                G.remove_edge(v,a2)

                if(G.in_degree(a1)==0):
                    at_idx=evict_node(a1,at_idx,local_mem)

                if(G.in_degree(a2)==0):
                    at_idx=evict_node(a2,at_idx,local_mem)

                at_eval[v]=True
        
            else:
                c_code = ccode(v)
                if(local_mem[v]==-1):
                    at_idx=store_node(v,at_idx,local_mem)
                    print("DENDRO_%d = %s;"%(local_mem[v],c_code)) 
                    
    else:
        at_eval[v] = True
        for u in descendents_list:
            if(at_eval[u] == True):
                if(local_mem[u]==-1):
                    c_code = ccode(u)
                    at_idx=store_node(u,at_idx,local_mem)
                    print("DENDRO_%d=%s;"%(local_mem[u],c_code))
                    
                    
                if(local_mem[v]==-1):
                    at_idx=store_node(v,at_idx,local_mem)
                    #print("\n// initialize reduce for %s"%v)
                    #print("\n// initialize reduction for ")
                    if(at_func[v] == sympy.core.add.Add):
                        print("DENDRO_%d = 0;\n"%(local_mem[v])) 
                    elif(at_func[v] == sympy.core.mul.Mul):
                        print("DENDRO_%d = 1;\n"%(local_mem[v])) 
                            
                        
                if(at_func[v] == sympy.core.add.Add):
                    print("DENDRO_%d += DENDRO_%d;"%(local_mem[v],local_mem[u])) 
                elif(at_func[v] == sympy.core.mul.Mul):
                    print("DENDRO_%d *= DENDRO_%d;"%(local_mem[v],local_mem[u])) 
                
                G.remove_edge(v,u)
                if G.in_degree(u)==0:
                    at_idx=evict_node(u,at_idx,local_mem)

            else:
                at_eval[v]=False
    
    nx.set_node_attributes(G, at_eval, "eval")
    return at_idx

def generate_code_nx_pq(ex, vnames, idx, iter=1000):
    """
    Generate the C++ code 
    random traversal of the network x graph, with priority queue. 
    in the beginging leaf nodes has the highest priority (lowest value top of the queue)
    """
    # print(ex)
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    
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
    print("|V|= %d |E| =%d"%(G.number_of_nodes(), G.number_of_edges()))
    optimal_traversal_cost=Infinity
    optimal_ts=list()
    
    ts=list(G.nodes())
    for s_id in range(iter):
        random.shuffle(ts)
        #print(ts)
        Gp = G.copy()
        at_eval       = nx.get_node_attributes(Gp,"eval")
        at_idx        = 0
        local_mem     = dict()
        local_mem["MAX_TEMP_VARS"]=0
        for node_id,v in enumerate(Gp.nodes):
            local_mem[v]=-1
            if Gp.out_degree(v) == 0:
                at_eval[v]=True

        nx.set_node_attributes(Gp,at_eval,"eval")
        
        W=queue.PriorityQueue()
        for v in ts:
            W.put(PrioritizedItem(node_priority(v,G, at_eval), v))
            
        with contextlib.redirect_stdout(io.StringIO()) as f:
            while len(W.queue) > 0:
                #v = W.pop()
                p_item = W.get()
                #print(p_item)
                v= p_item.item
                at_idx     = visit_node(Gp, v, W, local_mem,at_idx)
                at_eval    = nx.get_node_attributes(Gp,"eval")
                if(not at_eval[v]):
                    for i in range(len(W.queue)):
                        p_item=W.get()
                        p_item.priority = node_priority(p_item.item, G, at_eval)
                        W.put(p_item)
                    
                else:
                    c_code  = ccode(v)
                    if(local_mem[v]==-1):
                        at_idx=store_node(v, at_idx, local_mem)
                        print("DENDRO_%d=%s;"%(local_mem[v],c_code))
                        #print("if ( fabs((DENDRO_%d) - (%s))>1e-6) {printf(\"reduction error %s at DENDRO_%d=%%.8E expected=%%.8E \\n\",DENDRO_%d,%s);}"%(local_mem[v],c_code,v,local_mem[v],local_mem[v],c_code))

                    else:
                        #print("if ( fabs((DENDRO_%d) - (%s))>1e-6) {printf(\"reduction error %s at DENDRO_%d=%%.8E expected=%%.8E \\n\",DENDRO_%d,%s);}"%(local_mem[v],c_code,v,local_mem[v],local_mem[v],c_code))
                        if v in lexp:
                            print("%s=DENDRO_%d;"%(lname[lexp.index(v)],local_mem[v]))
                    
                    if(Gp.in_degree(v)==0):
                        at_idx=evict_node(v,at_idx,local_mem)

                nx.set_node_attributes(Gp,at_eval,"eval")
                
                for i in range(len(W.queue)):
                    p_item=W.get()
                    p_item.priority = node_priority(p_item.item, G, at_eval)
                    W.put(p_item)
                
        print("traversal id=%d temp_var_count=%d"%(s_id,local_mem["MAX_TEMP_VARS"]+1))
        if optimal_traversal_cost > local_mem["MAX_TEMP_VARS"]+1 : 
            optimal_traversal_cost=local_mem["MAX_TEMP_VARS"]+1
            optimal_ts=list(ts)
        
    return [optimal_ts,optimal_traversal_cost,G]

def generate_code_nx(ex, vnames, idx, iter=1000):
    """
    Generate the C++ code
    network x traversal of the graph, and 
    rely on the topological order to evaluate. 
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
    print("|V|= %d |E| =%d"%(G.number_of_nodes(), G.number_of_edges()))
    
    optimal_traversal_cost=Infinity
    optimal_ts=list()
    import numpy as np
    #for i , q in enumerate(all_tsorts):
    tgr=list(nx.topological_generations(G))
    print("topological generations = %d"%(len(tgr)))
    
    for s_id in range(iter):
        ts=list()
        for tg in tgr:
            random.shuffle(tg)
            ts+=tg
        
        #print(ts)
    # all_tsorts = nx.all_topological_sorts(G)
    # for s_id , q in enumerate(all_tsorts):
    #     ts=list(q)
        Gp = G.copy()
        at_eval       = nx.get_node_attributes(Gp,"eval")
        at_idx        = 0
        local_mem     = dict()
        local_mem["MAX_TEMP_VARS"]=0
        for node_id,v in enumerate(Gp.nodes):
            local_mem[v]=-1
            if Gp.out_degree(v) == 0:
                at_eval[v]=True

        nx.set_node_attributes(Gp,at_eval,"eval")
        W=list(ts)
        
        with contextlib.redirect_stdout(io.StringIO()) as f:
            while len(W)>0:
                v = W.pop()
                at_idx     = visit_node(Gp, v, W, local_mem, at_idx)
                at_eval    = nx.get_node_attributes(Gp,"eval")
                if(not at_eval[v]):
                    print(v, "is not ready invalid traversal order for the Work stack")
                    assert False
                    
                else:
                    c_code  = ccode(v)
                    if(local_mem[v]==-1):
                        at_idx=store_node(v, at_idx, local_mem)
                        print("DENDRO_%d=%s;"%(local_mem[v],c_code))
                        
                        #print("if ( fabs((DENDRO_%d) - (%s))>1e-6) {printf(\"reduction error %s at DENDRO_%d=%%.8E expected=%%.8E \\n\",DENDRO_%d,%s);}"%(local_mem[v],c_code,v,local_mem[v],local_mem[v],c_code))

                    else:
                        #print("if ( fabs((DENDRO_%d) - (%s))>1e-6) {printf(\"reduction error %s at DENDRO_%d=%%.8E expected=%%.8E \\n\",DENDRO_%d,%s);}"%(local_mem[v],c_code,v,local_mem[v],local_mem[v],c_code))
                        if v in lexp:
                            print("%s=DENDRO_%d;"%(lname[lexp.index(v)],local_mem[v]))
                    
                    if(Gp.in_degree(v)==0):
                        at_idx=evict_node(v,at_idx,local_mem)

                nx.set_node_attributes(Gp,at_eval,"eval")
            
        print("traversal id=%d temp_var_count=%d"%(s_id,local_mem["MAX_TEMP_VARS"]+1))
        if optimal_traversal_cost > local_mem["MAX_TEMP_VARS"]+1 : 
            optimal_traversal_cost=local_mem["MAX_TEMP_VARS"]+1
            optimal_ts=list(ts)

    return [optimal_ts,optimal_traversal_cost,G]
    
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

# igt  = dendro.sym_3x3("DENDRO_igt", "")
# R=dendro.sym_3x3("DENDRO_RIJ",'')
# CalGt = dendro.vec3("DENDRO_Gtk",'')

# C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
# for k in dendro.e_i:
#     C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

# C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
# for k in dendro.e_i:
#     C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

# C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
# for k in dendro.e_i:
#     C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

# dendro.inv_metric=igt
# dendro.C1=C1
# dendro.C2=C2
# dendro.C3=C3
# dendro.Ricci=R

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

B_rhs = [ (Gt_rhs[i] - eta_damp * B[i] +
        l3 * dendro.vec_j_ad_j(b, B[i]) -
        l4 * dendro.vec_j_ad_j(b, Gt[i]))
        for i in dendro.e_i ]

outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']

parser = argparse.ArgumentParser(description='minimize the local variable usage for a given computational graph')
parser.add_argument("-g", "--graph_name", help="graph output", type=str, default="DAG_BSSN.dat")
parser.add_argument("-ts", "--opt_ts_name", help="output the optimal traversal found", type=str, default="ts_opt.dat")
parser.add_argument("-m",  "--mode", help="output the optimal traversal found (0- tsort based or 1 - pq based)", type=int, default=0)
parser.add_argument("-i",  "--iteration", help="output the optimal traversal found", type=int, default="100")
args = parser.parse_args()
print(args)
if args.mode==0:
    [ts,ts_cost,G]=generate_code_nx(outs,vnames,"[pp]",args.iteration)
elif args.mode==1:
    [ts,ts_cost,G]=generate_code_nx_pq(outs,vnames,"[pp]",args.iteration)

print(ts_cost)
with open(args.opt_ts_name, "wb") as fp:   #Pickling
    pickle.dump(ts, fp)
nx.write_gpickle(G,args.graph_name)
