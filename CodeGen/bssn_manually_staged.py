'''
Try to stage BSSN computations manually. 
'''
from distutils.cygwinccompiler import CygwinCCompiler
import enum
from pyclbr import Function

from numpy import c_
import dendro
import contextlib 
import io
import importlib 
import nxgraph 
importlib.reload(nxgraph)
import networkx as nx
import sympy

from sympy import *

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

# Additional parameters for damping term
R0 = symbols('BSSN_ETA_R0')
ep1, ep2 = symbols('BSSN_ETA_POWER[0] BSSN_ETA_POWER[1]')

xi1, xi2, xi3 = symbols('BSSN_XI[0] BSSN_XI[1] BSSN_XI[2] ')
eta_damp=symbols("eta")
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

#d = dendro.set_first_derivative('grad')    # first argument is direction
#d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
#ad = dendro.set_advective_derivative('grad')  # first argument is direction
#kod = dendro.set_kreiss_oliger_dissipation('kograd')

# dendro.d   = lambda i,x : symbols("grad_%d_%s"%(i,x))
# dendro.d2  = lambda i,j,x : symbols("grad2_%d_%d_%s"%(min(i,j),max(i,j),x))
# dendro.d   = lambda i,x : Symbol("(deriv_evars->grad_%d_%s[d_pp])"%(i,str(x).split('[')[0]))
# dendro.d2  = lambda i,j,x : Symbol("deriv_evars->grad2_%d_%d_%s[d_pp]"%(min(i,j),max(i,j),str(x).split('[')[0]))
dendro.d   = lambda i,x : Symbol("grad_%d_%s"%(i,str(x).split('[')[0]))
dendro.d2  = lambda i,j,x : Symbol("grad2_%d_%d_%s"%(min(i,j),max(i,j),str(x).split('[')[0]))

dendro.ad  = dendro.d
dendro.kod = dendro.undef

d=dendro.d
d2=dendro.d2
ad=dendro.ad

dendro.set_metric(gt)
mi = [0, 1, 2, 4, 5, 8]
midx = ['00', '01', '02', '11', '12', '22']


def split_rhs_generation():
    MAX_TEMP_VARS=64
    CODEOUT_FOLDER="../BSSN_GR/cuda/scripts"
    temp_vars=["double DENDRO_%d;"%i for i in range(MAX_TEMP_VARS)]

    #point local variables. 
    igt  = dendro.sym_3x3("DENDRO_igt", "")
    R=dendro.sym_3x3("DENDRO_RIJ",'')
    CalGt = dendro.vec3("DENDRO_Gtk",'')

    C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

    C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

    C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

    _igt = dendro.get_inverse_metric()
    _igt = [_igt[i,j] for i in range(3) for j in range(i,3)]
    outs=[_igt]

    #store things on thread local variables. 
    precompute_code   = ""
    #declare_temp_vars = "".join(temp_vars)+"\n"
    declare_temp_vars = ""
    temp_vars  = [["DENDRO_igt%d"%i for i in range(6)],
                          ["DENDRO_C1_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_C2_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_C3_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_RIJ%d"%i for i in range(6)],
                          ["DENDRO_Gtk%d"%i for i in range(3)]
                         ] 

    for l in temp_vars:
        for v in l:
            declare_temp_vars +="double %s;\n"%v

    inverse_metric_code = ""
    christoffel_c1_code = ""
    christoffel_c2_code = ""
    Gt_from_metric_code = ""
    christoffel_c3_code = ""

    vnames=["DENDRO_igt"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
    inverse_metric_code= "\n{\n" +  f.getvalue() + "\n}\n"

    dendro.inv_metric=igt
    _C1 = dendro.get_first_christoffel()
    for k in dendro.e_i:
        _Ck = [_C1[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C1_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c1_code += "\n{\n" +  f.getvalue() + "\n}\n"
    dendro.C1 = C1

    _C2 = dendro.get_second_christoffel()
    for k in dendro.e_i:
        _Ck = [_C2[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C2_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c2_code += "\n{\n" +  f.getvalue() + "\n}\n"
    dendro.C2 = C2

    _C3 = dendro.get_complete_christoffel(chi)
    for k in dendro.e_i:
        _Ck = [_C3[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C3_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        
        christoffel_c3_code += "\n{\n" +  f.getvalue() + "\n}\n"
    dendro.C3 = C3

    [_R, _Rt, _Rphi, _CalGt] = dendro.compute_ricci(Gt, chi)
    _R = simplify(_R)
    _R = [_R[i,j] for i in range(3) for j in range(i,3)]

    outs=[_CalGt]
    vnames=["DENDRO_Gtk"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu(outs,vnames,'')
    Gt_from_metric_code  = "\n{\n" +  f.getvalue() + "\n}\n"

    # simplest variables, no connection coefficient / curvature tensor info needed. 
    a_rhs = l1*dendro.lie(b, a) - 2*a*K 
    b_rhs = [(Rational(3,4) * (lf0 + lf1*a) * B[i] + l2 * dendro.vec_j_ad_j(b, b[i])) for i in dendro.e_i ] 
    gt_rhs = dendro.lie(b, gt, weight) - 2*a*At 
    chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) 
    outs = [a_rhs, b_rhs, gt_rhs, chi_rhs]
    vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs']
    with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'[pp]')

    code_str =  "\n{\n" +  f.getvalue() + "\n}\n"
    precompute_code = ""
    with open(CODEOUT_FOLDER+"/a_beta_gt_chi.cpp",'w',encoding = 'utf-8') as code_file:
        code_file.write(precompute_code)
        code_file.write(code_str)    

    # computes each component seperately, 
    for i,m in enumerate(_R):
        outs=[m]
        vnames=["At_rhs%s"%midx[i]]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'[pp]')

        precompute_code = declare_temp_vars +  inverse_metric_code + christoffel_c1_code + christoffel_c2_code 
        c_code = f.getvalue()
        with open(CODEOUT_FOLDER+"/R_%s.cpp"%midx[i],'w',encoding = 'utf-8') as c_file:
            c_file.write(precompute_code)
            c_file.write(c_code)

    
    AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])
    At_rhs = chi*dendro.trace_free(a*R) - chi*dendro.trace_free(dendro.DiDj(a)) # + a*(K*At - 2*AikAkj.reshape(3, 3)) dendro.lie(b, At, weight) +
    load_ricci_tensor_code="".join(["DENDRO_RIJ%s=At_rhs%s[pp];\n"%(str(i),idx) for i,idx in enumerate(midx)])+"\n"
    At_rhs = [At_rhs[i,j] for i in range(3) for j in range(i,3)]
    precompute_code = declare_temp_vars +  inverse_metric_code + christoffel_c1_code + christoffel_c2_code + christoffel_c3_code + load_ricci_tensor_code
    with open(CODEOUT_FOLDER+"/At_rhs_tf.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(precompute_code)
        for i,m in enumerate(At_rhs):
            with contextlib.redirect_stdout(io.StringIO()) as f:
                dendro.generate_cpu([m],["At_rhs%s"%(midx[i])],'[pp]')
            c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
            c_file.write(c_code)

    precompute_code = declare_temp_vars +  inverse_metric_code
    At_rhs = dendro.lie(b, At, weight) + a*(K*At - 2*AikAkj.reshape(3, 3)) 
    At_rhs = [symbols("At_rhs%s%s[pp]"%(i,j)) + At_rhs[i,j] for i in range(3) for j in range(i,3)]
    with open(CODEOUT_FOLDER+"/At_rhs.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(precompute_code)
        for i,m in enumerate(At_rhs):
            with contextlib.redirect_stdout(io.StringIO()) as f:
                dendro.generate_cpu([m],["At_rhs%s"%(midx[i])],"[pp]")
            c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
            c_file.write(c_code)


    precompute_code = declare_temp_vars +  inverse_metric_code + christoffel_c1_code + christoffel_c2_code + christoffel_c3_code
    A_ijAuij = dendro.sqr(At)
    with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu_no_cse([expand(A_ijAuij)],["double AijAuij"],"")
    precompute_code+=f.getvalue()
    K_rhs = dendro.lie(b, K) + a*(K*K/3 + dendro.scalar("AijAuij",'')) 
    with open(CODEOUT_FOLDER+"/K_rhs1.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(precompute_code)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu([K_rhs],["K_rhs"],"[pp]")
        c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
        c_file.write(c_code)

    precompute_code = declare_temp_vars +  inverse_metric_code + christoffel_c1_code + christoffel_c2_code + christoffel_c3_code
    d_laplace=lambda a,chi : sum([ chi * dendro.inv_metric[i, j] * ( d2(i, j, a) - sum([C3[l, i, j] * d(l, a) for l in dendro.e_i]) ) for i, j in dendro.e_ij])
    K_rhs = dendro.scalar("K_rhs",'[pp]') - d_laplace(a,chi)
    with open(CODEOUT_FOLDER+"/K_rhs2.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(precompute_code)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu([K_rhs],["K_rhs"],"[pp]")
        c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
        c_file.write(c_code)




    At_UU = dendro.up_up(At)
    Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
            Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
            Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
            Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
            Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
            Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
            Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])
    Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]
    B_rhs = [ (dendro.vec3("Gt_rhs",'[pp]')[i] - eta_damp * B[i] +
            l3 * dendro.vec_j_ad_j(b, B[i]) -
            l4 * dendro.vec_j_ad_j(b, Gt[i]))
            for i in dendro.e_i ]


    precompute_code = declare_temp_vars +  inverse_metric_code + christoffel_c1_code + christoffel_c2_code + Gt_from_metric_code
    with open(CODEOUT_FOLDER+"/Gt_B_rhs.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(precompute_code)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu([Gt_rhs,B_rhs],["Gt_rhs","B_rhs"],"[pp]")
        c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
        c_file.write(c_code)


def rhs_reduced_temp_vars():
    MAX_TEMP_VARS=28
    CODEOUT_FOLDER="../BSSN_GR/cuda/scripts"
    temp_vars=["double DENDRO_%d;"%i for i in range(MAX_TEMP_VARS)]

    #point local variables. 
    igt  = dendro.sym_3x3("DENDRO_igt", "")
    R=dendro.sym_3x3("DENDRO_RIJ",'')
    CalGt = dendro.vec3("DENDRO_Gtk",'')

    C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

    C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

    C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

    _igt = dendro.get_inverse_metric()
    _igt = [_igt[i,j] for i in range(3) for j in range(i,3)]
    outs=[_igt]


    MAX_TEMP_VARS=64
    CODEOUT_FOLDER="../BSSN_GR/cuda/scripts"
    temp_vars=["double DENDRO_%d;"%i for i in range(MAX_TEMP_VARS)]

    #point local variables. 
    igt  = dendro.sym_3x3("DENDRO_igt", "")
    R=dendro.sym_3x3("DENDRO_RIJ",'')
    CalGt = dendro.vec3("DENDRO_Gtk",'')

    C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

    C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

    C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

    _igt = dendro.get_inverse_metric()
    _igt = [_igt[i,j] for i in range(3) for j in range(i,3)]
    outs=[_igt]

    #store things on thread local variables. 
    precompute_code   = ""
    declare_temp_vars = ""
    temp_vars  = [["DENDRO_igt%d"%i for i in range(6)],
                          ["DENDRO_C1_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_C2_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_C3_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_RIJ%d"%i for i in range(6)],
                          ["DENDRO_Gtk%d"%i for i in range(3)]
                         ] 
    inverse_metric_code = ""
    christoffel_c1_code = ""
    christoffel_c2_code = ""
    Gt_from_metric_code = ""
    christoffel_c3_code = ""
    ricci_code          = ""

    for l in temp_vars:
        for v in l:
            declare_temp_vars +="double %s;\n"%v

 
    vnames=["DENDRO_igt"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
    inverse_metric_code="{\n" + f.getvalue() + "\n}\n"
    
    dendro.inv_metric=igt
    _C1 = dendro.get_first_christoffel()
    for k in dendro.e_i:
        _Ck = [_C1[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C1_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c1_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C1 = C1

    _C2 = dendro.get_second_christoffel()
    for k in dendro.e_i:
        _Ck = [_C2[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C2_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c2_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C2 = C2

    _C3 = dendro.get_complete_christoffel(chi)
    for k in dendro.e_i:
        _Ck = [_C3[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C3_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        
        christoffel_c3_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C3 = C3

    [_R, _Rt, _Rphi, _CalGt] = dendro.compute_ricci(Gt, chi)
    _R = simplify(_R)
    _R = [_R[i,j] for i in range(3) for j in range(i,3)]

    for i,m in enumerate(_R):
        outs=[m]
        vnames=["DENDRO_RIJ%s"%str(i)]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')

        ricci_code += "{\n" + f.getvalue() + "\n}\n"

    outs=[_CalGt]
    vnames=["DENDRO_Gtk"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu(outs,vnames,'')
    Gt_from_metric_code  = "{\n" + f.getvalue() + "\n}\n"

    
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
    
    outs   = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
    vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']

    
    def component_wise(outs,vnames,precompute):
        print(precompute)
        for i, e in enumerate(outs):
            if type(e) == list:
                for j, ev in enumerate(e):
                    print("{")
                    dendro.generate_cpu([ev],[vnames[i]+str(j)],'[pp]')
                    print("}")
            elif type(e) == Matrix:
                for j, k in enumerate(mi):
                    print("{")
                    dendro.generate_cpu([e[k]],[vnames[i]+midx[j]],'[pp]')
                    print("}")
            else:
                print("{")
                dendro.generate_cpu([e],[vnames[i]],'[pp]')
                print("}")

    print(declare_temp_vars)     
    precompute_code=""
    component_wise([a_rhs, b_rhs, gt_rhs, chi_rhs],['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs'],precompute_code)

    precompute_code="" + inverse_metric_code+ christoffel_c1_code + christoffel_c2_code + christoffel_c3_code + ricci_code + Gt_from_metric_code
    component_wise([At_rhs, K_rhs, Gt_rhs, B_rhs],['At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs'],precompute_code)
    

def rhs_reduced_temp_vars_with_component_wise_parallelism():
    MAX_TEMP_VARS=28
    CODEOUT_FOLDER="../BSSN_GR/cuda/scripts"
    temp_vars=["double DENDRO_%d;"%i for i in range(MAX_TEMP_VARS)]

    #point local variables. 
    igt  = dendro.sym_3x3("DENDRO_igt", "")
    R=dendro.sym_3x3("DENDRO_RIJ",'')
    CalGt = dendro.vec3("DENDRO_Gtk",'')

    C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

    C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

    C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

    _igt = dendro.get_inverse_metric()
    _igt = [_igt[i,j] for i in range(3) for j in range(i,3)]
    outs=[_igt]


    MAX_TEMP_VARS=64
    CODEOUT_FOLDER="../BSSN_GR/cuda/scripts"
    temp_vars=["double DENDRO_%d;"%i for i in range(MAX_TEMP_VARS)]

    #point local variables. 
    igt  = dendro.sym_3x3("DENDRO_igt", "")
    R=dendro.sym_3x3("DENDRO_RIJ",'')
    CalGt = dendro.vec3("DENDRO_Gtk",'')

    C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

    C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

    C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

    _igt = dendro.get_inverse_metric()
    _igt = [_igt[i,j] for i in range(3) for j in range(i,3)]
    outs=[_igt]

    #store things on thread local variables. 
    precompute_code   = ""
    declare_temp_vars = ""
    temp_vars  = [["DENDRO_igt%d"%i for i in range(6)],
                          ["DENDRO_C1_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_C2_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_C3_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                          ["DENDRO_RIJ%d"%i for i in range(6)],
                          ["DENDRO_Gtk%d"%i for i in range(3)]
                         ] 
    inverse_metric_code = ""
    christoffel_c1_code = ""
    christoffel_c2_code = ""
    Gt_from_metric_code = ""
    christoffel_c3_code = ""
    ricci_code          = ""

    for l in temp_vars:
        for v in l:
            declare_temp_vars +="double %s;\n"%v

 
    vnames=["DENDRO_igt"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
    inverse_metric_code="{\n" + f.getvalue() + "\n}\n"
    
    dendro.inv_metric=igt
    _C1 = dendro.get_first_christoffel()
    for k in dendro.e_i:
        _Ck = [_C1[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C1_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c1_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C1 = C1

    _C2 = dendro.get_second_christoffel()
    for k in dendro.e_i:
        _Ck = [_C2[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C2_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c2_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C2 = C2

    _C3 = dendro.get_complete_christoffel(chi)
    for k in dendro.e_i:
        _Ck = [_C3[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C3_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        
        christoffel_c3_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C3 = C3

    [_R, _Rt, _Rphi, _CalGt] = dendro.compute_ricci(Gt, chi)
    _R = simplify(_R)
    _R = [_R[i,j] for i in range(3) for j in range(i,3)]

    for i,m in enumerate(_R):
        outs=[m]
        vnames=["DENDRO_RIJ%s"%str(i)]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')

        ricci_code += "{\n" + f.getvalue() + "\n}\n"

    outs=[_CalGt]
    vnames=["DENDRO_Gtk"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu(outs,vnames,'')
    Gt_from_metric_code  = "{\n" + f.getvalue() + "\n}\n"

    
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
    
    outs   = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
    vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']

    
    precompute_code = "" 
    with open(CODEOUT_FOLDER+"/block_id_0.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(precompute_code)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu([a_rhs, b_rhs, gt_rhs, chi_rhs], ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs'],"[pp]")
        c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
        c_file.write(c_code)

    precompute_code="" + inverse_metric_code+ christoffel_c1_code + christoffel_c2_code + christoffel_c3_code + ricci_code + Gt_from_metric_code
    with open(CODEOUT_FOLDER+"/block_id_1.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(declare_temp_vars)
        c_file.write(precompute_code)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu([At_rhs], ['At_rhs'],"[pp]")
        c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
        c_file.write(c_code)

    precompute_code="" + inverse_metric_code+ christoffel_c1_code + christoffel_c2_code + christoffel_c3_code + Gt_from_metric_code
    with open(CODEOUT_FOLDER+"/block_id_2.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(declare_temp_vars)
        c_file.write(precompute_code)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu([Gt_rhs,B_rhs],["Gt_rhs","B_rhs"],"[pp]")
        c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
        c_file.write(c_code)

    precompute_code="" + inverse_metric_code+ christoffel_c1_code + christoffel_c2_code + christoffel_c3_code + Gt_from_metric_code
    with open(CODEOUT_FOLDER+"/block_id_3.cpp",'w',encoding = 'utf-8') as c_file:
        c_file.write(declare_temp_vars)
        c_file.write(precompute_code)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu([K_rhs],["K_rhs"],"[pp]")
        c_code =  "\n{\n" +  f.getvalue() + "\n}\n"
        c_file.write(c_code)
    

def rhs_with_d_computations():

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


    #point local variables. 
    igt  = dendro.sym_3x3("DENDRO_igt", "")
    R=dendro.sym_3x3("DENDRO_RIJ",'')
    CalGt = dendro.vec3("DENDRO_Gtk",'')

    C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

    C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

    C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

    #point local variables. 
    igt  = dendro.sym_3x3("DENDRO_igt", "")
    R=dendro.sym_3x3("DENDRO_RIJ",'')
    CalGt = dendro.vec3("DENDRO_Gtk",'')

    C1 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C1[k,:,:]=dendro.sym_3x3("DENDRO_C1_k%d_"%k,'')

    C2 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C2[k,:,:]=dendro.sym_3x3("DENDRO_C2_k%d_"%k,'')

    C3 = MutableDenseNDimArray(range(27), (3, 3, 3))
    for k in dendro.e_i:
        C3[k,:,:]=dendro.sym_3x3("DENDRO_C3_k%d_"%k,'')

    _igt = dendro.get_inverse_metric()
    _igt = [_igt[i,j] for i in range(3) for j in range(i,3)]
    outs=[_igt]

    #store things on thread local variables. 
    precompute_code   = ""
    declare_temp_vars = ""
    temp_vars  = [["DENDRO_igt%d"%i for i in range(6)],
                            ["DENDRO_C1_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                            ["DENDRO_C2_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                            ["DENDRO_C3_k%d_%d"%(i,j) for i in range(3) for j in range(6)],
                            ["DENDRO_RIJ%d"%i for i in range(6)],
                            ["DENDRO_Gtk%d"%i for i in range(3)]
                            ] 
    inverse_metric_code = ""
    christoffel_c1_code = ""
    christoffel_c2_code = ""
    Gt_from_metric_code = ""
    christoffel_c3_code = ""
    ricci_code          = ""

    for l in temp_vars:
        for v in l:
            declare_temp_vars +="double %s;\n"%v


    vnames=["DENDRO_igt"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
    inverse_metric_code="{\n" + f.getvalue() + "\n}\n"

    dendro.inv_metric=igt
    _C1 = dendro.get_first_christoffel()
    for k in dendro.e_i:
        _Ck = [_C1[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C1_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c1_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C1 = C1

    _C2 = dendro.get_second_christoffel()
    for k in dendro.e_i:
        _Ck = [_C2[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C2_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        christoffel_c2_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C2 = C2

    _C3 = dendro.get_complete_christoffel(chi)
    for k in dendro.e_i:
        _Ck = [_C3[k,i,j] for i in range(3) for j in range(i,3)]
        outs=[_Ck]
        vnames=["DENDRO_C3_k%d_"%k]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')
        
        christoffel_c3_code += "{\n" + f.getvalue() + "\n}\n"
    dendro.C3 = C3

    [_R, _Rt, _Rphi, _CalGt] = dendro.compute_ricci(Gt, chi)
    _R = simplify(_R)
    _R = [_R[i,j] for i in range(3) for j in range(i,3)]

    for i,m in enumerate(_R):
        outs=[m]
        vnames=["DENDRO_RIJ%s"%str(i)]
        with contextlib.redirect_stdout(io.StringIO()) as f:
            dendro.generate_cpu(outs,vnames,'')

        ricci_code += "{\n" + f.getvalue() + "\n}\n"

    outs=[_CalGt]
    vnames=["DENDRO_Gtk"]
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu(outs,vnames,'')
    Gt_from_metric_code  = "{\n" + f.getvalue() + "\n}\n"


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

    
    a_derivs   = first_and_second_derivs("alpha",True)  + ko_derivs_only("alpha",False)
    b_derivs   = "\n".join([first_and_second_derivs("beta%d"%i,True)     + ko_derivs_only("beta%d"%i,False) for i in range(3) ])
    
    c_code     = ""
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu([a_rhs],['a_rhs'],'[pp]')
    c_code += f.getvalue()
    c_code +=apply_bc("a_rhs")
    c_code +=apply_ko("a_rhs")
    print("%s\n{\n%s\n}\n"%(a_derivs,c_code))

    c_code     = ""
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu([b_rhs],['b_rhs'],'[pp]')
    c_code += f.getvalue()
    c_code +=apply_bc("b_rhs")
    c_code +=apply_ko("b_rhs")
    print("%s\n{\n%s\n}\n"%(b_derivs, c_code))

    chi_derivs = first_and_second_derivs("chi",True)  + ko_derivs_only("chi",False)
    c_code     = ""
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu([chi_rhs],['chi_rhs'],'[pp]')
    c_code += f.getvalue()
    c_code +=apply_bc("chi_rhs")
    c_code +=apply_ko("chi_rhs")
    print("%s\n{\n%s\n}\n"%(chi_derivs, c_code))

    gt_derivs  = "\n".join([first_and_second_derivs("gt%d"%i,True) + ko_derivs_only("gt%d"%i,False) for i in range(6)   ])
    c_code     = ""
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu([gt_rhs],['gt_rhs'],'[pp]')
    c_code += f.getvalue()
    c_code +=apply_bc("gt_rhs")
    c_code +=apply_ko("gt_rhs")
    print("%s\n{\n%s\n}\n"%(gt_derivs, c_code))

    # print(declare_temp_vars)
    # print(inverse_metric_code + christoffel_c1_code + christoffel_c2_code + christoffel_c3_code + Gt_from_metric_code)
    K_derivs   = first_derivs_only("K",True)  + ko_derivs_only("K",False)
    c_code     = ""
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu([K_rhs],['K_rhs'],'[pp]')
    c_code +=declare_temp_vars
    c_code +=inverse_metric_code + christoffel_c1_code + christoffel_c2_code + christoffel_c3_code 
    c_code += f.getvalue()
    c_code +=apply_bc("K_rhs")
    c_code +=apply_ko("K_rhs")
    print("%s\n{\n%s\n}\n"%(K_derivs, c_code))

    Gt_derivs  = "\n".join([first_derivs_only("Gt%d"%i,True)    + ko_derivs_only("Gt%d"%i,False) for i in range(3) ])
    B_derivs   = "\n".join([first_derivs_only("B%d"%i,True)     + ko_derivs_only("B%d"%i,False) for i in range(3) ])
    
    c_code     = ""
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu([Gt_rhs,B_rhs],['Gt_rhs', 'B_rhs'],'[pp]')
    c_code +=declare_temp_vars
    c_code +=inverse_metric_code + christoffel_c1_code + christoffel_c2_code + christoffel_c3_code + Gt_from_metric_code
    c_code += f.getvalue()
    c_code +=apply_bc("Gt_rhs")
    c_code +=apply_ko("Gt_rhs")
    c_code +=apply_bc("B_rhs")
    c_code +=apply_ko("B_rhs")
    print("%s\n{\n%s\n}\n"%(Gt_derivs + B_derivs, c_code))

    At_derivs  = "\n".join([first_derivs_only("At%d"%i,True) + ko_derivs_only("At%d"%i,False) for i in range(6)])
    c_code     = ""
    with contextlib.redirect_stdout(io.StringIO()) as f:
        dendro.generate_cpu([At_rhs],['At_rhs'],'[pp]')
    c_code +=declare_temp_vars
    c_code +=inverse_metric_code + christoffel_c1_code + christoffel_c2_code + christoffel_c3_code  + ricci_code
    c_code += f.getvalue()
    c_code +=apply_bc("At_rhs")
    c_code +=apply_ko("At_rhs")
    print("%s\n{\n%s\n}\n"%(At_derivs, c_code))



#dendro.generate_cpu(outs,vnames, '[pp]')
#dendro.generate_cpu([At_rhs],['At_rhs'], '[pp]')
# outs   = [a_rhs, b_rhs]
# vnames = ['a_rhs', 'b_rhs']
#exp_graph.add_expressions(outs,vnames)
#G=exp_graph.composed_graph(verbose=True)


#split_rhs_generation()
#rhs_reduced_temp_vars()
#rhs_reduced_temp_vars_with_component_wise_parallelism()
#traversal_base_code_gen()
rhs_with_d_computations()