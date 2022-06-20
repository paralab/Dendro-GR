'''
 BSSN core variables . 
'''
import sys as sys
import dendro
from sympy import *

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

# Additional parameters for damping term
R0 = symbols('BSSN_ETA_R0')
ep1, ep2 = symbols('BSSN_ETA_POWER[0] BSSN_ETA_POWER[1]')

xi1, xi2, xi3 = symbols('BSSN_XI[0] BSSN_XI[1] BSSN_XI[2] ')

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

# specify the functions for computing first and second derivatives

# d = dendro.set_first_derivative('grad')    # first argument is direction
# d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
# ad = dendro.set_advective_derivative('agrad')  # first argument is direction
# kod = dendro.set_kreiss_oliger_dissipation('kograd')
dendro.d   = lambda i,x : symbols("grad_%d_%s"%(i,x))
dendro.d2  = lambda i,j,x : symbols("grad2_%d_%d_%s"%(min(i,j),max(i,j),x))
dendro.ad  = dendro.d
dendro.kod = dendro.undef

d  = dendro.d
ad = dendro.ad
kod= dendro.kod
d2 = dendro.d2

dendro.set_metric(gt)
igt = dendro.get_inverse_metric()


eta_func = R0*sqrt(sum([igt[i,j]*d(i,chi)*d(j,chi) for i,j in dendro.e_ij]))/((1-chi**ep1)**ep2)

'''
BSSN puncture gauge (HAD/ traditional BSSN puncture gaugue) with const eta damping 
'''
def bssn_puncture_gauge(eta_damp , isStaged=False , prefix=""):

    if(not isStaged):

        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)
        
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

        ###################################################################
        # generate code
        ###################################################################

        outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
        vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
        dendro.generate_cpu(outs, vnames, '[pp]')


    else:
        # note: these are just the symbolic vars that is being used to generate the
        # Gt_rhs by satges

        _Gt_rhs_s1  = dendro.vec3("Gt_rhs_s1_", "[pp]")
        _Gt_rhs_s2  = dendro.vec3("Gt_rhs_s2_", "[pp]")
        _Gt_rhs_s3  = dendro.vec3("Gt_rhs_s3_", "[pp]")
        _Gt_rhs_s4  = dendro.vec3("Gt_rhs_s4_", "[pp]")
        _Gt_rhs_s5  = dendro.vec3("Gt_rhs_s5_", "[pp]")
        _Gt_rhs_s6  = dendro.vec3("Gt_rhs_s6_", "[pp]")
        _Gt_rhs_s7  = dendro.vec3("Gt_rhs_s7_", "[pp]")
        _CalGt  = dendro.vec3("CalGt", "[pp]")
        _Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")


        # Gt_rhs staged vars that is being used to generate the code.
        At_UU  = dendro.sym_3x3("At_UU", "[pp]")
        CalGt  = dendro.vec3("CalGt", "[pp]")
        Gt_rhs_s1  = dendro.vec3("Gt_rhs_s1_", "[pp]")
        Gt_rhs_s2  = dendro.vec3("Gt_rhs_s2_", "[pp]")
        Gt_rhs_s3  = dendro.vec3("Gt_rhs_s3_", "[pp]")
        Gt_rhs_s4  = dendro.vec3("Gt_rhs_s4_", "[pp]")
        Gt_rhs_s5  = dendro.vec3("Gt_rhs_s5_", "[pp]")
        Gt_rhs_s6  = dendro.vec3("Gt_rhs_s6_", "[pp]")
        Gt_rhs_s7  = dendro.vec3("Gt_rhs_s7_", "[pp]")

        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)
        
        a_rhs = l1*dendro.lie(b, a) - 2*a*K 
        
        b_rhs = [(Rational(3,4) * (lf0 + lf1*a) * B[i] + l2 * dendro.vec_j_ad_j(b, b[i])) for i in dendro.e_i ] 
            
        gt_rhs = dendro.lie(b, gt, weight) - 2*a*At 

        chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) 

        AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

        At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) 

        K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) 

        At_UU = dendro.up_up(At)


        Gt_rhs_s1= ([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i])
        Gt_rhs_s2= ([sum(_CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i])
        Gt_rhs_s3= ([ _CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ])
        Gt_rhs_s4= ([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i])
        Gt_rhs_s5= ([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i])
        Gt_rhs_s6= ([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i])
        Gt_rhs_s7= ([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])

        Gt_rhs = Matrix(_Gt_rhs_s1) - \
            Matrix(_Gt_rhs_s2) + \
            Rational(2,3)*Matrix(_Gt_rhs_s3) + \
            Matrix(_Gt_rhs_s4) - \
            Matrix(_Gt_rhs_s5) + \
            Matrix(_Gt_rhs_s6) - \
            Matrix(_Gt_rhs_s7)

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        B_rhs = [ (Gt_rhs[i] - eta_damp * B[i] +
                l3 * dendro.vec_j_ad_j(b, B[i]) -
                l4 * dendro.vec_j_ad_j(b, Gt[i]))
                for i in dendro.e_i ]

        outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, CalGt, Gt_rhs_s1, Gt_rhs_s2, Gt_rhs_s3, Gt_rhs_s4, Gt_rhs_s5, Gt_rhs_s6, Gt_rhs_s7, Gt_rhs, B_rhs]
        vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'CalGt', 'Gt_rhs_s1_', 'Gt_rhs_s2_', 'Gt_rhs_s3_', 'Gt_rhs_s4_', 'Gt_rhs_s5_', 'Gt_rhs_s6_', 'Gt_rhs_s7_', 'Gt_rhs', 'B_rhs']
        
        ###################################################################
        # generate code
        ###################################################################

        numVars=len(outs)
        for i in range(0,numVars):
            dendro.generate_separate([outs[i]],[vnames[i]],'[pp]')

        

'''
 Uses Rochester puncture gauge.     
'''        
def bssn_rochester_puncture_gauge(eta_damp,isStaged=False,prefix=""):

    if(not isStaged):

        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

        a_rhs = l1*dendro.lie(b, a) - 2*a*K 
    
        b_rhs = [( xi2 * dendro.vec_j_ad_j(b, b[i])  + Rational(3,4) * xi3 * Gt[i] - eta_damp*b[i]) for i in dendro.e_i ] 
        
        gt_rhs = dendro.lie(b, gt, weight) - 2*a*At 

        chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) 

        AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

        At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) 

        K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) 

        At_UU = dendro.up_up(At)

        B_rhs = 0
        
        Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
                Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
                Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
                Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
                Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
                Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
                Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])
            

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        ###################################################################
        # generate code
        ###################################################################

        outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs]
        vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs']
        dendro.generate_cpu(outs, vnames, '[pp]')

    else:
        
        _Gt_rhs_s1  = dendro.vec3("Gt_rhs_s1_", "[pp]")
        _Gt_rhs_s2  = dendro.vec3("Gt_rhs_s2_", "[pp]")
        _Gt_rhs_s3  = dendro.vec3("Gt_rhs_s3_", "[pp]")
        _Gt_rhs_s4  = dendro.vec3("Gt_rhs_s4_", "[pp]")
        _Gt_rhs_s5  = dendro.vec3("Gt_rhs_s5_", "[pp]")
        _Gt_rhs_s6  = dendro.vec3("Gt_rhs_s6_", "[pp]")
        _Gt_rhs_s7  = dendro.vec3("Gt_rhs_s7_", "[pp]")
        _CalGt  = dendro.vec3("CalGt", "[pp]")
        _Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")


        # Gt_rhs staged vars that is being used to generate the code.
        At_UU  = dendro.sym_3x3("At_UU", "[pp]")
        CalGt  = dendro.vec3("CalGt", "[pp]")
        Gt_rhs_s1  = dendro.vec3("Gt_rhs_s1_", "[pp]")
        Gt_rhs_s2  = dendro.vec3("Gt_rhs_s2_", "[pp]")
        Gt_rhs_s3  = dendro.vec3("Gt_rhs_s3_", "[pp]")
        Gt_rhs_s4  = dendro.vec3("Gt_rhs_s4_", "[pp]")
        Gt_rhs_s5  = dendro.vec3("Gt_rhs_s5_", "[pp]")
        Gt_rhs_s6  = dendro.vec3("Gt_rhs_s6_", "[pp]")
        Gt_rhs_s7  = dendro.vec3("Gt_rhs_s7_", "[pp]")

        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

        a_rhs = l1*dendro.lie(b, a) - 2*a*K 
    
        b_rhs = [( xi2 * dendro.vec_j_ad_j(b, b[i])  + Rational(3,4) * xi3 * Gt[i] - eta_damp*b[i]) for i in dendro.e_i ] 
        
        gt_rhs = dendro.lie(b, gt, weight) - 2*a*At 

        chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) 

        AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

        At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) 

        K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) 

        At_UU = dendro.up_up(At)

        B_rhs = 0

        Gt_rhs_s1= ([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i])
        Gt_rhs_s2= ([sum(_CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i])
        Gt_rhs_s3= ([ _CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ])
        Gt_rhs_s4= ([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i])
        Gt_rhs_s5= ([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i])
        Gt_rhs_s6= ([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i])
        Gt_rhs_s7= ([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])

        Gt_rhs = Matrix(_Gt_rhs_s1) - \
            Matrix(_Gt_rhs_s2) + \
            Rational(2,3)*Matrix(_Gt_rhs_s3) + \
            Matrix(_Gt_rhs_s4) - \
            Matrix(_Gt_rhs_s5) + \
            Matrix(_Gt_rhs_s6) - \
            Matrix(_Gt_rhs_s7)

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, CalGt, Gt_rhs_s1, Gt_rhs_s2, Gt_rhs_s3, Gt_rhs_s4, Gt_rhs_s5, Gt_rhs_s6, Gt_rhs_s7, Gt_rhs]
        vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'CalGt', 'Gt_rhs_s1_', 'Gt_rhs_s2_', 'Gt_rhs_s3_', 'Gt_rhs_s4_', 'Gt_rhs_s5_', 'Gt_rhs_s6_', 'Gt_rhs_s7_', 'Gt_rhs']
        
        ###################################################################
        # generate code
        ###################################################################

        numVars=len(outs)
        for i in range(0,numVars):
            dendro.generate_separate([outs[i]],[vnames[i]],'[pp]')



def main():
    if(len(sys.argv)<4):
        print("Error in the bssn code generation.")
        print("usage: python3 bssn.py type[staged|unstaged] gauge[standard, rochester] eta_damp[const|func] prefix[folder parth for staged version]")
        sys.exit(0);

    if(sys.argv[1]=="staged"):
        
        print("//Codgen: generating staged version ")
        if(sys.argv[2]== "rochester"):
            print("//Codgen: using rochester gauge")
            if(sys.argv[3]=="func"):
                print("//Codgen: using eta func damping")
                bssn_rochester_puncture_gauge(eta_func,True,sys.argv[4])
            else:
                print("//Codgen: using eta const damping")
                bssn_rochester_puncture_gauge(eta,True,sys.argv[4])

        else:
            print("//Codgen: using standard gauge")
            if(sys.argv[3]=="func"):
                print("//Codgen: using eta func damping")
                bssn_puncture_gauge(eta_func,True,sys.argv[4])
            else:
                print("//Codgen: using eta const damping")
                bssn_puncture_gauge(eta,True,sys.argv[4])


    else:
        print("//Codgen: generating unstage version ")
        if(sys.argv[2]== "rochester"):
            print("//Codgen: using rochester gauge")
            if(sys.argv[3]=="func"):
                print("//Codgen: using eta func damping")
                bssn_rochester_puncture_gauge(eta_func,False,sys.argv[4])
            else:
                print("//Codgen: using eta const damping")
                bssn_rochester_puncture_gauge(eta,False,sys.argv[4])

        else:
            print("//Codgen: using standard gauge")
            if(sys.argv[3]=="func"):
                print("//Codgen: using eta func damping")
                bssn_puncture_gauge(eta_func,False,sys.argv[4])
            else:
                print("//Codgen: using eta const damping")
                bssn_puncture_gauge(eta,False,sys.argv[4])

        
if __name__ == "__main__":
    main()