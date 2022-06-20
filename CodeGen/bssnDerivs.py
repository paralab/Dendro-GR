"""
Contains derivative computation for BSSN formulation of ET equations. 
"""
D = ["alpha", "beta0", "beta1", "beta2",
     "B0", "B1", "B2",
     "chi", "Gt0", "Gt1", "Gt2", "K",
     "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
     "At0", "At1", "At2", "At3", "At4", "At5" ]

# advective derivatives
AD = [ "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
       "At0", "At1", "At2", "At3", "At4", "At5",
       "alpha", "beta0", "beta1", "beta2",
       "chi", "Gt0", "Gt1", "Gt2", "K",
       "B0", "B1", "B2"]

# Kries-Oliger derivatives. 
KO = [ "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
       "At0", "At1", "At2", "At3", "At4", "At5",
       "alpha", "beta0", "beta1", "beta2",
       "chi", "Gt0", "Gt1", "Gt2", "K",
       "B0", "B1", "B2"]

# second derivs required for RHS
DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi",
      "alpha", "beta0", "beta1", "beta2" ] 

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
PREFIX_AD  = ["agrad_0_",   "agrad_1_",     "agrad_2_"]
PREFIX_KOD = ["kograd_0_",  "kograd_1_",    "kograd_2_"]
PREFIX_DD  = ["grad2_0_0_", "grad2_0_1_",   "grad2_0_2_", "grad2_1_1_", "grad2_1_2_", "grad2_2_2_"]

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

#advective derivative in i direction
FUNC_AD_I=[]
for f in AD:
    for p in PREFIX_AD:
        FUNC_AD_I.append(p+f)


#Kriess-Oliger derivative in i direction
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



folderpath = "../scripts/"
filename   = folderpath + "bssnrhs_memalloc.h"
with open(filename, 'w') as out:
    for deriv in FUNC_D_I:
        out.write("  double *"+deriv+" = __mem_pool->allocate(n);\n")
    
    for deriv in FUNC_D_IJ:
        out.write("  double *"+deriv+" = __mem_pool->allocate(n);\n")


filename   = folderpath + "bssnrhs_dealloc.h"
with open(filename, 'w') as out:
    for deriv in FUNC_D_I:
        out.write("__mem_pool->free(" + deriv + ");\n")
    
    for deriv in FUNC_D_IJ:
        out.write("__mem_pool->free(" + deriv + ");\n")


filename   = folderpath + "bssnrhs_memalloc_adv.h"
with open(filename, 'w') as out:
    out.write("#ifdef BSSN_USE_ADVECTIVE_DERIVS\n")
    for deriv in FUNC_AD_I:
        out.write("  double *"+deriv+" = __mem_pool->allocate(n);\n")
    out.write("#else\n")
    for f in AD:
        for p in range(0,len(PREFIX_AD)):
            out.write("  double *"+(PREFIX_AD[p]+f)+ " = "+(PREFIX_D[p]+f)+";\n")
    out.write("#endif\n")

filename   = folderpath + "bssnrhs_dealloc_adv.h"
with open(filename, 'w') as out:
    out.write("#ifdef BSSN_USE_ADVECTIVE_DERIVS\n")
    for deriv in FUNC_AD_I:
        out.write("__mem_pool->free(" + deriv + ");\n")
    out.write("#endif\n")

filename   = folderpath + "bssnrhs_derivs.h"
with open(filename, 'w') as out:
    for var in D:
        if var in DD:
            out.write("  deriv_x(%s, %s, hx, sz, bflag);\n" %(PREFIX_D[0] + var ,var))
            out.write("  deriv_xx(%s, %s, hx, sz, bflag);\n" %(PREFIX_DD[0] + var ,var))

            out.write("  deriv_y(%s, %s, hy, sz, bflag);\n" %(PREFIX_D[1] + var ,var))
            out.write("  deriv_yy(%s, %s, hy, sz, bflag);\n" %(PREFIX_DD[3] + var ,var))

            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_D[2] + var ,var))
            out.write("  deriv_zz(%s, %s, hz, sz, bflag);\n" %(PREFIX_DD[5] + var ,var))

            # mixed 2nd derivs. 
            out.write("  deriv_y(%s, %s, hy, sz, bflag);\n" %(PREFIX_DD[1] + var , PREFIX_D[0] + var ))
            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_DD[2] + var , PREFIX_D[0] + var ))
            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_DD[4] + var , PREFIX_D[1] + var))
        else:
            out.write("  deriv_x(%s, %s, hx, sz, bflag);\n" %(PREFIX_D[0] + var ,var))
            out.write("  deriv_y(%s, %s, hy, sz, bflag);\n" %(PREFIX_D[1] + var ,var))
            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_D[2] + var ,var))



filename   = folderpath + "bssnrhs_derivs_adv.h"
with open(filename, 'w') as out:
    out.write("#ifdef BSSN_USE_ADVECTIVE_DERIVS\n")
    for var in AD:
        out.write("  adv_deriv_x(%s, %s, hx, sz, beta0, bflag);\n" %(PREFIX_AD[0] + var ,var))
        out.write("  adv_deriv_y(%s, %s, hy, sz, beta1, bflag);\n" %(PREFIX_AD[1] + var ,var))
        out.write("  adv_deriv_z(%s, %s, hz, sz, beta2, bflag);\n" %(PREFIX_AD[2] + var ,var))
    out.write("#endif\n")

filename   = folderpath + "bssnrhs_ko_derivs.h"
with open(filename, 'w') as out:
    for var in KO:
        out.write("  ko_deriv_x(%s, %s, hx, sz, bflag);\n" %(PREFIX_D[0] + var ,var))
        out.write("  ko_deriv_y(%s, %s, hy, sz, bflag);\n" %(PREFIX_D[1] + var ,var))
        out.write("  ko_deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_D[2] + var ,var))



## ============= constraint equations =============================

filename   = folderpath + "constraint_memalloc.h"
with open(filename, 'w') as out:
    for deriv in FUNC_CONS:
        out.write("  double *"+deriv+" = __mem_pool->allocate(n);\n")



filename   = folderpath + "constraint_dealloc.h"
with open(filename, 'w') as out:
    for deriv in FUNC_CONS:
        out.write("__mem_pool->free(" + deriv + ");\n")


filename   = folderpath + "constraint_derivs.h"
with open(filename, 'w') as out:
    for var in CONSTRAINT_D:
        if var in CONSTRAINT_DD:
            out.write("  deriv_x(%s, %s, hx, sz, bflag);\n" %(PREFIX_D[0] + var ,var))
            out.write("  deriv_xx(%s, %s, hx, sz, bflag);\n" %(PREFIX_DD[0] + var ,var))

            out.write("  deriv_y(%s, %s, hy, sz, bflag);\n" %(PREFIX_D[1] + var ,var))
            out.write("  deriv_yy(%s, %s, hy, sz, bflag);\n" %(PREFIX_DD[3] + var ,var))

            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_D[2] + var ,var))
            out.write("  deriv_zz(%s, %s, hz, sz, bflag);\n" %(PREFIX_DD[5] + var ,var))

            # mixed 2nd derivs. 
            out.write("  deriv_y(%s, %s, hy, sz, bflag);\n" %(PREFIX_DD[1] + var , PREFIX_D[0] + var ))
            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_DD[2] + var , PREFIX_D[0] + var ))
            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_DD[4] + var , PREFIX_D[1] + var))

        else:
            out.write("  deriv_x(%s, %s, hx, sz, bflag);\n" %(PREFIX_D[0] + var ,var))
            out.write("  deriv_y(%s, %s, hy, sz, bflag);\n" %(PREFIX_D[1] + var ,var))
            out.write("  deriv_z(%s, %s, hz, sz, bflag);\n" %(PREFIX_D[2] + var ,var))

