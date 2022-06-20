# first derivatives
# D = ["alpha", "beta0", "beta1", "beta2",
#      "B0", "B1", "B2",
#      "chi", "Gt0", "Gt1", "Gt2", "K",
#      "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
#      "At0", "At1", "At2", "At3", "At4", "At5" ]

# same order as in grDef.h
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

DX="device::__deriv644_x"
DY="device::__deriv644_y"
DZ="device::__deriv644_z"

DXX="device::__deriv644_xx"
DYY="device::__deriv644_yy"
DZZ="device::__deriv644_zz"

KDX = "device::__ko_deriv42_x"
KDY = "device::__ko_deriv42_y"
KDZ = "device::__ko_deriv42_z"


deriv_dec   = ""
folderpath = "../include/"
filename   = folderpath + "bssnrhs_evar_derivs.cuh"
with open(filename, 'w') as out:
    out.write("#pragma once\n")
    out.write("#include \"device.h\"\n")
    out.write("#include <stdio.h>\n")
    for deriv in FUNC_D_I:
        deriv_dec += "\tDEVICE_REAL *"+deriv+";\n" 
    
    for deriv in FUNC_D_IJ:
        deriv_dec += "\tDEVICE_REAL *"+deriv+";\n"

    for deriv in FUNC_KOD_I:
        deriv_dec += "\tDEVICE_REAL *"+deriv+";\n"

    out.write("struct BSSN_EVAR_DERIVS{\n"+deriv_dec+"};\n")


folderpath = "./"
filename   = folderpath + "bssnrhs_memalloc.h"
deriv_evar_struct = "deriv_evars"
with open(filename,'w') as out:
    num_derivs = len(FUNC_D_I) + len(FUNC_D_IJ)
    
    deriv_id=0
    for deriv in FUNC_D_I:
        out.write("double * "+deriv+" = deriv_base + %d * BLK_SZ;\n"%deriv_id)
        deriv_id+=1
    
    for deriv in FUNC_D_IJ:
        out.write("double * "+deriv+" = deriv_base + %d * BLK_SZ;\n"%deriv_id)
        deriv_id+=1

    for deriv in FUNC_KOD_I:
        out.write("double * "+deriv+" = deriv_base + %d * BLK_SZ;\n"%deriv_id)
        deriv_id+=1

# filename   = folderpath + "bssnrhs_dealloc.cuh"
# with open(filename, 'w') as out:
#     out.write("GPUDevice::device_free(deriv_base);\n")
#     # for deriv in FUNC_D_I:
#     #     out.write("GPUDevice::device_free(" + deriv_evar_struct+"->" + deriv + ");\n")

#     # for deriv in FUNC_D_IJ:
#     #     out.write("GPUDevice::device_free(" + deriv_evar_struct+"->" + deriv + ");\n")

#     # for deriv in FUNC_KOD_I:
#     #     out.write("GPUDevice::device_free(" + deriv_evar_struct+"->" + deriv + ");\n")

filename  = folderpath + "bssnrhs_deriv_x_dir.cuh"
with open(filename, 'w') as out:
    for var in D:
        out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID);\n" %(DX, deriv_evar_struct+"->"+PREFIX_D[0] + var ,var))
        out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID);\n" %(KDX ,deriv_evar_struct+"->"+PREFIX_KOD[0] + var ,var))
        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID);\n" %(DXX,deriv_evar_struct+"->"+PREFIX_DD[0] + var ,var))
        
            

filename  = folderpath + "bssnrhs_deriv_y_dir.cuh"
with open(filename, 'w') as out:
    for var in D:
        out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID );\n" %(DY,deriv_evar_struct+"->"+PREFIX_D[1] + var ,var))
        out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID );\n" %(KDY,deriv_evar_struct+"->"+PREFIX_KOD[1] + var ,var))
        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID);\n" %(DYY,deriv_evar_struct+"->"+PREFIX_DD[3] + var ,var))
            
            out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s + BLK_ID * BLK_SZ, blk + BLK_ID);\n" %(DY,deriv_evar_struct+"->"+PREFIX_DD[1] + var , deriv_evar_struct+"->" + PREFIX_D[0] + var ))

filename  = folderpath + "bssnrhs_deriv_z_dir.cuh"
with open(filename, 'w') as out:
    for var in D:
        out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID);\n" %(DZ,deriv_evar_struct+"->"+PREFIX_D[2] + var ,var))
        out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID);\n" %(KDZ,deriv_evar_struct+"->"+PREFIX_KOD[2] + var ,var))
        
        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s, blk + BLK_ID);\n" %(DZZ,deriv_evar_struct+"->"+PREFIX_DD[5] + var ,var))

            out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s + BLK_ID * BLK_SZ, blk + BLK_ID);\n" %(DZ,deriv_evar_struct+"->"+PREFIX_DD[2] + var , deriv_evar_struct+"->" + PREFIX_D[0] + var ))
            out.write("%s<pw,pencils,pencil_sz>(%s + BLK_ID * BLK_SZ, %s + BLK_ID * BLK_SZ, blk + BLK_ID);\n" %(DZ,deriv_evar_struct+"->"+PREFIX_DD[4] + var , deriv_evar_struct+"->" + PREFIX_D[1] + var))
        

# renaming the bssn derivs to compute the rhs. 
folderpath = "./"
filename   = folderpath + "bssnrhs_evar_derivs.h"
with open(filename, 'w') as out:
    for deriv in FUNC_D_I:
        out.write("const DEVICE_REAL %s = *(%s->%s + BLK_ID * BLK_SZ + pp);\n"%(deriv,deriv_evar_struct,deriv))
        
    for deriv in FUNC_D_IJ:
        out.write("const DEVICE_REAL %s = *(%s->%s + BLK_ID * BLK_SZ + pp);\n"%(deriv,deriv_evar_struct,deriv))

    for deriv in FUNC_KOD_I:
        out.write("const DEVICE_REAL %s = *(%s->%s + BLK_ID * BLK_SZ + pp);\n"%(deriv,deriv_evar_struct,deriv))

filename  = folderpath + "compute_bssnrhs_evar_derivs1.cuh"
with open(filename, 'w') as out:
    for var in D:
        SV=var
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(DX, deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(DY, deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_D[2] +  var  + " + Z_ID * BLK_SZ", SV))

        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(KDX,deriv_evar_struct+"->"+PREFIX_KOD[0] + var + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(KDY,deriv_evar_struct+"->"+PREFIX_KOD[1] + var + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(KDZ,deriv_evar_struct+"->"+PREFIX_KOD[2] + var + " + Z_ID * BLK_SZ", SV))
        
        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk + BLK_ID);\n" %(DXX,deriv_evar_struct+"->"+PREFIX_DD[0] + var + " + Z_ID * BLK_SZ", SV))
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk + BLK_ID);\n" %(DYY,deriv_evar_struct+"->"+PREFIX_DD[3] + var + " + Z_ID * BLK_SZ", SV))
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk + BLK_ID);\n" %(DZZ,deriv_evar_struct+"->"+PREFIX_DD[5] + var + " + Z_ID * BLK_SZ", SV))

    out.write("// mixed deriv pass \n")
    out.write("g.sync();\n")
    for var in DD:
        SV=deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ"
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(DY, deriv_evar_struct+"->"+PREFIX_DD[1] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_DD[2] +  var  + " + Z_ID * BLK_SZ", SV))
        
        SV=deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ"
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk + BLK_ID);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_DD[4] +  var  + " + Z_ID * BLK_SZ", SV))
            

DX  = "device::__blk_deriv644_x"
DY  = "device::__blk_deriv644_y"
DZ  = "device::__blk_deriv644_z"

DXX = "device::__blk_deriv644_xx"
DYY = "device::__blk_deriv644_yy"
DZZ = "device::__blk_deriv644_zz"

KDX = "device::__blk_ko_deriv42_x"
KDY = "device::__blk_ko_deriv42_y"
KDZ = "device::__blk_ko_deriv42_z"

LD  = "device::__load_blk_var__"
SV  = "su"
ST  = "__syncthreads();"

folderpath = "./"
filename   = folderpath + "bssnrhs_memalloc.cuh"
deriv_evar_struct = "deriv_evars"
deriv_id=0
with open(filename,'w') as out:
    for var in D:
        out.write(deriv_evar_struct+"->"+(PREFIX_D[0] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id))
        out.write(deriv_evar_struct+"->"+(PREFIX_D[1] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+1))
        out.write(deriv_evar_struct+"->"+(PREFIX_D[2] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+2))

        out.write(deriv_evar_struct+"->"+(PREFIX_KOD[0] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+3))
        out.write(deriv_evar_struct+"->"+(PREFIX_KOD[1] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+4))
        out.write(deriv_evar_struct+"->"+(PREFIX_KOD[2] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+5))
        deriv_id+=6

        if var in DD:
            out.write(deriv_evar_struct+"->"+(PREFIX_DD[0] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id))
            out.write(deriv_evar_struct+"->"+(PREFIX_DD[1] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+1))
            out.write(deriv_evar_struct+"->"+(PREFIX_DD[2] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+2))
            out.write(deriv_evar_struct+"->"+(PREFIX_DD[3] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+3))
            out.write(deriv_evar_struct+"->"+(PREFIX_DD[4] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+4))
            out.write(deriv_evar_struct+"->"+(PREFIX_DD[5] +  var)+" = deriv_base + %d * BATCHED_BLOCKS_SZ * BLK_SZ;\n"%(deriv_id+5))
            deriv_id+=6


folderpath = "./"
filename  = folderpath + "compute_bssnrhs_evar_derivs2.cuh"
with open(filename, 'w') as out:
    for i, var in enumerate(D):
        out.write("if(V_ID == %d) { \n"%i)
        out.write("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV,"var_in"))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDX,deriv_evar_struct+"->"+PREFIX_KOD[0] + var + " + Z_ID * BLK_SZ", SV))
        
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDY,deriv_evar_struct+"->"+PREFIX_KOD[1] + var + " + Z_ID * BLK_SZ", SV))
        
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_D[2] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDZ,deriv_evar_struct+"->"+PREFIX_KOD[2] + var + " + Z_ID * BLK_SZ", SV))

        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk);\n" %(DXX,deriv_evar_struct+"->"+PREFIX_DD[0] + var + " + Z_ID * BLK_SZ", SV))
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk);\n" %(DYY,deriv_evar_struct+"->"+PREFIX_DD[3] + var + " + Z_ID * BLK_SZ", SV))
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk);\n" %(DZZ,deriv_evar_struct+"->"+PREFIX_DD[5] + var + " + Z_ID * BLK_SZ", SV))

            out.write("%s\n"%ST)
            out.write("%s<DEVICE_REAL,pw,nx>(%s , %s, blk);\n"%(LD, SV,deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ"))
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, deriv_evar_struct+"->"+PREFIX_DD[1] +  var  + " + Z_ID * BLK_SZ", SV))
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_DD[2] +  var  + " + Z_ID * BLK_SZ", SV))

            out.write("%s\n"%ST)
            out.write("%s<DEVICE_REAL,pw,nx>(%s , %s, blk);\n"%(LD, SV,deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ"))
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_DD[4] +  var  + " + Z_ID * BLK_SZ", SV))

        out.write("return;\n")
        out.write("} // end of var %s \n"%var)

folderpath = "./"
filename  = folderpath + "compute_bssnrhs_evar_derivs.cuh"
with open(filename, 'w') as out:
    for var in D:
        out.write("%s<DEVICE_REAL,pw,nx>(%s , %s, blk);\n"%(LD, SV,var))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDX,deriv_evar_struct+"->"+PREFIX_KOD[0] + var + " + Z_ID * BLK_SZ", SV))

        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk);\n" %(DXX,deriv_evar_struct+"->"+PREFIX_DD[0] + var + " + Z_ID * BLK_SZ", SV))
        

        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDY,deriv_evar_struct+"->"+PREFIX_KOD[1] + var + " + Z_ID * BLK_SZ", SV))

        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk);\n" %(DYY,deriv_evar_struct+"->"+PREFIX_DD[3] + var + " + Z_ID * BLK_SZ", SV))

        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_D[2] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDZ,deriv_evar_struct+"->"+PREFIX_KOD[2] + var + " + Z_ID * BLK_SZ", SV))
        
        if var in DD:
            out.write("%s<pw,pencils,pencil_sz>(%s, %s, blk);\n" %(DZZ,deriv_evar_struct+"->"+PREFIX_DD[5] + var + " + Z_ID * BLK_SZ", SV))

        out.write("%s\n"%ST)

    out.write("// mixed deriv pass \n")
    out.write("%s\n"%ST)
    for var in DD:
        out.write("%s<DEVICE_REAL,pw,nx>(%s , %s, blk);\n"%(LD, SV,deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ"))

        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, deriv_evar_struct+"->"+PREFIX_DD[1] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_DD[2] +  var  + " + Z_ID * BLK_SZ", SV))

        out.write("%s\n"%ST)
        
        out.write("%s<DEVICE_REAL,pw,nx>(%s , %s, blk);\n"%(LD, SV,deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ"))
        out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, deriv_evar_struct+"->"+PREFIX_DD[4] +  var  + " + Z_ID * BLK_SZ", SV))
        out.write("%s\n"%ST)

folderpath = "./"
filename   = folderpath + "bssnrhs_evar_derivs_zid.h"
with open(filename, 'w') as out:
    for deriv in FUNC_D_I:
        out.write("const DEVICE_REAL * __restrict__ const %s = (%s->%s + Z_ID * BLK_SZ);\n"%(deriv,deriv_evar_struct,deriv))
        
    for deriv in FUNC_D_IJ:
        out.write("const DEVICE_REAL * __restrict__ const %s = (%s->%s + Z_ID * BLK_SZ);\n"%(deriv,deriv_evar_struct,deriv))

    for deriv in FUNC_KOD_I:
        out.write("const DEVICE_REAL * __restrict__ const %s = (%s->%s + Z_ID * BLK_SZ);\n"%(deriv,deriv_evar_struct,deriv))



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


folderpath = "./"
filename  = folderpath + "compute_bssnrhs_evar_derivs3.cuh"
with open(filename, 'w') as out:
    for i, var in enumerate(D):
        if var in DD:
            # load
            out.write("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV, var))
            #xx
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DXX, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[0] +  var, DSV))
            out.write("%s\n"%ST)

            # kox
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDX, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[0] +  var, DSV))
            out.write("%s\n"%ST)
            
            #x
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[0] +  var, DSV))
            out.write("%s\n"%ST)

            #xy
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, DDSV, DSV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[1] +  var, DDSV))
            out.write("%s\n"%ST)

            #xz
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DDSV, DSV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[2] +  var, DDSV))
            out.write("%s\n"%ST)

            #yy
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DYY, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[3] +  var, DSV))
            out.write("%s\n"%ST)

            # koy
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDY, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[1] +  var, DSV))
            out.write("%s\n"%ST)

            # y 
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[1] +  var, DSV))
            out.write("%s\n"%ST)

            # #yz
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DDSV, DSV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[4] +  var, DDSV))
            out.write("%s\n"%ST)

            #zz
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZZ, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_DD[5] +  var, DSV))
            out.write("%s\n"%ST)

            # koz
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDZ, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[2] +  var, DSV))
            out.write("%s\n"%ST)

            #z
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[2] +  var, DSV))
            out.write("%s\n"%ST)

            


        else:

            out.write("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV, var))
            
            # x
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[0] +  var, DSV))
            out.write("%s\n"%ST)

            # kox
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDX, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[0] +  var, DSV))
            out.write("%s\n"%ST)

            # y
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[1] +  var, DSV))
            out.write("%s\n"%ST)

            # koy
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDY, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[1] +  var, DSV))
            out.write("%s\n"%ST)

            # z
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_D[2] +  var, DSV))
            out.write("%s\n"%ST)

            # koz
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDZ, DSV, SV))
            out.write("%s\n"%ST)
            out.write("const DEVICE_REAL %s=%s[gidx];\n"%(PREFIX_KOD[2] +  var, DSV))
            out.write("%s\n"%ST)


folderpath = "./"
filename  = folderpath + "compute_bssnrhs_evar_derivs4.cuh"
with open(filename, 'w') as out:
    for i, var in enumerate(D):
        if var in DD:
            # load
            out.write("%s\n"%ST)
            out.write("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV, var))
            #xx
            dvar =  deriv_evar_struct+"->"+PREFIX_DD[0] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DXX, dvar, SV))
            

            # kox
            dvar =  deriv_evar_struct+"->"+PREFIX_KOD[0] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDX, dvar, SV))
            
            
            #x
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, DSV, SV))
            out.write("%s\n"%ST)
            dvar =  deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(SD, DSV, dvar))

            #xy
            dvar =  deriv_evar_struct+"->"+PREFIX_DD[1] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, dvar, DSV))
            
            #xz
            dvar =  deriv_evar_struct+"->"+PREFIX_DD[2] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, dvar, DSV))
            
            #yy
            dvar =  deriv_evar_struct+"->"+PREFIX_DD[3] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DYY, dvar, SV))
            
            #koy
            dvar =  deriv_evar_struct+"->"+PREFIX_KOD[1] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDY, dvar, SV))
            
            out.write("%s\n"%ST)
            #y 
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, DSV, SV))
            out.write("%s\n"%ST)
            dvar =  deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(SD, DSV, dvar))

            #yz
            dvar =  deriv_evar_struct+"->"+PREFIX_DD[4] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, dvar, DSV))
            
            #zz
            dvar =  deriv_evar_struct+"->"+PREFIX_DD[5] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZZ, dvar, SV))
            
            #koz
            dvar =  deriv_evar_struct+"->"+PREFIX_KOD[2] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDZ, dvar, SV))
            
            #z
            dvar =  deriv_evar_struct+"->"+PREFIX_D[2] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, dvar, SV))
            
        else:
            out.write("%s\n"%ST)
            out.write("%s<DEVICE_REAL,pw,nx>   (%s , %s, blk);\n"%(LD, SV, var))
            
            # x
            dvar =  deriv_evar_struct+"->"+PREFIX_D[0] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DX, dvar, SV))
            
            # kox
            dvar =  deriv_evar_struct+"->"+PREFIX_KOD[0] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDX, dvar, SV))
            
            # y
            dvar =  deriv_evar_struct+"->"+PREFIX_D[1] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DY, dvar, SV))
            
            # koy
            dvar =  deriv_evar_struct+"->"+PREFIX_KOD[1] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDY, dvar, SV))
            
            # z
            dvar =  deriv_evar_struct+"->"+PREFIX_D[2] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(DZ, dvar, SV))
            
            # koz
            dvar =  deriv_evar_struct+"->"+PREFIX_KOD[2] +  var  + " + Z_ID * BLK_SZ"
            out.write("%s<pw,pencils,pencil_sz>(%s , %s, blk);\n" %(KDZ, dvar, SV))
            
        

        