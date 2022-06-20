# @brief simple python script to generate nlsm parameter files for the convergence study. 

import argparse
import resource
import gc
import json as js
from pprint import pprint
import sys as sys


def main():

    if(len(sys.argv) == 0):
        print ("Error: nlsm template par file parameter file not specified")
        sys.exit(0)

    with open(sys.argv[1]) as f:
        params = js.load(f)

    with open(sys.argv[2]) as f:
        vparams = js.load(f)


    pOrder = [4,8]
    depth = [8,9,10]
    #cfl = [0.125,0.25,0.5]
    cfl = [0.25]
    grid =["wamr"]

    partition ="soc-kp"
    account = "soc-kp"
    nodes =10
    cores_per_node=28
    base_wtol = 1e-5


    print("#!/bin/bash")
    print("#SBATCH --time=72:00:00")
    print("#SBATCH --nodes=%d" %(nodes))
    print("#SBATCH -o nlsm-%j.out-%N ")
    print("#SBATCH --ntasks=%d" %(nodes*cores_per_node))
    print("#SBATCH --account=%s" %(partition) )
    print("#SBATCH --partition=%s" %(account))

    work_dir ="/uufs/chpc.utah.edu/common/home/u1011531/Research/Dendro-5.0/build"
    scratch_dir ="/scratch/kingspeak/serial/u1011531/nlsm"

    print("cd %s" %work_dir)
    print("make nlsmSolver -j28")
    print("cp NLSigma/nlsmSolver %s/." %scratch_dir)
    print("cp *.par.json %s/." %scratch_dir)
    print("cp *.vis.par.json %s/." %scratch_dir)
    print("cp *.py %s/." %scratch_dir)
    print("cd %s" %scratch_dir)
    
    vis_files=list()

    for p in pOrder:
        for g in grid:
            for d in depth:
                for c in cfl:
                    params["NLSM_ELE_ORDER"] =p

                    params["NLSM_MAXDEPTH"] = d
                    vparams["NLSM_MAXDEPTH"] = d
                    
                    params["NLSM_CFL_FACTOR"] = c
                    vparams["NLSM_CFL_FACTOR"] = c

                    params["NLSM_VTU_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/nlsm/vtu/nlsm_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    vparams["NLSM_VTU_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/nlsm/vtu/nlsm_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    vparams["NLSM_IMG_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/nlsm/img/nlsm_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    
                    params["NLSM_CHKPT_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/nlsm/cp/nlsm_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    
                    
                    params["NLSM_PROFILE_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/nlsm/nlsm_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"
                    
                    params["NLSM_WAVELET_TOL"] = base_wtol / ((1/2)**d)
                    vparams["NLSM_WAVELET_TOL"] = base_wtol / ((1/2)**d)

                    if (g=="box"):
                        params["NLSM_ENABLE_BLOCK_ADAPTIVITY"] = 1
                    else:
                        params["NLSM_ENABLE_BLOCK_ADAPTIVITY"] = 0
                        if(g == "wamr"):
                            params["NLSM_NLSM_REFINE_MODE"]=0
                        elif( g == "fr"):
                            params["NLSM_NLSM_REFINE_MODE"]=1
                        elif ( g == "wamr_fr"):
                            params["NLSM_NLSM_REFINE_MODE"]=2

                    
                    pfname = "nlsm_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"+"par.json" 
                    with open(pfname, 'w') as outfile:
                        js.dump(params, outfile,indent=4)

                    pfname = "nlsm_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"+".vis.par.json" 
                    vis_files.append(pfname)
                    with open(pfname, 'w') as outfile:
                        js.dump(vparams, outfile,indent=4)

                    print("\n")
                    print("mpirun -np $SLURM_NTASKS ./NLSigma/nlsmSolver "+pfname)
                    print("\n")
    

    for v in vis_files:
        print("mpirun -np 56 paraview-mesa pvbatch -- --force-offscreen-rendering nlsm_vis.py "+ v +" 0 600 50")

if __name__ == "__main__":
    main()



                



    
