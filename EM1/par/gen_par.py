# @brief simple python script to generate em1 parameter files for the convergence study. 

import argparse
import resource
import gc
import json as js
from pprint import pprint
import sys as sys


def main():

    if(len(sys.argv) == 0):
        print ("Error: em1 template par file parameter file not specified")
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
    print("#SBATCH -o em1-%j.out-%N ")
    print("#SBATCH --ntasks=%d" %(nodes*cores_per_node))
    print("#SBATCH --account=%s" %(partition) )
    print("#SBATCH --partition=%s" %(account))

    work_dir ="/uufs/chpc.utah.edu/common/home/u1011531/Research/Dendro-5.0/build"
    scratch_dir ="/scratch/kingspeak/serial/u1011531/em1"

    print("cd %s" %work_dir)
    print("make em1Solver -j28")
    print("cp EM1/em1Solver %s/." %scratch_dir)
    print("cp *.par.json %s/." %scratch_dir)
    print("cp *.vis.par.json %s/." %scratch_dir)
    print("cp *.py %s/." %scratch_dir)
    print("cd %s" %scratch_dir)
    
    vis_files=list()

    for p in pOrder:
        for g in grid:
            for d in depth:
                for c in cfl:
                    params["EM1_ELE_ORDER"] =p

                    params["EM1_MAXDEPTH"] = d
                    vparams["EM1_MAXDEPTH"] = d
                    
                    params["EM1_CFL_FACTOR"] = c
                    vparams["EM1_CFL_FACTOR"] = c

                    params["EM1_VTU_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em1/vtu/em1_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    vparams["EM1_VTU_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em1/vtu/em1_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    vparams["EM1_IMG_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em1/img/em1_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    
                    params["EM1_CHKPT_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em1/cp/em1_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    
                    
                    params["EM1_PROFILE_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em1/em1_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"
                    
                    params["EM1_WAVELET_TOL"] = base_wtol / ((1/2)**d)
                    vparams["EM1_WAVELET_TOL"] = base_wtol / ((1/2)**d)

                    if (g=="box"):
                        params["EM1_ENABLE_BLOCK_ADAPTIVITY"] = 1
                    else:
                        params["EM1_ENABLE_BLOCK_ADAPTIVITY"] = 0
                        if(g == "wamr"):
                            params["EM1_EM1_REFINE_MODE"]=0
                        elif( g == "fr"):
                            params["EM1_EM1_REFINE_MODE"]=1
                        elif ( g == "wamr_fr"):
                            params["EM1_EM1_REFINE_MODE"]=2

                    
                    pfname = "em1_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"+"par.json" 
                    with open(pfname, 'w') as outfile:
                        js.dump(params, outfile,indent=4)

                    pfname = "em1_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"+".vis.par.json" 
                    vis_files.append(pfname)
                    with open(pfname, 'w') as outfile:
                        js.dump(vparams, outfile,indent=4)

                    print("\n")
                    print("mpirun -np $SLURM_NTASKS ./EM1/em1Solver "+pfname)
                    print("\n")
    

    for v in vis_files:
        print("mpirun -np 56 paraview-mesa pvbatch -- --force-offscreen-rendering em1_vis.py "+ v +" 0 600 50")

if __name__ == "__main__":
    main()



                



    
