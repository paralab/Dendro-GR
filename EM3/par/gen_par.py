# @brief simple python script to generate em3 parameter files for the convergence study. 

import argparse
import resource
import gc
import json as js
from pprint import pprint
import sys as sys


def main():

    if(len(sys.argv) == 0):
        print ("Error: em3 template par file parameter file not specified")
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
    print("#SBATCH -o em3-%j.out-%N ")
    print("#SBATCH --ntasks=%d" %(nodes*cores_per_node))
    print("#SBATCH --account=%s" %(partition) )
    print("#SBATCH --partition=%s" %(account))

    work_dir ="/uufs/chpc.utah.edu/common/home/u1011531/Research/Dendro-5.0/build"
    scratch_dir ="/scratch/kingspeak/serial/u1011531/em3"

    print("cd %s" %work_dir)
    print("make em3Solver -j28")
    print("cp EM3/em3Solver %s/." %scratch_dir)
    print("cp *.par.json %s/." %scratch_dir)
    print("cp *.vis.par.json %s/." %scratch_dir)
    print("cp *.py %s/." %scratch_dir)
    print("cd %s" %scratch_dir)
    
    vis_files=list()

    for p in pOrder:
        for g in grid:
            for d in depth:
                for c in cfl:
                    params["EM3_ELE_ORDER"] =p

                    params["EM3_MAXDEPTH"] = d
                    vparams["EM3_MAXDEPTH"] = d
                    
                    params["EM3_CFL_FACTOR"] = c
                    vparams["EM3_CFL_FACTOR"] = c

                    params["EM3_VTU_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em3/vtu/em3_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    vparams["EM3_VTU_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em3/vtu/em3_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    vparams["EM3_IMG_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em3/img/em3_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    
                    params["EM3_CHKPT_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em3/cp/em3_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_" 
                    
                    
                    params["EM3_PROFILE_FILE_PREFIX"] = "/scratch/kingspeak/serial/u1011531/em3/em3_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"
                    
                    params["EM3_WAVELET_TOL"] = base_wtol / ((1/2)**d)
                    vparams["EM3_WAVELET_TOL"] = base_wtol / ((1/2)**d)

                    if (g=="box"):
                        params["EM3_ENABLE_BLOCK_ADAPTIVITY"] = 1
                    else:
                        params["EM3_ENABLE_BLOCK_ADAPTIVITY"] = 0
                        if(g == "wamr"):
                            params["EM3_EM3_REFINE_MODE"]=0
                        elif( g == "fr"):
                            params["EM3_EM3_REFINE_MODE"]=1
                        elif ( g == "wamr_fr"):
                            params["EM3_EM3_REFINE_MODE"]=2

                    
                    pfname = "em3_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"+"par.json" 
                    with open(pfname, 'w') as outfile:
                        js.dump(params, outfile,indent=4)

                    pfname = "em3_" +g+"_d"+str(d)+"_cfl"+str(c)+"_e"+str(p)+"_"+".vis.par.json" 
                    vis_files.append(pfname)
                    with open(pfname, 'w') as outfile:
                        js.dump(vparams, outfile,indent=4)

                    print("\n")
                    print("mpirun -np $SLURM_NTASKS ./EM3/em3Solver "+pfname)
                    print("\n")
    

    for v in vis_files:
        print("mpirun -np 56 paraview-mesa pvbatch -- --force-offscreen-rendering em3_vis.py "+ v +" 0 600 50")

if __name__ == "__main__":
    main()



                



    
