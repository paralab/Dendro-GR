#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Lonestar6 AMD Milan nodes
#
#   *** MPI Job in Normal Queue ***
# 
# Last revised: October 22, 2021
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch milan.mpi.slurm" on a Lonestar6 login node.
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do NOT use mpirun or mpiexec.
#
#   -- Max recommended MPI ranks per Milan node: 128
#      (start small, increase gradually).
#
#   -- If you're running out of memory, try running
#      fewer tasks per node to give each task more memory.
#
#----------------------------------------------------

#SBATCH -J dgr            # Job name
#SBATCH -o .dgr.o%j       # Name of stdout output file
#SBATCH -e .dgr.e%j       # Name of stderr error file
#SBATCH -p gpu-a100         # Queue (partition) name
#SBATCH -N 16              # Total # of nodes 
#SBATCH -n 32              # Total # of mpi tasks
#SBATCH -t 01:30:00       # Run time (hh:mm:ss)
#SBATCH --mail-type=all   # Send email at begin and end of job
##SBATCH -A myproject     # Project/Allocation name (req'd if you have more than 1)
##SBATCH --mail-user=username@tacc.utexas.edu

# Any other commands must follow all #SBATCH directives...
module list
pwd
date

make bssnSolverCUDA -j4
#ibrun -np 1 ./BSSN_GR/bssnSolverCUDA q1_r2.2.par.json 1
# Launch MPI code... 
ibrun ./BSSN_GR/bssnSolverCUDA q1_r2.2.par.json 1
