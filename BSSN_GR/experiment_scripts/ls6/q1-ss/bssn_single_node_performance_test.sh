#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Lonestar6 AMD Milan nodes
#
#   *** Serial Job in Normal Queue***
# 
# Last revised: October 22, 2021
#
# Notes:
#
#  -- Copy/edit this script as desired.  Launch by executing
#     "sbatch milan.serial.slurm" on a Lonestar6 login node.
#
#  -- Serial codes run on a single node (upper case N = 1).
#       A serial code ignores the value of lower case n,
#       but slurm needs a plausible value to schedule the job.
#
#  -- Use TACC's launcher utility to run multiple serial 
#       executables at the same time, execute "module load launcher" 
#       followed by "module help launcher".
#----------------------------------------------------

#SBATCH -J bssn           # Job name
#SBATCH -o .single_node.o%j       # Name of stdout output file
#SBATCH -e .single_node.e%j       # Name of stderr error file
#SBATCH -p gpu-a100          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=all    # Send email at begin and end of job
##SBATCH -A st       # Project/Allocation name (req'd if you have more than 1)
##SBATCH --mail-user=username@tacc.utexas.edu

# Any other commands must follow all #SBATCH directives...
module list
pwd
date

# Launch serial code...
#./myprogram         # Do not use ibrun or any other MPI launcher
echo "unzip-scatter vs unzip-gather"
export OMP_NUM_THREADS=1
./run_meshgpu_tests 14 1e-4 0.1 6 0
./run_meshgpu_tests 14 1e-5 0.1 6 0
./run_meshgpu_tests 14 1e-6 0.1 6 0
./run_meshgpu_tests 14 1e-7 0.1 6 0
./run_meshgpu_tests 14 1e-8 0.1 6 0


export OMP_NUM_THREADS=128
echo "sympy_cse"
./BSSN_GR/rhsgpu_test_algebraic_sympy_cse 256 10 0
./BSSN_GR/rhsgpu_test_algebraic_sympy_cse 512 10 0
./BSSN_GR/rhsgpu_test_algebraic_sympy_cse 1024 10 0
./BSSN_GR/rhsgpu_test_algebraic_sympy_cse 2048 10 0
./BSSN_GR/rhsgpu_test_algebraic_sympy_cse 4096 10 0
./BSSN_GR/rhsgpu_test_algebraic_sympy_cse 8192 10 0
./BSSN_GR/rhsgpu_test_algebraic_sympy_cse 16384 10 0


echo "nx_cse"
./BSSN_GR/rhsgpu_test_algebraic_nx_cse 256 10 0
./BSSN_GR/rhsgpu_test_algebraic_nx_cse 512 10 0
./BSSN_GR/rhsgpu_test_algebraic_nx_cse 1024 10 0
./BSSN_GR/rhsgpu_test_algebraic_nx_cse 2048 10 0
./BSSN_GR/rhsgpu_test_algebraic_nx_cse 4096 10 0
./BSSN_GR/rhsgpu_test_algebraic_nx_cse 8192 10 0
./BSSN_GR/rhsgpu_test_algebraic_nx_cse 16384 10 0

echo "fused_with_dorder"
./BSSN_GR/rhsgpu_test_fused_key_expr_sympy_cse  256 10 0
./BSSN_GR/rhsgpu_test_fused_key_expr_sympy_cse  512 10 0
./BSSN_GR/rhsgpu_test_fused_key_expr_sympy_cse 1024 10 0
./BSSN_GR/rhsgpu_test_fused_key_expr_sympy_cse 2048 10 0
./BSSN_GR/rhsgpu_test_fused_key_expr_sympy_cse 4096 10 0
./BSSN_GR/rhsgpu_test_fused_key_expr_sympy_cse 8192 10 0
./BSSN_GR/rhsgpu_test_fused_key_expr_sympy_cse 16384 10 0




