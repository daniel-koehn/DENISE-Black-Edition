#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name="DENISE"
# SBATCH --res=HACKATHON_TEAMS
#SBATCH --gres=gpu:p100:1
#SBATCH --nodes=1
#SBATCH --tasks=36
#SBATCH --time=00:01:00
#SBATCH --err=JOB.%j.err
#SBATCH --output=JOB.%j.out
#SBATCH --exclusive

/usr/bin/time mpirun -n 1 ../bin/denise DENISE_marm_OBC.inp FWI_workflow_marmousi.inp

#mpirun -n 1 nvprof -o resultNV.nvprof --cpu-profiling on ../bin/denise DENISE_marm_OBC.inp FWI_workflow_marmousi.inp | tee profiling_CPU.txt

#mpirun -n 1 nvprof --print-summary-per-gpu --profile-from-start on -o resultNVi_%p.nvprof -f --log-file log_%p.txt  ../bin/denise DENISE_marm_OBC.inp FWI_workflow_marmousi.inp

 
