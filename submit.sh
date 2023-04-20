#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --partition=batch
#SBATCH --mem=3G
#SBATCH --job-name ising-sim 
#SBATCH --output logs/ising-sim-%A-%a.out
#SBATCH --error logs/ising-sim-%A-%a.out
#SBATCH --constraint cascade
#SBATCH --array=1-21

module load julia/1.8.5

kT=`head -n $SLURM_ARRAY_TASK_ID temperatures.txt | tail -1`
julia --threads $SLURM_CPUS_PER_TASK --project=. src/zerofield_susceptibility.jl $kT

echo "Done!"
