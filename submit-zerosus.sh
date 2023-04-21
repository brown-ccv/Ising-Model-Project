#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH -c 48
#SBATCH --exclusive
#SBATCH --partition=batch
#SBATCH --mem=20G
#SBATCH --job-name ising-sim-threaded
#SBATCH --output logs/threaded/ising-sim-%A-%a.out
#SBATCH --error logs/threaded/ising-sim-%A-%a.out
#SBATCH --constraint cascade
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cpaniaguam@brown.edu

module load julia/1.8.5

# kT=`head -n $SLURM_ARRAY_TASK_ID temperatures.txt | tail -1`
julia --threads $SLURM_CPUS_PER_TASK --project=. zerosus.jl --nspins 50 --mcsteps 1000000 --plt true

echo "Done!"
