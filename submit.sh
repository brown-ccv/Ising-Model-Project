#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=01:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --partition=batch
#SBATCH --mem=60G
#SBATCH --job-name ising-sim 
#SBATCH --output logs/ising-sim-%j.out
#SBATCH --error logs/ising-sim-%j.out
#SBATCH --constraint cascade

module load julia/1.8.5

julia --threads 48 --project=. src/zerofield_susceptibility.jl

echo "Done!"
