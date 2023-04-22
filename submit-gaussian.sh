#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=80:00:00
#SBATCH --nodes=1
#SBATCH -c 24
#SBATCH --partition=batch
#SBATCH --mem=5G
#SBATCH --job-name ising-sim-threaded
#SBATCH --output logs/threaded/ising-sim-%A-%a.out
#SBATCH --error logs/threaded/ising-sim-%A-%a.out
#SBATCH --constraint cascade
#SBATCH --mail-type=ALL
#SBATCH --array=1-10

module load julia/1.8.5

#specify the standard deviation of gaussian field by setting sd to desired value below:

julia --threads $SLURM_CPUS_PER_TASK --project=. finalcodes/gaussianRF_sim.jl --nspins 300 --mcsteps 1_000_000 --sd 0.8

echo "Done!"
