#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH -c 24
#SBATCH --partition=batch
#SBATCH --mem=5G
#SBATCH --job-name ising-sim-threaded
#SBATCH --output logs/threaded/ising-sim-%A-%a.out
#SBATCH --error logs/threaded/ising-sim-%A-%a.out
#SBATCH --constraint cascade
#SBATCH --mail-type=ALL

module load julia/1.8.5

julia --threads $SLURM_CPUS_PER_TASK --project=. finalcodes/zerofield_sim.jl --nspins 50 --mcsteps 1000 --plt true

echo "Done!"
