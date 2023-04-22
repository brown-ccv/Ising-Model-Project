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
#SBATCH --array=1-4

module load julia/1.8.5

hf=`head -n $SLURM_ARRAY_TASK_ID field_values.txt | tail -1`
julia --threads $SLURM_CPUS_PER_TASK --project=. finalcodes/randomsign_sim.jl --nspins 300 --mcsteps 1_000_000 --hf $hf

echo "Done!"
