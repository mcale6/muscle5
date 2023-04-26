#!/bin/bash
#SBATCH --job-name=muscle_MSA
#SBATCH --output=msa_%j.out
#SBATCH --error=msa_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Load necessary modules
module load Python/3.8.6-GCCcore-10.2.0

# Change to the directory containing the Python script and input file
cd /home/adaddi/data/muscle5/

# Run the Python script with the specified arguments
python muscle_MSA.py /home/adaddi/scratch/muscle_resampling/c_crbn_mmseqs2_uniref_env_org.a3m c_crbn_uniref_env diversified
