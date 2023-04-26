#!/bin/bash
#SBATCH --job-name=muscle_MSA
#SBATCH --output=msa_%j.out
#SBATCH --error=msa_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G

module load mamba intel
source activate main_env

A3M_FILE=$1
NAME=$2
MODE=$3

python /home/adaddi/data/muscle5/muscle_MSA.py --a3m_file "$A3M_FILE" --name "$NAME" --mode "$MODE"
