#!/bin/bash
#
#SBATCH --job-name=SMARTSLAB-fit-LENS
#SBATCH --output=slurm.log
#SBATCH --partition=sib2.9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-user=tienyiah@uci.edu
#
#
#SBATCH --array=1-100%60

srun julia $SLURM_SUBMIT_DIR/single_x.jl --selected-index=$SLURM_ARRAY_TASK_ID

