#!/bin/bash
#
#SBATCH --job-name=fit-NKMLM
#SBATCH --output=slurm.log
#SBATCH --partition=nes2.8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-user=tienyiah@uci.edu
#
#
#SBATCH --array=1-100%10

srun julia $SLURM_SUBMIT_DIR/single_x.jl --selected-index=$SLURM_ARRAY_TASK_ID

