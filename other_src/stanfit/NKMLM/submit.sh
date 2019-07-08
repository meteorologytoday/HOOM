#!/bin/bash
#
#SBATCH --job-name=SMARTSLAB-fit-LENS
#SBATCH --output=slurm.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-user=tienyiah@uci.edu
#
#
#SBATCH --array=1-100%80

julia $( dirname $0 )/single_x.jl --selected-index=1

