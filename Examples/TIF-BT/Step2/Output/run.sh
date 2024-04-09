#!/bin/bash -l

#SBATCH --job-name opt-8
#SBATCH -p troisi
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time 0-24:00:00
module load apps/gaussian/16

g16 oligomer-8.com

