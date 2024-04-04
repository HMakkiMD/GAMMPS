#!/bin/bash -l

#SBATCH --job-name DPPBTZ-charge
#SBATCH -p troisi
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time 0-24:00:00

module load apps/gaussian/16

g16 DPPBTZ-charge.com

