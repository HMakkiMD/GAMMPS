#!/bin/bash -l

#SBATCH --job-name DPPBTZ-tor1
#SBATCH -p troisi
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time 3-00:00:00

module load apps/gaussian/16

g16 DPPBTZ-torsion-1.gjf

