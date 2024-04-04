#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 72:00:00
#SBATCH -J TIF-BT
#SBATCH -p "CPUnode"

module load apps/gaussian/16
module load apps/anaconda3/2022.10
source activate myenv

python3 QC_calculation.py input_variables.inp > output_QC.txt
