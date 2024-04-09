#!/bin/bash -l

#SBATCH --job-name "PolymerName"
#SBATCH -p "CPUnode"
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --time 1-00:00:00

module load apps/anaconda3/2022.10
module load apps/vmd/1.9.4a38
module load apps/gromacs_cuda/2022.0
source activate myenv

python3 RUbuilder.py
python3 SCbuilder.py
python3 PObuilder.py
