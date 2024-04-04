#!/bin/bash -l

#SBATCH --job-name TIF-BT
#SBATCH -p "CPUnode"
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time 3-00:00:00

module load apps/anaconda3/2022.10
module load apps/vmd/1.9.4a38
module load apps/gromacs_cuda/2022.0
source activate myenv

python3 DOSindex.py
python3 DOSinput.py
python3 DOSinput_soup.py
