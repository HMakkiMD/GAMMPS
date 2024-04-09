#!/bin/bash -l

#SBATCH --job-name "PolymerName"
#SBATCH -p "CPU_node"
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time 1-00:00:00

module load apps/gromacs_cuda/2022.0
module load apps/anaconda3/2022.10
module load apps/vmd/1.9.4a38
source activate myenv

python3 OligomerBuilder.py
while qstat -u hmakki | grep opt; do echo "I am waiting"; sleep 60; done
python3 ChargeTorsionInput.py
while qstat -u hmakki | grep -tor; do echo "I am waiting"; sleep 3600; done
python3 Oeff.py
while qstat -u hmakki | grep corr; do echo "I am waiting"; sleep 30; done
python3 TorsionCorrection.py
while qstat -u hmakki | grep corr; do echo "I am waiting"; sleep 30; done
python3 FF_test.py
