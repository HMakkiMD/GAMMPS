#!/bin/bash -l

#SBATCH --job-name corr_1
#SBATCH -p gputroisi
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --gres gpu:1
#SBATCH --time 1-00:00:00

module load apps/gromacs_cuda/2022.0

module load apps/openbabel/2.4.1/gcc-4.8.5+eigen-3.0.5
export GMX_MAXBACKUP=-1

for i in {0..36}; do 
gmx editconf -f $i.gro -box 10 -o $i.gro 

gmx grompp -f ../minim.mdp -c $i.gro -p $i.top -o $i.tpr 
gmx mdrun -deffnm $i -nt 6
echo -e '1\n2\n3\n4\n5\n6\n7\n8\n9\n' | gmx energy -f $i.edr -o $i.xvg
done 

grep '[0-9]' 0.xvg | tail -1 > energy.txt 

for i in {1..36}; do 
grep '[0-9]' $i.xvg | tail -1 >> energy.txt 
done 
rm *.trr *.tpr *.cpt *.edr *.log