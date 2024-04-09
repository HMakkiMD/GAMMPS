#!/bin/bash -l

for a in {0..8}; do 
sed -i "s/oligomer-[0-9]*/oligomer-${a}/g" run.sh
sed -i "s/opt-[0-9]*/opt-${a}/g" run.sh
sbatch run.sh
done 

