#!/bin/bash

for a in {1..3}
do
  sed -i "s/torsion-[0-9]*/torsion-${a}/g" run-torsion.sh
  sed -i "s/dtt-tor[0-9]*/dtt-tor${a}/g" run-torsion.sh
  sbatch run-torsion.sh
done

