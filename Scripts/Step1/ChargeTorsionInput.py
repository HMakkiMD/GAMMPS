#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import subprocess
import copy as cp
import time
import re
start = time.time()

from parameters import *


'''
Extracting the single point energy from Gaussian log file. 
Just give the path to open the log file that contains optimization by Gaussian 
and the code will give "ENERGYSUB", which is a list contain the exess value of energy 
for each configuration as compared to the configuration with the lowest energy. Thus,
the index of the min of ENERGYSUB is the most stable configuration.
'''

parent_dir=cp.deepcopy(PATH_OUTPUT) #the Gaussian log files should be generated in the same folder as oligomer.xyz and .com files
moleculename = OLIGOMERNAME
path_output = os.path.join(parent_dir, moleculename)
path_fragments=cp.deepcopy(PATH_FRAGMENTS)
theory=cp.deepcopy(THEORY)

fragmentlist=cp.deepcopy(FRAGMENTLIST)
fragmentlist.insert(0,fragmentlist[-1]) # add capping fragments (to the first end) to make oligomer
fragmentlist.insert(len(fragmentlist),fragmentlist[1]) # add capping fragments (to the second end) to make oligomer

file_num_list=[]

word = 'SCF Done'
ENERGYLAST=[]
for olig_num in range(int(2**(len(fragmentlist)-1))):
    with open(path_output+'/oligomer-'+str(olig_num+1)+'.log', 'r') as fp:
        # read all lines in a 
        lines = fp.readlines()
        ENERGY=[]
        for line in lines:
            # check if string present on a current line
            if (line.find(word) != -1):
                energy=line.split()
                ENERGY.append(float(energy[4])*27.2114) #Hartree to ev
        ENERGYLAST.append(ENERGY[-1])

ENERGYSUB=[x-min(ENERGYLAST) for x in ENERGYLAST]

with open(path_output+'.log', "a") as f:
    f.write("\n\nStep 1-2: Calculating point charges and torsional potentials for oligomer with minimum energy level \n\n")
    f.write("\nMaximum energy difference between different isomers in eV is: "+str(round(max(ENERGYSUB),2))+"\n")
    f.write("\nMinimum energy difference between different isomers in eV is: "+str(round(sorted(ENERGYSUB)[1],2))+"\n")
    f.write("\nOligomer number with minimum energy is: "+ str(ENERGYSUB.index(min(ENERGYSUB))+1)+"\n")
    f.write("\nOligomer number with maximum energy is: "+ str(ENERGYSUB.index(max(ENERGYSUB))+1)+"\n")
with open(path_output+'.html', "a") as f:
    f.write('\n\n<p><font size = "+1">\nCalculating point charges and torsional potentials for oligomer with minimum energy level </font> </p>\n\n')
    f.write("\nMaximum energy difference between different isomers in eV is: "+str(round(max(ENERGYSUB),2))+"<br>\n")
    f.write("\nMinimum energy difference between different isomers in eV is: "+str(round(sorted(ENERGYSUB)[1],2))+"<br>\n")
    f.write("\nOligomer number with minimum energy is: "+ str(ENERGYSUB.index(min(ENERGYSUB))+1)+"<br>\n")
    f.write("\nOligomer number with maximum energy is: "+ str(ENERGYSUB.index(max(ENERGYSUB))+1)+"<br><br>\n")

'''

Making input file for point charge and torsional potential calculations.

'''

# input for point charge calculation
subprocess.run(['obabel','oligomer-'+str(ENERGYSUB.index(min(ENERGYSUB))+1)+'.log','-O', 'oligomer-opt-'+str(ENERGYSUB.index(min(ENERGYSUB))+1)+'.xyz'], cwd = path_output)

with open(path_output+'/oligomer-opt-'+str(ENERGYSUB.index(min(ENERGYSUB))+1)+'.xyz', "r") as f:
        coordlines=f.readlines()[2:]
with open(path_output+'/'+moleculename+'-charge.com', "w") as f:
        f.write("%nproc=10\n")
        f.write("%mem=9GB\n")
        f.write("%chk="+moleculename+"-chelpg.chk\n")
        f.write("#p opt "+theory[1]+" pop=chelpg\n\n")
        f.write("Gaussian optimization calculation on oligomer"+str(ENERGYSUB.index(min(ENERGYSUB))+1)+"\n\n")
        f.write("0,1\n")
        f.writelines(coordlines)
        f.write("\n")

fragmentlist=cp.deepcopy(FRAGMENTLIST)
fragmentlist.insert(0,fragmentlist[-1]) # add capping fragments (to the first end) to make oligomer
fragmentlist.insert(len(fragmentlist),fragmentlist[1]) # add capping fragments (to the second end) to make oligomer

# finding torsion pairs 

fragment_pairs=[]
fragment_pairs_index=[]
for i in range(len(fragmentlist)-1):
    fragment_pairs_index.append([i,i+1])


torsion_pairs=[]
for i in range(len(fragment_pairs_index)):
    torsion_pairs.append(fragment_pairs_index[i][0])
    torsion_pairs.append(fragment_pairs_index[i][0]+1)
# torsion pairs are found --> e.g., [0,1,3,4] means two torsions between 0,1 and 3,4

# finding atom number of torsions for the oligomer
# ****** fragments are assumed to be symmetric so that IDT-BT and BT-IDT have a similar torsional potential ******
fragmentlist_natom=[]
for i,each in enumerate(fragmentlist):
    with open(path_fragments+each+'.xyz') as f:
        if (i==0):
            fragmentlist_natom.append(int(f.readlines()[0])-1)
        else:
            fragmentlist_natom.append(int(f.readlines()[0])-2)
fragmentlist_cumsum_second=np.cumsum(fragmentlist_natom)
fragmentlist_cumsum_first=np.insert(fragmentlist_cumsum_second[:-1], 0, 0)
torsion_atoms=[]
for i in range(0,len(torsion_pairs),2):
    torsion_atoms.append([fragmentlist_cumsum_first[torsion_pairs[i]]+2,fragmentlist_cumsum_first[torsion_pairs[i]]+1,fragmentlist_cumsum_second[torsion_pairs[i+1]],fragmentlist_cumsum_second[torsion_pairs[i+1]]-1])

subprocess.run(['obabel','oligomer-opt-'+str(ENERGYSUB.index(min(ENERGYSUB))+1)+'.xyz','-O', 'oligomer-opt-'+str(ENERGYSUB.index(min(ENERGYSUB))+1)+'.gzmat'], cwd = path_output)

with open(path_output+'/oligomer-opt-'+str(ENERGYSUB.index(min(ENERGYSUB))+1)+'.gzmat', "r") as f:
        gzmat=f.readlines()[6:]
for i in range(len(torsion_atoms)):
    with open(path_output+'/'+moleculename+'-torsion-'+str(i+1)+'.gjf', "w") as f:
            f.write("%nproc=10\n")
            f.write("%chk="+moleculename+'-torsion-'+str(i+1)+'.chk\n')
            f.write("#"+theory[1]+" SCF=tight Opt=ModRedundant\n\n")
            f.write("Torsion calculation on oligomer"+str(ENERGYSUB.index(min(ENERGYSUB))+1)+"\n\n")
            f.write("0,1\n")
            f.writelines(gzmat)
            f.write(str(torsion_atoms[i][0])+' '+str(torsion_atoms[i][1])+' '+str(torsion_atoms[i][2])+' '+str(torsion_atoms[i][3])+' S 36 10.0'+"\n\n")
with open(path_output+'.log', "a") as f:
    f.write("\nNumber of torsion paris for DFT scan: "+str(len(torsion_atoms))+"\n")
    for i in range(0,len(torsion_pairs),2):
        f.write("\nTorsion pair No. "+str(int(i/2)+1)+": "+fragmentlist[torsion_pairs[i]]+"-"+fragmentlist[torsion_pairs[i+1]]+"\n")
        f.write("\nTorsion atoms of this torsion "+str(torsion_atoms[int(i/2)][0])+' '+str(torsion_atoms[int(i/2)][1])+' '+str(torsion_atoms[int(i/2)][2])+' '+str(torsion_atoms[int(i/2)][3])+"\n")
    end = time.time()
    f.write("\n\nRuntime of this step was "+str(round(end - start,1))+" s\n")
with open(path_output+'.html', "a") as f:
    f.write("\nNumber of torsion paris for DFT scan: "+str(len(torsion_atoms))+"<br>\n")
    for i in range(0,len(torsion_pairs),2):
        f.write("\nTorsion pair No. "+str(int(i/2)+1)+": "+fragmentlist[torsion_pairs[i]]+"-"+fragmentlist[torsion_pairs[i+1]]+"<br>\n")
        f.write("\nTorsion atoms of this torsion "+str(torsion_atoms[int(i/2)][0])+' '+str(torsion_atoms[int(i/2)][1])+' '+str(torsion_atoms[int(i/2)][2])+' '+str(torsion_atoms[int(i/2)][3])+"<br>\n")
    end = time.time()
    f.write("\n\n<b>The level of theory used for atomic charge and torsional potential calculation is "+theory[1]+" <br><br>\n")
    f.write("\n\nRuntime of this step was "+str(round(end - start,1))+" s<br><br>\n")
with open(path_output+'/run-charge.sh', 'w') as f:
    f.write("#!/bin/bash -l\n\n#SBATCH --job-name "+moleculename+"-charge\n#SBATCH -p "+CPUNODE+"\n#SBATCH -N 1\n#SBATCH -n 10\n#SBATCH --time 0-24:00:00\n\n")
    f.write("module load apps/gaussian/16\n\n")
    f.write("g16 "+moleculename+"-charge.com\n\n")
if (SUBMITJOB=='yes'):
    subprocess.run(['sbatch','run-charge.sh'], cwd = path_output)

with open(path_output+'/edit-torsion.sh', 'w') as f:
    f.write("#!/bin/bash\n\n")
    f.write("for a in {1.."+str(len(torsion_atoms))+"}\n")
    f.write('do\n  sed -i "s/torsion-[0-9]*/torsion-${a}/g" run-torsion.sh\n  sbatch run-torsion.sh\ndone\n\n')

with open(path_output+'/run-torsion.sh', 'w') as f:
    f.write("#!/bin/bash -l\n\n#SBATCH --job-name "+moleculename+"-tor1\n#SBATCH -p "+CPUNODE+"\n#SBATCH -N 1\n#SBATCH -n 10\n#SBATCH --time 3-00:00:00\n\n")
    f.write("module load apps/gaussian/16\n\n")
    f.write("g16 "+moleculename+"-torsion-1.gjf\n\n")
if (SUBMITJOB=='yes'):
    subprocess.run(['chmod','+x','edit-torsion.sh'], cwd = path_output)
    subprocess.run(['./edit-torsion.sh'], cwd = path_output)

