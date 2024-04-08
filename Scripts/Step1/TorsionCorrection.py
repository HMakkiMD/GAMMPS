#!/usr/bin/env python
# coding: utf-8

'''
This is the fourth code to run in Step 1.

This code performs the correction for intermonomer torsional potential. 
It calculates the ΔV = V_QM − V_MM where V_QM represents the total potential energy of the oligomer at each scan point from DFT scan
and V_MM is the total potential energy minus the targeted torsional potential obtained from the generated force field at each scan point.
For each scan point, the code generates input files for energy minimization by GROMACS with restraints on the targeted torsion.
The magnitude of the restraint should be provided in parameters.py.

'''

import subprocess
import numpy as np
from matplotlib import pyplot as plt
import os
import time
import re
import copy as cp
from scipy.interpolate import CubicSpline
import collections
start = time.time()

def dihedral_calculator(point_1, point_2, point_3, point_4):
    a = np.array(point_1)
    b = np.array(point_2)
    c = np.array(point_3)
    d = np.array(point_4)
    b0 = -1.0*(b - a)
    b1 = c - b
    b2 = d - c

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

from parameters import *

parent_dir=cp.deepcopy(PATH_OUTPUT) #the Gaussian log files should be generated in the same folder as oligomer.xyz and .com files
moleculename = OLIGOMERNAME
path_output = os.path.join(parent_dir, moleculename)
path_output_torsion = os.path.join(parent_dir, moleculename+'/Torsion/') # where the charge and torsion files will be find
path_output_oeff=os.path.join(parent_dir, moleculename+'/') # where the force field files will be coppied
fragmentlist=cp.deepcopy(FRAGMENTLIST)
fragmentlist.insert(0,fragmentlist[-1]) # add capping fragments (to the first end) to make oligomer
fragmentlist.insert(len(fragmentlist),fragmentlist[1]) # add capping fragments (to the second end) to make oligomer
path_fragments=cp.deepcopy(PATH_FRAGMENTS)
col=cp.deepcopy(COLOR)


file_num_list=[]

for i in range(len(fragmentlist)-1):
    file_num_list.append(i+1)

LEGEND=[]
for i in file_num_list:
	LEGEND.append('torsion'+str(i))
	LEGEND.append('torsion'+str(i)+'-table')

torsion_pairs=[]
for i in range(len(fragmentlist)-1):
	torsion_pairs.append(i)
	torsion_pairs.append(i+1)

fragmentlist_natom=[]
for i,each in enumerate(fragmentlist):
    with open(path_fragments+each+'.xyz') as f:
        if (i==0):
            fragmentlist_natom.append(int(f.readlines()[0])-1)
        else:
            fragmentlist_natom.append(int(f.readlines()[0])-2)
fragmentlist_cumsum_second=np.cumsum(fragmentlist_natom)
fragmentlist_cumsum_first=np.insert(fragmentlist_cumsum_second[:-1], 0, 0)
all_torsion_atoms=[]
for i in range(0,len(torsion_pairs),2):
    all_torsion_atoms.append([fragmentlist_cumsum_first[torsion_pairs[i]]+2,fragmentlist_cumsum_first[torsion_pairs[i]]+1,fragmentlist_cumsum_second[torsion_pairs[i+1]],fragmentlist_cumsum_second[torsion_pairs[i+1]]-1])



TORSION_ATOMS=[]
DIHEDRAL_ANGLES=[]
SORTED_ENERGIES=[]

for i_index, i in enumerate(file_num_list):
	with open(path_output_oeff+moleculename+"-torsion-"+str(i)+".gjf", "r") as f:
		lines = (line.rstrip() for line in f) # All lines including the blank ones
		lines = list(line for line in lines if line) # Non-blank lines
		torsion_atoms=[int(lines[-1].split()[0]), int(lines[-1].split()[1]), int(lines[-1].split()[2]), int(lines[-1].split()[3])]
	TORSION_ATOMS.append(torsion_atoms)
	dihedral_angles=[]
	ang_rng=[0,36]
	#increment
	incr=1
	dhd_rngs=range(ang_rng[0],ang_rng[1]+1,incr)
	for j in dhd_rngs:
		with open(path_output_torsion+moleculename+"-torsion-"+str(i)+'-'+str(j)+'.xyz', 'r') as f:
			lines=f.readlines()
			coord1=[float(lines[torsion_atoms[0]+1].split()[1]),float(lines[torsion_atoms[0]+1].split()[2]),float(lines[torsion_atoms[0]+1].split()[3])]
			coord2=[float(lines[torsion_atoms[1]+1].split()[1]),float(lines[torsion_atoms[1]+1].split()[2]),float(lines[torsion_atoms[1]+1].split()[3])]
			coord3=[float(lines[torsion_atoms[2]+1].split()[1]),float(lines[torsion_atoms[2]+1].split()[2]),float(lines[torsion_atoms[2]+1].split()[3])]
			coord4=[float(lines[torsion_atoms[3]+1].split()[1]),float(lines[torsion_atoms[3]+1].split()[2]),float(lines[torsion_atoms[3]+1].split()[3])]
		dihedral_angles.append(dihedral_calculator(coord1,coord2,coord3,coord4))
	DIHEDRAL_ANGLES.append(dihedral_angles)
	with open(path_output_torsion+'/'+moleculename+'-'+str(i)+'-FF/energy.txt', "r") as f:
		lines=f.readlines()
		energies=[]
		for line in lines:
			energies.append(float(line.split()[9])-float(line.split()[8])) #potential energy minus dihedral restraint
	energies=np.array(energies)+abs(min(energies))
	
	
	
	dihedral_angles_sorted = np.array(dihedral_angles[:-1]).argsort()
	sorted_dihedral_angles = np.array(dihedral_angles[:-1])[dihedral_angles_sorted]
	sorted_energies = np.array(energies[:-1])[dihedral_angles_sorted]
	sorted_energies_norm=sorted_energies-min(sorted_energies)

	#select the positive or negative (or both) angles are good for fitting
	if HALF_FITTING=='yes':
		if sorted_dihedral_angles[list(sorted_energies_norm).index(min(sorted_energies_norm))]>=0:
			for each in range(int(len(sorted_energies_norm)/2)):
				sorted_energies_norm[each]=sorted_energies_norm[35-each]
		else:
			for each in range(18,len(sorted_energies_norm)):
				sorted_energies_norm[each]=sorted_energies_norm[35-each]

	sorted_dihedral_angles_extended=np.insert(sorted_dihedral_angles,0,min(sorted_dihedral_angles)-10, axis=0)
	sorted_dihedral_angles_extended=np.append(sorted_dihedral_angles_extended,max(sorted_dihedral_angles)+10)
	sorted_energies_extended=np.insert(sorted_energies_norm,0,sorted_energies_norm[-1], axis=0)
	sorted_energies_extended=np.append(sorted_energies_extended,sorted_energies_norm[0])

	cs = CubicSpline(sorted_dihedral_angles_extended, sorted_energies_extended)
	xs=np.linspace(-180., 180., 361)
	#xs=np.linspace(-540., 540., 1081)
	#print(sorted_dihedral_angles_extended)
	
	
	SORTED_ENERGIES.append(sorted_energies_extended)
	#col = (np.random.random(), np.random.random(), np.random.random())
	plt.figure(1)
	plt.plot(sorted_dihedral_angles,sorted_energies_norm, c=col[i], linestyle='none', marker='D')
	plt.plot(xs,cs(xs), c=col[i], linestyle='dashed')


plt.xlabel("dihedral angle [degree]", fontsize=14) 
plt.ylabel("energy [kJ/mol]", fontsize=14)
lgd=plt.legend(LEGEND,
	           bbox_to_anchor=(1, 0.5), prop={'size': 12})

plt.savefig(path_output_torsion+'FF.jpeg', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')

with open(path_output+'.html', "a") as f:
    f.write('<aside class="figures">\n')
    f.write('<figure>\n<img src="'+moleculename+'/Torsion/FF.jpeg" style="width: 500px; height: 300px; object-fit: cover;" alt="V<sup>FF</sup>">\n<figcaption>V<sup>FF</sup></figcaption>\n</figure>')
    f.write('</aside>\n')

for enr_index, energy in enumerate(SORTED_ENERGIES):
	str_tors=subprocess.check_output(["sed", "-n", "/N-N= /,/@/p", path_output_torsion+moleculename+"-torsion-"+str(enr_index+1)+".log"])
	newl_tors=str_tors.decode('utf-8').replace("\n","").replace(" ","")
	tors_no_suffix_prefix=newl_tors.split("HF=")[1].split("RMSD=")[0].rstrip("\\").rstrip("/")
	energies_tors=tors_no_suffix_prefix.split(",")
	energies_tors_eV=np.array(energies_tors, float)*27.2114 # unit in eV 
	energies_tors_kJmol=np.array(energies_tors, float)*2625.5 # unit in kJ/mol
	energies_tors_eV_norm=energies_tors_eV-min(energies_tors_eV)
	energies_tors_kJmol_norm=energies_tors_kJmol-min(energies_tors_kJmol)
	dihedral_angles_sorted = np.array(DIHEDRAL_ANGLES[enr_index][:-1]).argsort()
	sorted_dihedral_angles = np.array(DIHEDRAL_ANGLES[enr_index][:-1])[dihedral_angles_sorted]
	sorted_energies_tors_kJmol_norm = np.array(energies_tors_kJmol_norm[:-1])[dihedral_angles_sorted]
	
	#select the positive or negative (or both) angles are good for fitting
	if HALF_FITTING=='yes':
		if sorted_dihedral_angles[list(sorted_energies_tors_kJmol_norm).index(min(sorted_energies_tors_kJmol_norm))]>=0:
			for each in range(int(len(sorted_energies_tors_kJmol_norm)/2)):
				sorted_energies_tors_kJmol_norm[each]=sorted_energies_tors_kJmol_norm[35-each]
		else:
			for each in range(18,len(sorted_energies_tors_kJmol_norm)):
				sorted_energies_tors_kJmol_norm[each]=sorted_energies_tors_kJmol_norm[35-each]

	sorted_dihedral_angles_extended=np.insert(sorted_dihedral_angles,0,min(sorted_dihedral_angles)-10, axis=0)
	sorted_dihedral_angles_extended=np.append(sorted_dihedral_angles_extended,max(sorted_dihedral_angles)+10)
	sorted_energies_tors_kJmol_norm_extended=np.insert(sorted_energies_tors_kJmol_norm,0,sorted_energies_tors_kJmol_norm[-1], axis=0)
	sorted_energies_tors_kJmol_norm_extended=np.append(sorted_energies_tors_kJmol_norm_extended,sorted_energies_tors_kJmol_norm[0])
	v_corr=np.array(sorted_energies_tors_kJmol_norm_extended)-np.array(energy)
	v_corr=np.array(v_corr)+abs(min(v_corr))

	
	
	cs = CubicSpline(sorted_dihedral_angles_extended, v_corr)
	xs=np.linspace(-180., 180., 361)
	plt.figure(2)
	plt.plot(sorted_dihedral_angles_extended,v_corr, c=col[enr_index+1], linestyle='none', marker='D')
	plt.plot(xs,cs(xs), c=col[enr_index+1], linestyle='dashed')
	table=[]
	for j in range(len(xs)):
		table.append(str(xs[j])+'    '+str(round(cs(xs)[j],3))+'     '+str(-round(cs(xs, 1)[j],3)))
	with open(path_output_torsion+'table_corr_d'+str(enr_index+1)+'.xvg', 'w') as f:
		f.writelines("%s\n" % l for l in table)
	gromacs_command='gmx' # gmx command based on the module loaded on Barkla
	minim_mdp='minim.mdp' # name of the minimization mdp file
	box_size='10' # box size for minimization of 1 monomer
	energy_list='1\\n2\\n3\\n4\\n5\\n6\\n7\\n8\\n9\\n10\\n' # enery contributions to be calculated by the energy command of the gromcas

    # writing table for corrected potential
	
	with open(path_output_torsion+'/'+moleculename+'-'+str(enr_index+1)+'-FF/'+'run_corr.sh', 'w') as f:
		f.write("#!/bin/bash -l\n\n")
		f.write("#SBATCH --job-name corr_chk_"+str(enr_index+1)+"\n#SBATCH -p "+GPUNODE+"\n#SBATCH -N 1\n#SBATCH -n 6\n#SBATCH --gres gpu:1\n#SBATCH --time 1-00:00:00\n\n")
		f.write("module load apps/gromacs_cuda/2022.0\n\n")
		f.write("module load apps/openbabel/2.4.1/gcc-4.8.5+eigen-3.0.5\n")
		f.write("export GMX_MAXBACKUP=-1\n\n")
		f.write("for i in {0..36}; do \n")
		f.write(gromacs_command+" editconf -f $i.gro -box "+box_size+" -o $i.gro \n\n")
		f.write(gromacs_command+" grompp -f ../"+minim_mdp+" -c $i.gro -p $i-tab.top -o $i.tpr \n")
		f.write(gromacs_command+" mdrun -deffnm $i -nt 6 -tableb ../table_corr_d"+str(enr_index+1)+".xvg\n")
		f.write("echo -e '"+energy_list+"' | "+gromacs_command+" energy -f $i.edr -o $i.xvg\n")
		f.write("done \n\n")
		f.write("grep '[0-9]' 0.xvg | tail -1 > energy_corr.txt \n\n")
		f.write("for i in {1..36}; do \n")
		f.write("grep '[0-9]' $i.xvg | tail -1 >> energy_corr.txt \n")
		f.write("done \n")
		f.write("rm *.trr *.tpr *.cpt *.edr *.log")
	if (SUBMITJOB=='yes'):
		subprocess.run(['sbatch','run_corr.sh'],cwd = path_output_torsion+'/'+moleculename+'-'+str(enr_index+1)+'-FF/')

plt.xlabel("dihedral angle [degree]", fontsize=14) 
plt.ylabel("energy [kJ/mol]", fontsize=14)
lgd=plt.legend(LEGEND,
	           bbox_to_anchor=(1, 0.5), prop={'size': 12})

plt.savefig(path_output_torsion+'CORR.jpeg', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')

end = time.time()

with open(path_output+'.log', "a") as f:
    f.write("\nCorrected torsional potential parameters were calculated \n")
    f.write("\nMD scan over 36 frames with corrected torsional potential have been submitted\n\nRuntime of this step was "+str(round(end - start,1))+" s\n")


with open(path_output+'.html', "a") as f:
	f.write("\nCorrected torsional potential parameters were calculated <br><br>\n")
	f.write("\nMD scan over 36 frames with corrected torsional potential have been submitted\n\nRuntime of this step was "+str(round(end - start,1))+" s<br><br>\n")
	f.write('<aside class="figures">\n')
	f.write('<figure>\n<img src="'+moleculename+'/Torsion/CORR.jpeg" style="width: 500px; height: 300px; object-fit: cover;" alt="V<sup>CORR</sup>">\n<figcaption>V<sup>CORR</sup></figcaption>\n</figure>')
	f.write('</aside>\n')
	f.write("\nMD scan over 36 frames with corrected torsional potential have been submitted\n\n<br><br>Runtime of this step was "+str(round(end - start,1))+" s<br><br>\n")
	f.write("\nV<sup>DFT</sup> vs V<sup>TOTAL</sup> are shown below <b><b>")



