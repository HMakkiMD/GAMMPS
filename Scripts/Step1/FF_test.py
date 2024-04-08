#!/usr/bin/env python
# coding: utf-8

'''
This is the fifth code to run in Step 1.

This code checks if the total energy from the force field at each scan point matches the DFT-calculated torsional potential.
It plots both values (V_tot and V_DFT) at each scan point for all torsions in the oligomer structure separately which can be seen
in the html log file generated alongside the output files.

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
SORTED_DIHEDRAL_ANGLES=[]
SORTED_ENERGIES=[]
for i in file_num_list:
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
	
	with open(path_output_torsion+'/'+moleculename+'-'+str(i)+'-FF/energy_corr.txt', "r") as f:
		lines=f.readlines()
		energies_total=[]
		for line in lines:
			energies_total.append(float(line.split()[10])-float(line.split()[9])) #potential energy minus dihedral restraint
	#energies_total=np.array(energies_total)+abs(min(energies_total))
	energies_total_norm=np.array(energies_total[:-1])-min(energies_total[:-1])


	str_tors=subprocess.check_output(["sed", "-n", "/N-N= /,/@/p", path_output_torsion+moleculename+"-torsion-"+str(i)+".log"])
	newl_tors=str_tors.decode('utf-8').replace("\n","").replace(" ","")
	tors_no_suffix_prefix=newl_tors.split("HF=")[1].split("RMSD=")[0].rstrip("\\").rstrip("/")
	energies_tors=tors_no_suffix_prefix.split(",")
	energies_tors_eV=np.array(energies_tors, float)*27.2114 # unit in eV 
	energies_tors_kJmol=np.array(energies_tors, float)*2625.5 # unit in kJ/mol
	energies_tors_eV_norm=energies_tors_eV-min(energies_tors_eV)
	energies_tors_kJmol_norm=energies_tors_kJmol-min(energies_tors_kJmol)
	

	dihedral_angles_sorted = np.array(dihedral_angles[:-1]).argsort()
	sorted_dihedral_angles = np.array(dihedral_angles[:-1])[dihedral_angles_sorted]
	sorted_energies_tors_kJmol_norm = np.array(energies_tors_kJmol_norm[:-1])[dihedral_angles_sorted]
	
	if HALF_FITTING=='yes':
		if sorted_dihedral_angles[list(sorted_energies_tors_kJmol_norm).index(min(sorted_energies_tors_kJmol_norm))]>=0:
			for each in range(int(len(sorted_energies_tors_kJmol_norm)/2)):
				sorted_energies_tors_kJmol_norm[each]=sorted_energies_tors_kJmol_norm[35-each]
		else:
			for each in range(18,len(sorted_energies_tors_kJmol_norm)):
				sorted_energies_tors_kJmol_norm[each]=sorted_energies_tors_kJmol_norm[35-each]
	
	
	
	#col = (np.random.random(), np.random.random(), np.random.random())
	plt.figure(i)
	plt.plot(dihedral_angles[:-1],energies_total_norm, c=col[i], linestyle='none', marker='D', markerfacecolor='none')
	plt.plot(sorted_dihedral_angles,sorted_energies_tors_kJmol_norm, c=col[i], linestyle='none', marker='o', markerfacecolor='none')


	plt.xlabel("dihedral angle [degree]", fontsize=14) 
	plt.ylabel("energy [kJ/mol]", fontsize=14)
	plt.xticks(np.arange(-180, 190, step=30))
	lgd=plt.legend(['torsion'+str(i)+'_total','torsion'+str(i)+'-dft',],
		           bbox_to_anchor=(1, 0.5), prop={'size': 12})

	plt.savefig(path_output_torsion+'TOTAL_'+str(i)+'.jpeg', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')


	with open(path_output+'.html', "a") as f:
	    f.write('<aside class="figures">\n')
	    f.write('<figure>\n<img src="'+moleculename+'/Torsion/TOTAL_'+str(i)+'.jpeg" style="width: 500px; height: 300px; object-fit: cover;" alt="V<sup>TOTAL_'+str(i)+'</sup>">\n<figcaption>V<sup>TOTAL</sup> vs V<sup>DFT</sup> for torsion '+str(i)+'</figcaption>\n</figure>')
	    f.write('</aside>\n')

end = time.time()



