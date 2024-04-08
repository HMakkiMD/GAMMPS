#!/usr/bin/env python
# coding: utf-8

'''
This is the first code to run in this step.

This code generates the index files for GROMACS to generate input files for DOS and localization length calculations.
The output of this step is a prerequisite for the second and third codes of Step 3.

'''


import numpy as np
import subprocess
import copy as cp
import subprocess
import os
import time

from parameters import *

parent_dir=cp.deepcopy(PATH_OUTPUT) #the Gaussian log files should be generated in the same folder as oligomer.xyz and .com files
moleculename = OLIGOMERNAME
path_output = os.path.join(parent_dir, moleculename)
path_output_oeff=os.path.join(parent_dir, moleculename+'/') # where the force field files will be coppied
path_output_dos=os.path.join(parent_dir, moleculename+'/DOS/') # where the force field files will be coppied
degree_polymer=cp.deepcopy(DP)

if (SUBMITJOB=='yes'):
    subprocess.run(['mkdir','DOS'],cwd=path_output_oeff)


with open(path_output_oeff+moleculename+'_RU.xyz', 'r') as f:
    backbone=int(f.readlines()[0])
with open(path_output_oeff+moleculename+'_RU_SC.xyz', 'r') as f:
    backandside=int(f.readlines()[0])
sidechain=backandside-backbone

#for melt polymers
for degree in degree_polymer:
    with open(path_output_oeff+moleculename+'_RU_SC_'+degree+'.xyz', 'r') as f:
        total=int(f.readlines()[0])

    charged=[]
    for chain in range(nochainsinbox[degree]):
        for i in range(int(degree)):
            for j in range(backandside):
                if(j//backbone==0):
                    charged.append(chain*total+int((i)*backandside+(j+1)))
        charged.append((chain+1)*total-1)
        charged.append((chain+1)*total)
    with open(path_output_dos+moleculename+'charge_'+degree+'.ndx', 'w') as f:
        f.write('[ backbone_'+degree+']\n')
        for i in range(len(charged)//15):
            f.writelines("%s " % l for l in charged[i*15:(i+1)*15])
            f.write('\n')
        f.writelines("%s " % l for l in charged[(i+1)*15:(i+1)*15+1+len(charged)%15])

#for soups
for degree in degree_polymer:
    if (int(degree)>2):
        with open(path_output_oeff+moleculename+'_RU_SC_'+degree+'.xyz', 'r') as f:
            total_pol=int(f.readlines()[0])

        charged=[]
        for i in range(int(degree)):
            for j in range(backandside):
                if(j//backbone==0):
                    charged.append(int((i)*backandside+(j+1)))
        charged.append(total_pol-1)
        charged.append(total_pol)
        with open(path_output_oeff+moleculename+'_RU_SC_1.xyz', 'r') as f:
            total_soup=int(f.readlines()[0])
        for chain in range(nochainsinbox['1']):
            for i in range(int('1')):
                for j in range(backandside):
                    if(j//backbone==0):
                        charged.append(chain*total_soup+total_pol+int((i)*backandside+(j+1)))
            charged.append(total_pol+(chain+1)*total_soup-1)
            charged.append(total_pol+(chain+1)*total_soup)
        with open(path_output_dos+moleculename+'charge_'+degree+'_1.ndx', 'w') as f:
            f.write('[ backbone_'+degree+']\n')
            for i in range(len(charged)//15):
                f.writelines("%s " % l for l in charged[i*15:(i+1)*15])
                f.write('\n')
            f.writelines("%s " % l for l in charged[(i+1)*15:(i+1)*15+1+len(charged)%15])

        charged=[]
        for i in range(int(degree)):
            for j in range(backandside):
                if(j//backbone==0):
                    charged.append(int((i)*backandside+(j+1)))
        charged.append(total_pol-1)
        charged.append(total_pol)
        with open(path_output_oeff+moleculename+'_RU_SC_2.xyz', 'r') as f:
            total_soup=int(f.readlines()[0])
        for chain in range(nochainsinbox['2']):
            for i in range(int('2')):
                for j in range(backandside):
                    if(j//backbone==0):
                        charged.append(chain*total_soup+total_pol+int((i)*backandside+(j+1)))
            charged.append(total_pol+(chain+1)*total_soup-1)
            charged.append(total_pol+(chain+1)*total_soup)
        with open(path_output_dos+moleculename+'charge_'+degree+'_2.ndx', 'w') as f:
            f.write('[ backbone_'+degree+']\n')
            for i in range(len(charged)//15):
                f.writelines("%s " % l for l in charged[i*15:(i+1)*15])
                f.write('\n')
            f.writelines("%s " % l for l in charged[(i+1)*15:(i+1)*15+1+len(charged)%15])

    

