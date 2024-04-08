#!/usr/bin/env python
# coding: utf-8

'''
This is the first code to run in this step.

This code generates the coordinate (.xyz) and force field (.itp and .top according to GROMACS file formats) files
for the repeat unit structure specified in parameters.py.

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

# creat pair, angle, and dihedral lists
def angle_dihedral(bondlist):
    addbond=bondlist.copy()
    bonds=bondlist.copy()
    pairs=[]
    angles=[]
    editbond=[]
    dihedrals=[]
    addangle=[]
    editbondlines=['none']*3
    addanglelines=['none']*4
    adddihedral=[]
    adddihedrallines=['none']*5

    for each in addbond:
        addbondline = each.split()
        for item in bonds:
            bondlines = item.split()
            if addbondline[0]==bondlines[0]:
                if  bondlines[1]!=addbondline[1]:
                    addanglelines[0]=str(bondlines[1])
                    addanglelines[1]=str(bondlines[0])
                    addanglelines[2]=str(addbondline[1])
                    
                    addangleline='{0: <8}'.format(addanglelines[0])+'{0: <8}'.format(addanglelines[1])+'{0: <8}'.format(addanglelines[2])
                    addangle.append(addangleline)
                for seconditem in bonds:
                    secondbondlines = seconditem.split()
                    if bondlines[1]==secondbondlines[0] and secondbondlines[1]!=bondlines[0] and secondbondlines[0]!=addbondline[1]:
                                adddihedrallines[0]=str(secondbondlines[1])
                                adddihedrallines[1]=str(secondbondlines[0])
                                adddihedrallines[2]=str(bondlines[0])
                                adddihedrallines[3]=str(addbondline[1])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)
                    if bondlines[1]==secondbondlines[1] and secondbondlines[0]!=bondlines[0] and secondbondlines[1]!=addbondline[1]:
                                adddihedrallines[0]=str(secondbondlines[0])
                                adddihedrallines[1]=str(secondbondlines[1])
                                adddihedrallines[2]=str(bondlines[0])
                                adddihedrallines[3]=str(addbondline[1])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)
            if addbondline[0]==bondlines[1]:
                if  bondlines[0]!=addbondline[1]:
                    addanglelines[0]=str(bondlines[0])
                    addanglelines[1]=str(bondlines[1])
                    addanglelines[2]=str(addbondline[1])
                    
                    addangleline='{0: <8}'.format(addanglelines[0])+'{0: <8}'.format(addanglelines[1])+'{0: <8}'.format(addanglelines[2])
                    addangle.append(addangleline)
                for seconditem in bonds:
                    secondbondlines = seconditem.split()
                    if bondlines[0]==secondbondlines[0] and secondbondlines[1]!=bondlines[1] and secondbondlines[0]!=addbondline[1]:
                                adddihedrallines[0]=str(secondbondlines[1])
                                adddihedrallines[1]=str(secondbondlines[0])
                                adddihedrallines[2]=str(bondlines[1])
                                adddihedrallines[3]=str(addbondline[1])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)
                    if bondlines[0]==secondbondlines[1] and secondbondlines[0]!=bondlines[1] and secondbondlines[1]!=addbondline[1]:
                                adddihedrallines[0]=str(secondbondlines[0])
                                adddihedrallines[1]=str(secondbondlines[1])
                                adddihedrallines[2]=str(bondlines[1])
                                adddihedrallines[3]=str(addbondline[1])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)
            if addbondline[1]==bondlines[0]:
                if  bondlines[1]!=addbondline[0]:
                    addanglelines[0]=str(bondlines[1])
                    addanglelines[1]=str(bondlines[0])
                    addanglelines[2]=str(addbondline[0])
                    
                    addangleline=addangleline='{0: <8}'.format(addanglelines[0])+'{0: <8}'.format(addanglelines[1])+'{0: <8}'.format(addanglelines[2])
                    addangle.append(addangleline)
                for seconditem in bonds:
                    secondbondlines = seconditem.split()
                    if bondlines[1]==secondbondlines[0] and secondbondlines[1]!=bondlines[0] and secondbondlines[0]!=addbondline[0]:
                                adddihedrallines[0]=str(secondbondlines[1])
                                adddihedrallines[1]=str(secondbondlines[0])
                                adddihedrallines[2]=str(bondlines[0])
                                adddihedrallines[3]=str(addbondline[0])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)
                    if bondlines[1]==secondbondlines[1] and secondbondlines[0]!=bondlines[0] and secondbondlines[1]!=addbondline[0]:
                                adddihedrallines[0]=str(secondbondlines[0])
                                adddihedrallines[1]=str(secondbondlines[1])
                                adddihedrallines[2]=str(bondlines[0])
                                adddihedrallines[3]=str(addbondline[0])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)
            if addbondline[1]==bondlines[1]:
                if bondlines[0]!=addbondline[0]:
                    addanglelines[0]=str(bondlines[0])
                    addanglelines[1]=str(bondlines[1])
                    addanglelines[2]=str(addbondline[0])
                    
                    addangleline=addangleline='{0: <8}'.format(addanglelines[0])+'{0: <8}'.format(addanglelines[1])+'{0: <8}'.format(addanglelines[2])
                    addangle.append(addangleline)
                for seconditem in bonds:
                    secondbondlines = seconditem.split()
                    if bondlines[0]==secondbondlines[0] and secondbondlines[1]!=bondlines[1] and secondbondlines[0]!=addbondline[0]:
                                adddihedrallines[0]=str(secondbondlines[1])
                                adddihedrallines[1]=str(secondbondlines[0])
                                adddihedrallines[2]=str(bondlines[1])
                                adddihedrallines[3]=str(addbondline[0])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)
                    if bondlines[0]==secondbondlines[1] and secondbondlines[0]!=bondlines[1] and secondbondlines[1]!=addbondline[0]:
                                adddihedrallines[0]=str(secondbondlines[0])
                                adddihedrallines[1]=str(secondbondlines[1])
                                adddihedrallines[2]=str(bondlines[1])
                                adddihedrallines[3]=str(addbondline[0])
                                
                                adddihedralline='{0: <8}'.format(adddihedrallines[0])+'{0: <8}'.format(adddihedrallines[1])+'{0: <8}'.format(adddihedrallines[2])+'{0: <8}'.format(adddihedrallines[3])
                                adddihedral.append(adddihedralline)


    angle1lines=['none']*5
    toremove=[]
    import itertools
    #removing repitition of angles
    addangle=list(dict.fromkeys(addangle))
    for each in addangle:
        angle1 = each.split()
        index = addangle.index(each)
        for line in addangle[index:]:
            angle2 = line.split()
            if angle1[0]==angle2[2] and angle1[1]==angle2[1] and angle1[2]==angle2[0]:
                    toremove.append(index)            

    def delete_multiple_element(list_object, indices):
        indices = sorted(indices, reverse=True)
        for idx in indices:
            if idx < len(list_object):
                list_object.pop(idx)
    delete_multiple_element(addangle, toremove)


    dihed1lines=['none']*5
    toremove=[]


    #removing repitition of dihedrals
    adddihedral=list(dict.fromkeys(adddihedral))
    for each in adddihedral:
        dihed1 = each.split()
        index = adddihedral.index(each)
        for line in adddihedral[index:]:
            dihed2 = line.split()
            if dihed1[0]==dihed2[3] and dihed1[1]==dihed2[2] and dihed1[2]==dihed2[1] and dihed1[3]==dihed2[0]:
                    toremove.append(index)            

    def delete_multiple_element(list_object, indices):
        indices = sorted(indices, reverse=True)
        for idx in indices:
            if idx < len(list_object):
                list_object.pop(idx)
    delete_multiple_element(adddihedral, toremove)

    addpair=[]
    addpairlines=['none']*3
    #add new pairs
    for each in adddihedral:
        dihed1 = each.split()
        addpairlines[0]=dihed1[0]
        addpairlines[1]=dihed1[3]
        
        addpairline='{0: <8}'.format(addpairlines[0])+'{0: <8}'.format(addpairlines[1])
        addpair.append(addpairline)
    return addpair, addangle, adddihedral
#end find pair, angles, dihedrals
#calculate angle between three points
def angle_calculator(point_1, point_2, point_3):
    import numpy as np
    a = np.array(point_1)
    b = np.array(point_2)
    c = np.array(point_3)
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)
#end of angle calculator
#calculate dihedral angle of four points
def dihedral_calculator(point_1, point_2, point_3, point_4):
    import numpy as np
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
# extract xyz file from torsion log file
def find_input_orientation(fold_log, f_name, dihedral_l, folder_to_wr, write_to_file): 

    log_file=fold_log+f_name+".log"
    
    coordinates = []
    indexes = []
    indexes_2 = []
    stationary_points = []
    atom_count = 0
    atom_list = []
    
    
        
    with open(log_file, "r") as f:
        lines = f.readlines()
        
    #read the number of atoms in the molecule
    for i in range(0, len(lines)): 
        
        if "NAtoms=" in lines[i]:
            split1 = [elna for elna in lines[i].split(" ") if elna!=""]
            
            atom_count = int(split1[1])
            
            break
    
    #read the number of scan geometries done
    for i2 in range(0,len(lines)):
        
        if "Symbolic Z-matrix:" in lines[i2]:
            for j in range(i2+2,i2+atom_count+2):
                split2 = [elsz for elsz in lines[j].rstrip("\n").split(" ") if elsz!=""]
                
                atom_list.append(split2[0])
               
    
    #find indexes for the input orientations and stationary point geometries to be found
    for i3 in range(0, len(lines)):
        
        if "Standard orientation:" in lines[i3]:
            indexes.append(int(i3))
            
            
        if "Optimized Parameters" in lines[i3]:
            stationary_points.append(i3)
            
    
    
    #find the indexes of the input orientation tags for the stationary point geometries
    for i4 in range(0,len(stationary_points)):
        for j2 in range(0,len(indexes)-1):
            
            if indexes[j2] < stationary_points[i4] < indexes[j2+1]:
                indexes_2.append(indexes[j2])
    
    
    #read the lines containing the coordinates in the input orientation sections with indexes in indexes_2
    for i5 in indexes_2: 
        
        coordinates_sub = []
        
        for j3 in range(i5+5,i5+5+atom_count):
            
            split3 = lines[j3].rstrip("\n").split(" ")
            spl=[el for el in split3 if el !=""]
            
            coordinates_sub.append([atom_list[int(spl[0])-1],spl[3],spl[4],spl[5]])
            
        coordinates.append(coordinates_sub)
  
    #write to file (optional) in .xyz format
    if write_to_file: 
        
        for idhd in range(len(dihedral_l)):
            
            file_name = folder_to_wr + f_name.rstrip("tors") +'-'+ str(dihedral_l[idhd]) + ".xyz"

            crd_to_file = [crdline[0] + " "*4 + crdline[1] + " " + crdline[2]  + " " + crdline[3] for crdline in coordinates[idhd]]
            
            with open(file_name, "w") as k: 
                k.write(str(atom_count) + "\n" + f_name.rstrip("tors") +'-'+ str(dihedral_l[idhd]) + "\n")
                for lcrds in crd_to_file:
                    
                    k.write(lcrds+"\n")
                
                

    return (True)

'''
extract torsional potential from Gaussian log file of torsional scan. 
For each log file exists in Torsion folder, one txt file (actual value of potential at the exact dihedral angle) 
and one xvg file (a table from -180 to 180, in the format of GROMACS table torsion potential) will be generated.
'''

from parameters import *

parent_dir=cp.deepcopy(PATH_OUTPUT) #the Gaussian log files should be generated in the same folder as oligomer.xyz and .com files
moleculename = OLIGOMERNAME
path_output = os.path.join(parent_dir, moleculename)
path_output_torsion = os.path.join(parent_dir, moleculename+'/Torsion/') # where the charge and torsion files will be find
path_output_oeff=os.path.join(parent_dir, moleculename+'/') # where the force field files will be coppied
fragmentlist=cp.deepcopy(FRAGMENTLIST)
fragmentlist_temp=cp.deepcopy(FRAGMENTLIST)
fragmentlist_temp.insert(0,fragmentlist_temp[-1]) # add capping fragments (to the first end) to make oligomer
fragmentlist_temp.insert(len(fragmentlist_temp),fragmentlist_temp[1]) # add capping fragments (to the second end) to make oligomer
with open(PATH_FRAGMENTS+fragmentlist_temp[0]+'.xyz', 'r') as f:
    begin=f.readline()
with open(PATH_FRAGMENTS+fragmentlist_temp[-1]+'.xyz', 'r') as f:
    end=f.readline()
path_fragments=cp.deepcopy(PATH_FRAGMENTS)
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
        fragmentlist_natom.append(int(f.readlines()[0])-2)
fragmentlist_cumsum_second=np.cumsum(fragmentlist_natom)
fragmentlist_cumsum_first=np.insert(fragmentlist_cumsum_second[:-1], 0, 0)

all_torsion_atoms=[]
for i in range(0,len(torsion_pairs),2):
    all_torsion_atoms.append([fragmentlist_cumsum_first[torsion_pairs[i]]+2,fragmentlist_cumsum_first[torsion_pairs[i]]+1,fragmentlist_cumsum_second[torsion_pairs[i+1]],fragmentlist_cumsum_second[torsion_pairs[i+1]]-1])

print(all_torsion_atoms)
DIHEDRAL_ANGLES=[]
for i_index,i in enumerate(file_num_list):
    ang_rng=[0,36]
    #increment
    incr=1
    dhd_rngs=range(ang_rng[0],ang_rng[1]+1,incr)
    find_input_orientation(path_output_torsion, moleculename+"-torsion-"+str(i), dhd_rngs, path_output_torsion, True)
    
    dihedral_angles=[]
    for j in dhd_rngs:
        with open(path_output_torsion+moleculename+"-torsion-"+str(i)+'-'+str(j)+'.xyz', 'r') as f:
            lines=f.readlines()
            coord1=[float(lines[all_torsion_atoms[i_index][0]+1].split()[1]),float(lines[all_torsion_atoms[i_index][0]+1].split()[2]),float(lines[all_torsion_atoms[i_index][0]+1].split()[3])]
            coord2=[float(lines[all_torsion_atoms[i_index][1]+1].split()[1]),float(lines[all_torsion_atoms[i_index][1]+1].split()[2]),float(lines[all_torsion_atoms[i_index][1]+1].split()[3])]
            coord3=[float(lines[all_torsion_atoms[i_index][2]+1].split()[1]),float(lines[all_torsion_atoms[i_index][2]+1].split()[2]),float(lines[all_torsion_atoms[i_index][2]+1].split()[3])]
            coord4=[float(lines[all_torsion_atoms[i_index][3]+1].split()[1]),float(lines[all_torsion_atoms[i_index][3]+1].split()[2]),float(lines[all_torsion_atoms[i_index][3]+1].split()[3])]
        dihedral_angles.append(dihedral_calculator(coord1,coord2,coord3,coord4))
    DIHEDRAL_ANGLES.append(dihedral_angles) 
    str_tors=subprocess.check_output(["sed", "-n", "/N-N= /,/@/p", path_output_torsion+moleculename+"-torsion-"+str(i)+".log"])
    newl_tors=str_tors.decode('utf-8').replace("\n","").replace(" ","")
    tors_no_suffix_prefix=newl_tors.split("HF=")[1].split("RMSD=")[0].rstrip("\\").rstrip("/")
    energies_tors=tors_no_suffix_prefix.split(",")
    energies_tors_eV=np.array(energies_tors, float)*27.2114 # unit in eV 
    energies_tors_kJmol=np.array(energies_tors, float)*2625.5 # unit in kJ/mol
    energies_tors_eV_norm=energies_tors_eV-min(energies_tors_eV)
    energies_tors_kJmol_norm=energies_tors_kJmol-min(energies_tors_kJmol)
    energies_tors_dihedangle_kJmol_norm=np.append(np.expand_dims(dihedral_angles,axis=1), np.expand_dims(energies_tors_kJmol-min(energies_tors_kJmol),axis=1), 1)
    energies_tors_lines=[]

    for each in energies_tors_dihedangle_kJmol_norm:
        energies_tors_lines.append(str(round(each[0],1))+'    '+str(round(each[1],2)))
    with open(path_output_torsion+moleculename+"-torsion-"+str(i)+".txt", "w") as f:
        f.write("#"+str(all_torsion_atoms[i_index])+"\n")
        f.writelines("%s\n" % l for l in energies_tors_lines)
    # sort angle vs potential so that cubic spline is possible
    dihedral_angles_sorted = np.array(dihedral_angles[:-1]).argsort()
    sorted_dihedral_angles = np.array(dihedral_angles[:-1])[dihedral_angles_sorted]
    sorted_energies_tors_kJmol_norm = np.array(energies_tors_kJmol_norm[:-1])[dihedral_angles_sorted]
    
    # cubic spline fit on sorted data
    cs = CubicSpline(sorted_dihedral_angles, sorted_energies_tors_kJmol_norm)
    xs=np.linspace(-180., 180., 361)
    
    table=[]
    for j in range(len(xs)):
        table.append(str(xs[j])+'    '+str(round(cs(xs)[j],3))+'     '+str(-round(cs(xs, 1)[j],3)))

    # writing table potential
    with open(path_output_torsion+'table_d'+str(i)+'.xvg', 'w') as f:
        f.writelines("%s\n" % l for l in table)
    
    # plotting both actual data (10 degree spacing) and cubic spline fit (1 degree spacing)
    col=cp.deepcopy(COLOR)
    plt.plot(sorted_dihedral_angles,sorted_energies_tors_kJmol_norm, c=col[i], linestyle='none', marker='D')
    plt.plot(xs,cs(xs), c=col[i], linestyle='dashed')


plt.xlabel("dihedral angle [degree]", fontsize=14) 
plt.ylabel("energy [kJ/mol]", fontsize=14)
lgd=plt.legend(LEGEND,
	           bbox_to_anchor=(1, 0.5), prop={'size': 12})

plt.savefig(path_output_torsion+'dft.eps', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')

'''
Writing force field file for the oligomer.
'''
# extracting charges and rounding them so that the sum is absolute zero
word1=' ESP charges:'
word2=' Sum'
copy = False
with open(path_output_oeff+moleculename+'-charge.log', 'r') as fp:
    lines = fp.readlines()
    for line in lines:
        if (line.startswith(word1)):
            copy = True
            charges=[]
            continue
        elif (line.startswith(word2)):
            copy = False
            continue
        elif copy:
            charges.append(line.split())
subprocess.run(['obabel',path_output_oeff+moleculename+'-charge.log','-O', moleculename+'.xyz'], cwd = path_output_torsion)
subprocess.run(['mv',moleculename+'.xyz',path_output_oeff],cwd = path_output_torsion)
charges.pop(0)
CHARGES=[column[2] for column in charges]
CHARGES_float=[float(item) for item in CHARGES]

CHARGES_float_trunc=CHARGES_float[(int(begin)-1):-(int(end)-1)]
leftover=0-sum(CHARGES_float_trunc)
add_value=leftover/len(CHARGES_float_trunc)
CHARGES_float_corr=[item+add_value for item in CHARGES_float_trunc]

''' this function sometimes does not work
def round_retain_sum(x):
    x = x*1000 # We want 3 decimal precision
    N = np.round(np.sum(x)).astype(int)
    y = x.astype(int)
    M = np.sum(y)
    K = N - M 
    z = y-x 
    if K!=0:
        idx = np.argpartition(z,K)[:K]
        y[idx] += 1     
    return y/1000.
CHARGES_round=round_retain_sum(np.array(CHARGES_float))
'''
CHARGES_round = []
rest = 0.0
for n in CHARGES_float_corr:
    new_n = round(n + rest,3)
    rest += n - new_n
    CHARGES_round.append( new_n )

#end of read esp charge

atom_type=[]
atom=[]
bond=[]
'''
individual_fragments_code=['A','B','C','D','E','F','G']
individual_fragments_name=['IDT','T','BTT','DPP','BT','BTZ','BPD']
residue_name=[]
with open(path+'.xyz', 'r') as f:#the fragmentlist should be in the second line of xyz file
    fragmentlist=[ch for ch in f.readlines()[1] if ch != '\n' if ch != ' ']
for i in fragmentlist:
    residue_name.append(individual_fragments_name[individual_fragments_code.index(i)])
'''
counter_atom_type=0
counter_atom=0
resnr=1
cgnr=1
bond_type=1 #harmonic
kb=320000.0 # bond constant
pair_type=1
angle_type=1 #harmonic
ka=500.0 # angle constant
dihedral_type=1 #proper dihedral with Phi (dihedral angle) and a force constant (kJ/mol)
kd=10.0 #
multiplicity=1

flat_dihedral_value='{0: <6}'.format(3)+'{0: <8}'.format(30.334)+'{0: <8}'.format(0.000)+'{0: <8}'.format(-30.334)+'{0: <8}'.format(0.000)+'{0: <8}'.format(0.000)+'{0: <8}'.format(0.000) 
fitted_dihedral_funct={}
fitted_dihedral_atoms=[]
for i in file_num_list:
	fitted_dihedral_funct["torsion{0}".format(i)] = '{0: <6}'.format(8)+'{0: <8}'.format(i+1)+'{0: <8}'.format(1) # i+1 because we remove the first fragment in the monomer





improper_type=4 #improper dihedral with Phi (dihedral angle) and a force constant (kJ/mol)
kim=10.0
multiplicity_improper=2
improper_value = 180.0 # all impropers are flat
for counter_frag, fragname in enumerate(fragmentlist):
    coordinate=[]
    with open(path_fragments+fragname+'.xyz', 'r') as f:
            lines=f.readlines()[2:-2]
            number_of_atoms=len(lines)+2
            for line in lines:
                coordinate.append(line.split()[1:4])
            COORDINATE=np.array(coordinate).astype(float)
    with open(path_fragments+fragname+'.mol', 'r') as f:
        lines=f.readlines()[4+number_of_atoms:-1]
        for counter_bond, line in enumerate(lines):
            bond_lines=line.split()
            if (int(bond_lines[0])==number_of_atoms or int(bond_lines[1])==number_of_atoms or int(bond_lines[0])==number_of_atoms-1 or int(bond_lines[1])==number_of_atoms-1):
                continue
            bond_line='{0: <8}'.format(str(int(bond_lines[0])+counter_atom_type))+'{0: <8}'.format(str(int(bond_lines[1])+counter_atom_type))+'{0: <6}'.format(str(bond_type))+'{0: <10}'.format("{:.4f}".format(np.linalg.norm(COORDINATE[int(line.split()[0])-1]-COORDINATE[int(line.split()[1])-1])*0.1))+'{0: <10}'.format(str(kb))
            bond.append(bond_line)
    with open(path_fragments+fragname+'.atp', 'r') as f:
        lines=f.readlines()[2:-2]
        for line in lines:
            atom_typelines=line.split()
            atom_typeline ='{0: <12}'.format('oeff_'+str(counter_atom_type+1))+'{0: <6}'.format(atom_typelines[0])+'{0: <6}'.format(atom_typelines[1])+'{0: <12}'.format(atom_typelines[2])+'{0: <12}'.format(atom_typelines[3])+'{0: <6}'.format(atom_typelines[4])+'{0: >15}'.format(atom_typelines[5])+'{0: >15}'.format(atom_typelines[6])
            atom_line ='{0: <6}'.format(str(counter_atom_type+1))+'{0: <12}'.format('oeff_'+str(counter_atom_type+1))+'{0: <6}'.format(str(resnr))+'{0: <8}'.format(fragmentlist[counter_frag])+'{0: <6}'.format(atom_typelines[0])+'{0: <6}'.format(str(cgnr))+'{0: >12}'.format(CHARGES_round[counter_atom_type])+'{0: >12}'.format(atom_typelines[2])
            atom_type.append(atom_typeline)
            atom.append(atom_line)
            counter_atom_type+=1
# add inter-fragment bonds
inter_bond_atoms_1=[]
inter_bond_atoms_2=[]
counter_atom=0
if(len(fragmentlist)>1):
    for counter_frag, fragname in enumerate(fragmentlist):
        if (counter_frag==0):
            with open(path_fragments+fragname+'.xyz', 'r') as f:
                    lines=f.readlines()[2:-2]
                    inter_bond_atoms_1.append(counter_atom+1)
                    for each in lines:
                        counter_atom+=1
        else:
            with open(path_fragments+fragname+'.xyz', 'r') as f:
                    lines=f.readlines()[2:-2]
                    if (counter_frag!=len(fragmentlist)-1):
                        inter_bond_atoms_1.append(counter_atom+1)
                    for each in lines:
                        counter_atom+=1
                    inter_bond_atoms_2.append(counter_atom)
#making corrdinate file for repeat unit
with open(path_output_oeff+moleculename+'.xyz', 'r') as f:
    lines_RU=f.readlines()[2+(int(begin)-1):-(int(end)-1)]
with open(path_output_oeff+moleculename+'_RU.xyz', 'w') as f:
    f.write(str(len(lines_RU))+"\n\n")
    f.writelines("%s" % l for l in lines_RU)

oligomer_coordinate=[]
with open(path_output_oeff+moleculename+'_RU.xyz', 'r') as f:
    lines=f.readlines()[2:]
    for each in lines:
        coordinate=each.split()
        oligomer_coordinate.append([float(coordinate[1]), float(coordinate[2]), float(coordinate[3])])   

bond_inter=[]
bond_addedinter=cp.deepcopy(bond)
OLIGOMER_COORDINATE=np.array(oligomer_coordinate)
for i in range(len(inter_bond_atoms_1)):
    bond_inter.append([inter_bond_atoms_1[i], inter_bond_atoms_2[i]])
for i in range(len(bond_inter)):
    bond_line='{0: <8}'.format(str(int(bond_inter[i][0])))+'{0: <8}'.format(str(int(bond_inter[i][1])))+'{0: <6}'.format(str(bond_type))+'{0: <10}'.format("{:.4f}".format(np.linalg.norm(OLIGOMER_COORDINATE[int(bond_inter[i][0])-1]-OLIGOMER_COORDINATE[int(bond_inter[i][1])-1])*0.1))+'{0: <10}'.format(str(kb))
    bond_addedinter.append(bond_line)
# end add inter-fragment bonds

bondlist_temp=[]
for item in bond_addedinter:
	bondlist_temp.append([int(item.split()[0]),int(item.split()[1])])

# find sp3 carbons
temp=[]
for each in bondlist_temp:
    temp.append(each[0])
    temp.append(each[1])  
sp3_carbons=[item for item, count in collections.Counter(temp).items() if count > 3]

pairs=[s + '{0: <6}'.format(str(pair_type)) for s in angle_dihedral(bond)[0]]
angles=[s + '{0: <6}'.format(str(angle_type)) for s in angle_dihedral(bond)[1]]
# calculate the angle between three points in angle
for i, each in enumerate(angles):
    angle_value=angle_calculator(oligomer_coordinate[int(each.split()[0])-1],oligomer_coordinate[int(each.split()[1])-1],oligomer_coordinate[int(each.split()[2])-1])
    angles[i]=each+'{0: <8}'.format(str(round(angle_value,1)))+'{0: <8}'.format(str(ka))
    
dihedrals=[s + '{0: <6}'.format(str(dihedral_type)) for s in angle_dihedral(bond)[2]]


          
for i, each in enumerate(dihedrals):
    dihedral_value=dihedral_calculator(oligomer_coordinate[int(each.split()[0])-1], oligomer_coordinate[int(each.split()[1])-1], oligomer_coordinate[int(each.split()[2])-1], oligomer_coordinate[int(each.split()[3])-1])
    #if (abs(dihedral_value)>179.0 or abs(dihedral_value)<1.0):
    if (int(each.split()[0]) in sp3_carbons or int(each.split()[1]) in sp3_carbons or int(each.split()[2]) in sp3_carbons or int(each.split()[3]) in sp3_carbons):
        dihedrals[i]=each+'{0: <8}'.format(str(round(dihedral_value-180,1)))+'{0: <8}'.format(str(kd))+'{0: <8}'.format(str(multiplicity))
    else:
        dihedrals[i]=each+flat_dihedral_value


# add inter-fragment pair angle dihedral (all dihedrals added)
pairs_addinter=[s + '{0: <6}'.format(str(pair_type)) for s in angle_dihedral(bond_addedinter)[0]]
angles_addinter=[s + '{0: <6}'.format(str(angle_type)) for s in angle_dihedral(bond_addedinter)[1]]
for i, each in enumerate(angles_addinter):
    angle_value=angle_calculator(oligomer_coordinate[int(each.split()[0])-1],oligomer_coordinate[int(each.split()[1])-1],oligomer_coordinate[int(each.split()[2])-1])
    angles_addinter[i]=each+'{0: <8}'.format(str(round(angle_value,1)))+'{0: <8}'.format(str(ka))


#dihedrals_addinter=[s + '{0: <6}'.format(str(dihedral_type)) for s in angle_dihedral(bond_addedinter)[2]]
dihedrals_addinter=angle_dihedral(bond_addedinter)[2]
flat_bond_inter = [item for sublist in bond_inter for item in sublist]

tobe_del=[]
for i, each in enumerate(dihedrals_addinter):
    line=each.split()
    #if ((int(line[1]) in flat_bond_inter and int(line[2]) in flat_bond_inter) or (int(line[0]) in flat_bond_inter and int(line[1]) in flat_bond_inter) or (int(line[2]) in flat_bond_inter and int(line[3]) in flat_bond_inter)):
    if (int(line[1]) in flat_bond_inter and int(line[2]) in flat_bond_inter):
        check=False
        for j, item in enumerate(all_torsion_atoms):
            if((int(line[0])==all_torsion_atoms[j][0] and int(line[3])==all_torsion_atoms[j][3]) or (int(line[0])==all_torsion_atoms[j][3] and int(line[3])==all_torsion_atoms[j][0])):
                for every_index, every in enumerate(all_torsion_atoms):
                	#***** for asymmetric fragments the condition should be changed ****
                	if(int(line[1])==every[1] and int(line[2])==every[2]):
                		dihedrals_addinter[i]=each+fitted_dihedral_funct['torsion'+str(every_index+1)]
                check=True
        if(check==False):
            tobe_del.append(i)
    else:
        dihedral_value=dihedral_calculator(oligomer_coordinate[int(each.split()[0])-1], oligomer_coordinate[int(each.split()[1])-1], oligomer_coordinate[int(each.split()[2])-1], oligomer_coordinate[int(each.split()[3])-1])
        #if (abs(dihedral_value)>179.0 or abs(dihedral_value)<1.0):
        if (int(each.split()[0]) in sp3_carbons or int(each.split()[1]) in sp3_carbons or int(each.split()[2]) in sp3_carbons or int(each.split()[3]) in sp3_carbons):
            dihedrals_addinter[i]=each+'{0: <6}'.format(str(dihedral_type))+'{0: <8}'.format(str(round(dihedral_value-180,1)))+'{0: <8}'.format(str(kd))+'{0: <8}'.format(str(multiplicity))
        else:
            dihedrals_addinter[i]=each+flat_dihedral_value
dihedrals_addinter_mod=[]
for i, each in enumerate(dihedrals_addinter):
    if (i in tobe_del):
        continue
    else:
        dihedrals_addinter_mod.append(dihedrals_addinter[i])
    
#dihedrals_inter=list(set(dihedrals_addinter_mod) - set(dihedrals))
# end add inter-fragment pair angle dihedral (all dihedrals added)



# writing itp and top files
with open(path_output_oeff+moleculename+'_RU.itp', "w") as f:
    f.write(";\n; GENERATED BY GAMMPS\n; Troisi Lab @ University of Liverpool\n")
    f.write("[ atomtypes ]\n")
    f.writelines("%s\n" % l for l in atom_type)
    f.write("[ moleculetype ]\n;Name     nrexcl\n")
    f.writelines(moleculename+'     3\n')
    f.write("[ atoms ]\n")
    f.writelines("%s\n" % l for l in atom)
    f.write("[ bonds ]\n")
    f.writelines("%s\n" % l for l in bond_addedinter)
#    f.write("[ pairs ]\n")
#    f.writelines("%s\n" % l for l in pairs_addinter)
    f.write("[ angles ]\n")
    f.writelines("%s\n" % l for l in angles_addinter)
    f.write("[ dihedrals ]\n")
    f.writelines("%s\n" % l for l in dihedrals_addinter_mod)


with open(path_output_oeff+moleculename+'_RU.top', "w") as f:
    f.write("[ defaults ]\n")
    f.write('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
    f.write('1     3     yes     0.5     0.5\n\n')
    f.write('#include "'+moleculename+'_RU.itp"\n\n')
    f.write('[ system ]\n; Name\n'+moleculename+'\n\n')
    f.write('[ molecules ]\n'+moleculename+'  1\n')	
end = time.time()


with open(path_output_oeff+'/RU.tcl', "w") as f:
    f.write('mol new '+moleculename+'_RU.xyz type {xyz}\n\nmol modstyle 0 0 CPK 1.000000 0.300000 120.000000 120.000000\n\nmol modcolor 0 0 Name\n\ndisplay resetview\n\naxes location Off\n\n')
    f.write("color Display Background white\n\ndisplay projection orthographic\n\n"+"render Tachyon "+moleculename+"_RU "+VMD_RENDER+" -aasamples 12 %s -format TARGA -o %s.tga -res 1000 1000"+"\nexit\n")
subprocess.run(['vmd','-dispdev','text','-e', 'RU.tcl'], cwd = path_output_oeff)
subprocess.run(['magick', moleculename+'_RU.tga', moleculename+'_RU.jpeg'], cwd = path_output_oeff)

with open(path_output+'.html', "a") as f:
    f.write('\n\n<p><font size = "+2">\nStep 2: Building up repeat unit and polymer coordinate and force field files </font> </p>\n\n')
    f.write('<aside class="figures">\n')
    f.write('<figure>\n<img src="'+moleculename+'/'+moleculename+'_RU.jpeg" style="width: 500px; height: 300px; object-fit: cover;" alt="RU">\n<figcaption>repeat unit</figcaption>\n</figure>')
    f.write('</aside>\n')

