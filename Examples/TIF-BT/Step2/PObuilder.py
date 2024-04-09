#!/usr/bin/env python
# coding: utf-8

import numpy as np
import copy as cp
import math
from scipy.spatial.transform import Rotation as R
import collections
import subprocess
from matplotlib import pyplot as plt
import os
import time
import re


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix
def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy
def inbetweenlines(a,b):
    copy = False
    inbetween=[]
    for line in lines:
        if line.strip() == a:
            copy = True
            continue
        elif line.strip() == b:
            copy = False
            continue
        elif copy:
            inbetween.append(line.split())
    return inbetween


start = time.time()
from parameters import *

#from OligomerBuilder import fiveringmembers
fiveringmembers=FIVERINGS#five rings should be treated differently

sixringconnect=SIXRINGCONNECT#fragments which have a six-atom ring at the connection with another fragment

parent_dir=cp.deepcopy(PATH_OUTPUT) #the Gaussian log files should be generated in the same folder as oligomer.xyz and .com files
moleculename = OLIGOMERNAME
path_output = os.path.join(parent_dir, moleculename)
path_output_oeff=os.path.join(parent_dir, moleculename+'/') # where the force field files will be coppied
fragmentlist=cp.deepcopy(FRAGMENTLIST)
degree_polymer=cp.deepcopy(DP)
tact=cp.deepcopy(TACTICITY)

for degree in degree_polymer:
    if (tact[0]=='syn'):
        sequence_alignement=['up','down']*int(int(degree)/2) #syndiotactic
    elif(tact[0]=='iso'):
        sequence_alignement=['up','up']*int(int(degree)/2) #isotactic
    #dimer_fragmentlist=[]
    H_C_bond_length=1.090
    C_C_bond_length=1.5290
    bond_type=1
    bond_stiffness=320000.0
    angle_type=1
    angle_stiffness=500.0
    dihedral_type=3
    fitted_dihedral_funct_corr=['8', '1', '1'] #table
    flat_dihedral_type=3
    flat_dihedral=['30.334', '0.000', '-30.334', '0.000', '0.000', '0.000']

    improper_type=4
    improper_value=0.0
    improper_stiffness=10.0
    multiplicity=2
    #for i in range(len(sequence_fragmentlist)-1):
    #    dimer_fragmentlist.append(sequence_fragmentlist[i][-1]+sequence_fragmentlist[i+1][0])
    #dimer_fragmentlist=list(dict.fromkeys(dimer_fragmentlist))  # all fragment connections
    inter_bond_info=[] #will contain n (number of monomers in polmer) list, each contains the connecting points of each monomer
    number_atom=[]#number of atoms of the conjugated part of each monomer in the sequence
    number_atom_total=[]#total number of atoms of each monomer including side chains
    for i in range(int(degree)):
        with open(PATH_FRAGMENTS+fragmentlist[0]+'.xyz', 'r') as f:
            temp=int(f.readline())-2
            inter_bond_info.append([temp])
        with open(PATH_FRAGMENTS+fragmentlist[-1]+'.xyz', 'r') as f:
            temp=int(f.readline())-2
        with open(path_output_oeff+moleculename+'_RU.xyz', 'r') as f:
            inter_bond_info[i].append(int(f.readline())-temp+1)
        with open(path_output_oeff+moleculename+'_RU.xyz', 'r') as f:
            number_atom.append(int(f.readline()))
        with open(path_output_oeff+moleculename+'_RU_SC.xyz', 'r') as f:
            number_atom_total.append(int(f.readline()))
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
    coordinates=[]# calculate anlge values for inter-agnle parameters for monomer-monomer connection

    with open(path_output_oeff+moleculename+'.xyz', 'r') as f:
        coords=f.readlines()[2:]
        for coord in coords:
            coordinates.append([float(coord.split()[1]), float(coord.split()[2]),float(coord.split()[3])])
    fragmentlist_olig=cp.deepcopy(fragmentlist)
    fragmentlist_olig.insert(0,fragmentlist_olig[-1]) # add capping fragments (to the first end) to make oligomer
    fragmentlist_olig.insert(len(fragmentlist_olig),fragmentlist_olig[1]) # add capping fragments (to the second end) to make oligomer
    with open(PATH_FRAGMENTS+fragmentlist_olig[-1]+'.xyz', 'r') as f:
        last_fragment_size=int(f.readline())
    with open(PATH_FRAGMENTS+fragmentlist_olig[-2]+'.xyz', 'r') as f:
        second_last_fragment_size=int(f.readline())

    a=[coordinates[-2], coordinates[-3], coordinates[-4]]
    b=[coordinates[-(last_fragment_size-1)-(second_last_fragment_size-1)+1],coordinates[-(last_fragment_size-1)-(second_last_fragment_size-1)+2],coordinates[-(last_fragment_size-1)-(second_last_fragment_size-1)+3]]
    angle_value1=angle_calculator(b[1],b[0],a[0])
    angle_value2=angle_calculator(b[2],b[0],a[0])
    angle_value3=angle_calculator(a[1],a[0],b[0])
    angle_value4=angle_calculator(a[2],a[0],b[0])
    bond_value=np.linalg.norm(np.array(a[0])-np.array(b[0])) 
    count=0
    coordinates=[]
    coordinates_noatomname=[]
    coordinates_atomname=[]
    with open(path_output_oeff+moleculename+'_RU_SC.xyz', 'r') as f: #for the first monomer
        coords=f.readlines()[2:]
        for coord in coords:
            coordinates.append(coord.split())
            coordinates_noatomname.append([float(coord.split()[1]), float(coord.split()[2]),float(coord.split()[3])])
            coordinates_atomname.append([coord.split()[0]])
        coordinates_atomname_numpy=np.array(coordinates_atomname,dtype=object)
        coordinates_noatomname_numpy=np.array(coordinates_noatomname)
        vector=coordinates_noatomname_numpy[inter_bond_info[count][1]-1]-coordinates_noatomname_numpy[inter_bond_info[count][0]-1]
        vector_unit=vector/np.linalg.norm(vector)
        rot_mat=rotation_matrix_from_vectors(vector_unit,np.array([1.0, 0.0, 0.0]))#align the molecule towards x axis
        coordinates_noatomname_numpy_rotx=np.empty([len(coordinates_noatomname_numpy),3])
        for i, each in enumerate(coordinates_noatomname_numpy):
            coordinates_noatomname_numpy_rotx[i]=np.matmul(rot_mat,each)
        #coordinates_noatomname_numpy_moved_rotx=np.matmul(rot_mat,coordinates_noatomname_numpy_moved)
        coordinates_noatomname_numpy_moved_rotx=coordinates_noatomname_numpy_rotx-coordinates_noatomname_numpy_rotx[inter_bond_info[count][1]-1]

        rot_temp=[]
        sum_coor_xy=[]
        vector=coordinates_noatomname_numpy_moved_rotx[inter_bond_info[count][0]-1]-coordinates_noatomname_numpy_moved_rotx[inter_bond_info[count][1]-1]
        vector_unit=vector/np.linalg.norm(vector)
        
        
        for i in range(360):#rotate to have the conjugated atoms in the xy plane
            rotation_degrees=float(i)
            rotation_radians = np.radians(rotation_degrees)
            rotation_vector = rotation_radians * vector_unit
            rotation = R.from_rotvec(rotation_vector)
            rotated_vec = rotation.apply(coordinates_noatomname_numpy_moved_rotx)
            rot_temp.append(rotated_vec)
            sum_coor_xy.append(abs(rotated_vec[inter_bond_info[count][1]-1:inter_bond_info[count][1]-1+3,2]).sum())
        coordinates_noatomname_numpy_moved_rotz=cp.deepcopy(rot_temp[sum_coor_xy.index(min(sum_coor_xy))])

        
        rot_temp=[]
        sum_coor_z=[]
        for i in range(360):#rotate the conjugated atoms in the xy plane to have connecting carbon in the front ([0.0, 0.0, 0.0])
            rotation_degrees=float(i)
            rotation_radians = np.radians(rotation_degrees)
            rotation_vector = rotation_radians * np.array([0.0, 0.0, 1.0])
            rotation = R.from_rotvec(rotation_vector)
            rotated_vec = rotation.apply(coordinates_noatomname_numpy_moved_rotz)
            rot_temp.append(rotated_vec)
            sum_coor_z.append(rotated_vec[inter_bond_info[count][1]-1:inter_bond_info[count][1]-1+1,0].sum())
        coordinates_noatomname_numpy_aligned=rot_temp[sum_coor_z.index(min(sum_coor_z))]
        if (ROTATEFIRSTMONOMER=='yes'):#fix the tacticity
            rotation_radians = np.radians(180)
            rotation_vector = rotation_radians * vector_unit
            rotation = R.from_rotvec(rotation_vector)
            coordinates_noatomname_numpy_aligned = rotation.apply(coordinates_noatomname_numpy_aligned)
        polymer_coordinate_noatomname=cp.deepcopy(coordinates_noatomname_numpy_aligned)
        polymer_coordinate_atomname=cp.deepcopy(coordinates_atomname_numpy)
        
    
    for count in range(int(degree)-1):
        coordinates=[]
        coordinates_noatomname=[]
        coordinates_atomname=[]
        
            
        
        with open(path_output_oeff+moleculename+'_RU_SC.xyz', 'r') as f:
            coords=f.readlines()[2:]
            for coord in coords:
                coordinates.append(coord.split())
                coordinates_noatomname.append([float(coord.split()[1]), float(coord.split()[2]),float(coord.split()[3])])
                coordinates_atomname.append([coord.split()[0]])
            coordinates_atomname_numpy=np.array(coordinates_atomname,dtype=object)
            coordinates_noatomname_numpy=np.array(coordinates_noatomname)
            
            vector=coordinates_noatomname_numpy[inter_bond_info[count][0]-1]-coordinates_noatomname_numpy[inter_bond_info[count][1]-1]
            vector_unit=vector/np.linalg.norm(vector)
            rot_mat=rotation_matrix_from_vectors(vector_unit,np.array([1.0, 0.0, 0.0]))#align the molecule towards x axis
            coordinates_noatomname_numpy_rotx=np.empty([len(coordinates_noatomname_numpy),3])
            for i, each in enumerate(coordinates_noatomname_numpy):
                coordinates_noatomname_numpy_rotx[i]=np.matmul(rot_mat,each)
            if (coordinates_noatomname_numpy_rotx[inter_bond_info[count][0]-1][0]>coordinates_noatomname_numpy_rotx[inter_bond_info[count][1]-1][0]):
                coordinates_noatomname_numpy_moved_rotx=coordinates_noatomname_numpy_rotx-coordinates_noatomname_numpy_rotx[inter_bond_info[count][0]-1]-[bond_value, 0, 0]
            elif (coordinates_noatomname_numpy_rotx[inter_bond_info[count][0]-1][0]<coordinates_noatomname_numpy_rotx[inter_bond_info[count][1]-1][0]):
                coordinates_noatomname_numpy_moved_rotx=coordinates_noatomname_numpy_rotx-coordinates_noatomname_numpy_rotx[inter_bond_info[count][0]-1]+[bond_value, 0, 0]
            #coordinates_noatomname_numpy_moved_rotx=np.matmul(rot_mat,coordinates_noatomname_numpy_moved)
            vector=coordinates_noatomname_numpy_moved_rotx[inter_bond_info[count][1]-1]-coordinates_noatomname_numpy_moved_rotx[inter_bond_info[count][0]-1]
            vector_unit=vector/np.linalg.norm(vector)
            rot_temp=[]
            sum_coor_xy=[]
            
            
            for i in range(360):#rotate to have the conjugated atoms in the xy plane
                rotation_degrees=float(i)
                rotation_radians = np.radians(rotation_degrees)
                rotation_vector = rotation_radians * vector_unit
                rotation = R.from_rotvec(rotation_vector)
                rotated_vec = rotation.apply(coordinates_noatomname_numpy_moved_rotx)
                rot_temp.append(rotated_vec)
                sum_coor_xy.append(abs(rotated_vec[inter_bond_info[count][1]-1:inter_bond_info[count][1]-1+3,2]).sum())
            coordinates_noatomname_numpy_moved_rotz=cp.deepcopy(rot_temp[sum_coor_xy.index(min(sum_coor_xy))])

            rot_temp=[]
            sum_coor_z=[]
            
            for i in range(360):#rotate the conjugated atoms in the xy plane to have connecting carbon in the back ([1.5, 0.0, 0.0])
                rotation_degrees=float(i)
                rotation_radians = np.radians(rotation_degrees)
                rotation_vector = rotation_radians * np.array([0.0, 0.0, 1.0])
                rotation = R.from_rotvec(rotation_vector)
                rotated_vec = rotation.apply(coordinates_noatomname_numpy_moved_rotz)
                rot_temp.append(rotated_vec)
                sum_coor_z.append(rotated_vec[inter_bond_info[count][1]-1:inter_bond_info[count][1]-1+1,0].sum())
            coordinates_noatomname_numpy_aligned=rot_temp[sum_coor_z.index(max(sum_coor_z))]
            
            if(sequence_alignement[count]=='down'): #fix the tacticity
                rotation_radians = np.radians(180)
                rotation_vector = rotation_radians * vector_unit
                rotation = R.from_rotvec(rotation_vector)
                coordinates_noatomname_numpy_aligned = rotation.apply(coordinates_noatomname_numpy_aligned)  
            
            polymer_coordinate_noatomname=np.append(polymer_coordinate_noatomname,coordinates_noatomname_numpy_aligned,axis=0)
            polymer_coordinate_atomname=np.append(polymer_coordinate_atomname,coordinates_atomname_numpy,axis=0)
            polymer_coordinate_noatomname=polymer_coordinate_noatomname-polymer_coordinate_noatomname[np.cumsum(number_atom_total[:count+1])[-1]+inter_bond_info[count][1]-1]
            
    #adding H to the ends
    #redist_flag=input('Would you like the capping groups set as united atom CH3? (yes/no) ')
    #if (redist_flag=='yes'):
    ''' if one wants to cap the polymer with C3 super atom
    polymer_coordinate_atomname=np.append(polymer_coordinate_atomname,np.array([['C']]),axis=0)
    polymer_coordinate_atomname=np.append(polymer_coordinate_atomname,np.array([['C']]),axis=0)

    polymer_coordinate_noatomname=np.append(polymer_coordinate_noatomname,np.array([polymer_coordinate_noatomname[inter_bond_info[0][0]-1]-np.array([C_C_bond_length, 0.0, 0.0])]),axis=0) #first H
    polymer_coordinate_noatomname=np.append(polymer_coordinate_noatomname,np.array([polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]-1]+np.array([C_C_bond_length, 0.0, 0.0])]),axis=0)
    '''
    #else: cap the polymer with H atoms
    polymer_coordinate_atomname=np.append(polymer_coordinate_atomname,np.array([['H']]),axis=0)
    polymer_coordinate_atomname=np.append(polymer_coordinate_atomname,np.array([['H']]),axis=0)

    polymer_coordinate_noatomname=np.append(polymer_coordinate_noatomname,np.array([polymer_coordinate_noatomname[inter_bond_info[0][0]-1]-np.array([H_C_bond_length, 0.0, 0.0])]),axis=0) #first H
    if (int(degree)>1):
        polymer_coordinate_noatomname=np.append(polymer_coordinate_noatomname,np.array([polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]-1]+np.array([H_C_bond_length, 0.0, 0.0])]),axis=0)
    else:
        polymer_coordinate_noatomname=np.append(polymer_coordinate_noatomname,np.array([polymer_coordinate_noatomname[inter_bond_info[-1][-1]-1]+np.array([H_C_bond_length, 0.0, 0.0])]),axis=0)
    coordinates_concat=np.append(polymer_coordinate_atomname,polymer_coordinate_noatomname,axis=1)
    #add bonds to FF
    addbond=[]
    number_atom_total_temp=np.insert(number_atom_total,0,0)
    for count in range(int(degree)-1):
        addbond.append([str(inter_bond_info[count][1]+np.cumsum(number_atom_total_temp)[count]), str(inter_bond_info[count+1][0]+np.cumsum(number_atom_total)[count]), str(bond_type), str(round(bond_value*0.1,4)), str(bond_stiffness)])
    # add anlges to FF
    addangle=[]
    for each in addbond:
        addangle.append([str(int(each[0])+1), str(each[0]), str(each[1]), str(angle_type), str(round(angle_value1,1)) , str(angle_stiffness)])
        addangle.append([str(int(each[0])+2), str(each[0]), str(each[1]), str(angle_type), str(round(angle_value2,1)) , str(angle_stiffness)])
        addangle.append([str(each[0]), str(each[1]), str(int(each[1])-1), str(angle_type), str(round(angle_value3,1)) , str(angle_stiffness)]) 
        #should be corrected in case of 5-ring fragments
        addangle.append([str(each[0]), str(each[1]), str(int(each[1])-2), str(angle_type), str(round(angle_value4,1)) , str(angle_stiffness)])
    # add dihedral to FF

    # writing FF file
    atom=[]
    bond=[]
    pair=[]
    angle=[]
    dihedral=[]
    improper=[]#for capping hydrogen
    atoms_monomer=0
    for i in range(int(degree)):
        with open(path_output_oeff+moleculename+'_RU_SC.itp', 'r') as f:
            lines=f.readlines()
            atomtype=inbetweenlines('[ atomtypes ]', '[ moleculetype ]')
            atom_monomer=inbetweenlines('[ atoms ]', '[ bonds ]')
            for item in atom_monomer:
                item[0]=str(int(item[0])+atoms_monomer)
            atom.extend(atom_monomer)
            bond_monomer=inbetweenlines('[ bonds ]', '[ pairs ]')
            for item in bond_monomer:
                item[0]=str(int(item[0])+atoms_monomer)
                item[1]=str(int(item[1])+atoms_monomer)
            bond.extend(bond_monomer)
            pair_monomer=inbetweenlines('[ pairs ]', '[ angles ]')
            for item in pair_monomer:
                item[0]=str(int(item[0])+atoms_monomer)
                item[1]=str(int(item[1])+atoms_monomer)
            pair.extend(pair_monomer)
            angle_monomer=inbetweenlines('[ angles ]', '[ dihedrals ]')
            for item in angle_monomer:
                item[0]=str(int(item[0])+atoms_monomer)
                item[1]=str(int(item[1])+atoms_monomer)
                item[2]=str(int(item[2])+atoms_monomer)
            angle.extend(angle_monomer)
            dihedral_monomer=inbetweenlines('[ dihedrals ]', '')
            for item in dihedral_monomer:
                item[0]=str(int(item[0])+atoms_monomer)
                item[1]=str(int(item[1])+atoms_monomer)
                item[2]=str(int(item[2])+atoms_monomer)
                item[3]=str(int(item[3])+atoms_monomer)
            dihedral.extend(dihedral_monomer)
        with open(path_output_oeff+moleculename+'_RU_SC.xyz', 'r') as f:
            atoms_monomer+=int(f.readline())

    adddihedral=[]
    for each in addbond:# the ones which addbond is in the middle
        #should be corrected in case of 5-ring fragments
        adddihedral.append([str(int(each[0])+1),str(each[0]), str(each[1]), str(int(each[1])-1), fitted_dihedral_funct_corr[0], fitted_dihedral_funct_corr[1], fitted_dihedral_funct_corr[2]])
    temp1=[item[0] for item in addbond]#all firs atoms in inter-fragment bonds in one list
    temp2=[item[1] for item in addbond]#all second atoms in inter-fragment bonds in one list
    test=[]
    for each in angle:# the ones which addbond is in the middle
        if(each[0] in temp1 and each[1] not in temp2):
            adddihedral.append([temp2[temp1.index(each[0])], each[0], each[1], each[2], str(flat_dihedral_type), flat_dihedral[0], flat_dihedral[1], flat_dihedral[2], flat_dihedral[3], flat_dihedral[4], flat_dihedral[5]])
        if(each[2] in temp1 and each[1] not in temp2):
            adddihedral.append([temp2[temp1.index(each[2])], each[2], each[1], each[0], str(flat_dihedral_type), flat_dihedral[0], flat_dihedral[1], flat_dihedral[2], flat_dihedral[3], flat_dihedral[4], flat_dihedral[5]])
        if(each[0] in temp2 and each[1] not in temp1):
            adddihedral.append([temp1[temp2.index(each[0])], each[0], each[1], each[2], str(flat_dihedral_type), flat_dihedral[0], flat_dihedral[1], flat_dihedral[2], flat_dihedral[3], flat_dihedral[4], flat_dihedral[5]])
        if(each[2] in temp2 and each[1] not in temp1):
            adddihedral.append([temp1[temp2.index(each[2])], each[2], each[1], each[0], str(flat_dihedral_type), flat_dihedral[0], flat_dihedral[1], flat_dihedral[2], flat_dihedral[3], flat_dihedral[4], flat_dihedral[5]])

    bond.extend(addbond)
    angle.extend(addangle)
    dihedral.extend(adddihedral)
    # add capping hydrogen to FF
    #if (redist_flag=='yes'):
    ''' if polymer is capped with C3 super atom
    atomtype.append(['CU_cap', 'C', '1', '15.03500', '0.000', 'A', '3.96000e-01', '6.06680e-01'])
    atom.append([str(atoms_monomer+1), 'CU_cap', '1', 'END', 'C', '1', '0.00', '15.03500'])
    atom.append([str(atoms_monomer+2), 'CU_cap', '1', 'END', 'C', '1', '0.00', '15.03500'])
    bond.append([str(inter_bond_info[0][0]),str(atoms_monomer+1), str(bond_type), str(round(C_C_bond_length*0.1,4)), str(bond_stiffness) ])  
    bond.append([str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(atoms_monomer+2), str(bond_type), str(round(C_C_bond_length*0.1,4)), str(bond_stiffness) ])        
    angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-1),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[inter_bond_info[0][0]-1-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-2-1]))/2,2)),str(angle_stiffness)])
    angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-2),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[inter_bond_info[0][0]-1-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-2-1]))/2,2)),str(angle_stiffness)])
    angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2-1]))/2,2)),str(angle_stiffness)])
    angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2-1]))/2,2)),str(angle_stiffness)])
    improper.append([str(inter_bond_info[0][0]),str(atoms_monomer+1),str(inter_bond_info[0][0]-1),str(inter_bond_info[0][0]-2), str(improper_type), str(improper_value), str(improper_stiffness), str(multiplicity)])
    improper.append([str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(atoms_monomer+2),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2), str(improper_type), str(improper_value), str(improper_stiffness), str(multiplicity)])
    '''
    #else: if polymer is capped with H atom
    atomtype.append(['h_cap', 'H', '1', '10.000', '0.000', 'A', '0.00000e-01', '0.00000e-01'])
    atom.append([str(atoms_monomer+1), 'h_cap', '1', 'END', 'H', '1', '0.00', '10.000'])
    atom.append([str(atoms_monomer+2), 'h_cap', '1', 'END', 'H', '1', '0.00', '10.000'])
    bond.append([str(inter_bond_info[0][0]),str(atoms_monomer+1), str(bond_type), str(round(H_C_bond_length*0.1,4)), str(bond_stiffness) ])  
    if (int(degree)>1):
        bond.append([str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(atoms_monomer+2), str(bond_type), str(round(H_C_bond_length*0.1,4)), str(bond_stiffness) ])        
    else:
        bond.append([str(inter_bond_info[-1][-1]),str(atoms_monomer+2), str(bond_type), str(round(H_C_bond_length*0.1,4)), str(bond_stiffness) ])        

    '''smart but old version of calculating angles
    angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-1),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[inter_bond_info[0][0]-1-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-2-1]))/2,2)),str(angle_stiffness)])
    angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-2),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[inter_bond_info[0][0]-1-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-1],polymer_coordinate_noatomname[inter_bond_info[0][0]-2-1]))/2,2)),str(angle_stiffness)])
    angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2-1]))/2,2)),str(angle_stiffness)])
    angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2),str(angle_type),str(round((360-angle_calculator(polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]-1],polymer_coordinate_noatomname[np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2-1]))/2,2)),str(angle_stiffness)])
    '''
    if (fragmentlist[0] not in sixringconnect):
        angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-1),str(angle_type),str(126.75),str(angle_stiffness)])
        if (fragmentlist[0] not in fiveringmembers):
            angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-2),str(angle_type),str(126.75),str(angle_stiffness)])
        else:
            angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-4),str(angle_type),str(126.75),str(angle_stiffness)])
    else:
        angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-1),str(angle_type),str(120.00),str(angle_stiffness)])
        angle.append([str(atoms_monomer+1), str(inter_bond_info[0][0]),str(inter_bond_info[0][0]-2),str(angle_type),str(120.00),str(angle_stiffness)])

    if (fragmentlist[-1] not in sixringconnect):
        if (int(degree)>1):
            angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1),str(angle_type),str(126.75),str(angle_stiffness)])
            angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2),str(angle_type),str(126.75),str(angle_stiffness)])
        else:
            angle.append([str(atoms_monomer+2), str(inter_bond_info[-1][-1]),str(inter_bond_info[-1][-1]+1),str(angle_type),str(126.75),str(angle_stiffness)])
            angle.append([str(atoms_monomer+2), str(inter_bond_info[-1][-1]),str(inter_bond_info[-1][-1]+2),str(angle_type),str(126.75),str(angle_stiffness)])
    else:
        if (int(degree)>1):
            angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1),str(angle_type),str(120.00),str(angle_stiffness)])
            angle.append([str(atoms_monomer+2), str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2),str(angle_type),str(120.00),str(angle_stiffness)])
        else:
            angle.append([str(atoms_monomer+2), str(inter_bond_info[-1][-1]),str(inter_bond_info[-1][-1]+1),str(angle_type),str(120.00),str(angle_stiffness)])
            angle.append([str(atoms_monomer+2), str(inter_bond_info[-1][-1]),str(inter_bond_info[-1][-1]+2),str(angle_type),str(120.00),str(angle_stiffness)])
    improper.append([str(inter_bond_info[0][0]),str(atoms_monomer+1),str(inter_bond_info[0][0]-1),str(inter_bond_info[0][0]-2), str(improper_type), str(improper_value), str(improper_stiffness), str(multiplicity)])
    if (int(degree)>1):
        improper.append([str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]),str(atoms_monomer+2),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+1),str(np.cumsum(number_atom_total)[-2]+inter_bond_info[-1][-1]+2), str(improper_type), str(improper_value), str(improper_stiffness), str(multiplicity)])
    else:
        improper.append([str(inter_bond_info[-1][-1]),str(atoms_monomer+2),str(inter_bond_info[-1][-1]+1),str(inter_bond_info[-1][-1]+2), str(improper_type), str(improper_value), str(improper_stiffness), str(multiplicity)])
    #coordinates_concat=np.append(coordinates_atomname_numpy, coordinates_noatomname_numpy_moved_rotz,axis=1)
    XYZ=[]#writing coordinate info in correct xyz format
    for each in coordinates_concat:
        XYZ.append('{0: <12}'.format(each[0])+'{0: <12}'.format(round(each[1],5))+'{0: <12}'.format(round(each[2],5))+'{0: <12}'.format(round(each[3],5)))
    with open(path_output_oeff+moleculename+'_RU_SC_'+str(degree)+'.xyz', 'w') as f:
        f.write(str(len(XYZ))+"\n\n")
        f.writelines("%s\n" % l for l in XYZ)

    atomtype_=[]#writing atomtypes info in correct GROMACS itp format
    for each in atomtype:
        atomtypeline ='{0: <12}'.format(each[0])+'{0: <6}'.format(each[1])+'{0: <6}'.format(each[2])+'{0: <12}'.format(each[3])+'{0: <12}'.format(each[4])+'{0: <6}'.format(each[5])+'{0: >15}'.format(each[6])+'{0: >15}'.format(each[7])
        atomtype_.append(atomtypeline)
    atom_=[]#writing atom info in correct GROMACS itp format
    for each in atom:
        atom_line ='{0: <6}'.format(each[0])+'{0: <12}'.format(each[1])+'{0: <6}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <6}'.format(each[5])+'{0: >12}'.format(each[6])+'{0: >12}'.format(each[7])
        atom_.append(atom_line)
    bond_=[]#writing bond info in correct GROMACS itp format
    for each in bond:
        bondline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <6}'.format(each[2])+'{0: <10}'.format(each[3])+'{0: <10}'.format(each[4])
        bond_.append(bondline)
    pair_=[]#writing bond info in correct GROMACS itp format
    for each in pair:
        pairline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <6}'.format(each[2])
        pair_.append(pairline)
    angle_=[]#writing angle info in correct GROMACS itp format
    for each in angle:
        angleline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <6}'.format(each[3])+'{0: <10}'.format(each[4])+'{0: <10}'.format(each[5])
        angle_.append(angleline)
    dihedral_=[]#writing angle info in correct GROMACS itp format
    for each in dihedral:
        if (len(each)==8):
            dihedralline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <8}'.format(each[5])+'{0: <8}'.format(each[6])+'{0: <8}'.format(each[7])
            dihedral_.append(dihedralline)
        elif (len(each)==7):
            dihedralline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <8}'.format(each[5])+'{0: <8}'.format(each[6])
            dihedral_.append(dihedralline)
        elif (len(each)==11):
            dihedralline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <8}'.format(each[5])+'{0: <8}'.format(each[6])+'{0: <8}'.format(each[7])+'{0: <8}'.format(each[8])+'{0: <8}'.format(each[9])+'{0: <8}'.format(each[10])
            dihedral_.append(dihedralline)
    improper_=[]
    for each in improper:
        improperline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <8}'.format(each[5])+'{0: <8}'.format(each[6])+'{0: <8}'.format(each[7])
        improper_.append(improperline)
    with open(path_output_oeff+moleculename+'_RU_SC_'+str(degree)+'.itp', 'w') as f:
        f.write(";\n; GENERATED BY GAMMPS\n; Troisi Lab @ University of Liverpool\n")
        f.write("[ atomtypes ]\n")
        f.writelines("%s\n" % l for l in atomtype_)
        f.write("[ moleculetype ]\n;Name     nrexcl\n")
        f.writelines(moleculename+'_RU_SC_'+str(degree)+'     3\n')
        f.write("[ atoms ]\n")
        f.writelines("%s\n" % l for l in atom_)
        f.write("[ bonds ]\n")
        f.writelines("%s\n" % l for l in bond_)
        f.write("[ pairs ]\n")
        f.writelines("%s\n" % l for l in pair_) #only for side chains
        f.write("[ angles ]\n")
        f.writelines("%s\n" % l for l in angle_)
        f.write("[ dihedrals ]\n")
        f.writelines("%s\n" % l for l in dihedral_)
        f.write("[ dihedrals ]\n")
        f.writelines("%s\n" % l for l in improper_)
    if (int(degree)==1):
        with open(path_output_oeff+moleculename+'_RU_SC_'+str(degree)+'_soup.itp', 'w') as f:
            f.write(";\n; GENERATED BY GAMMPS\n; Troisi Lab @ University of Liverpool\n")
            f.write("[ moleculetype ]\n;Name     nrexcl\n")
            f.writelines(moleculename+'_RU_SC_'+str(degree)+'     3\n')
            f.write("[ atoms ]\n")
            f.writelines("%s\n" % l for l in atom_)
            f.write("[ bonds ]\n")
            f.writelines("%s\n" % l for l in bond_)
            f.write("[ pairs ]\n")
            f.writelines("%s\n" % l for l in pair_) #only for side chains
            f.write("[ angles ]\n")
            f.writelines("%s\n" % l for l in angle_)
            f.write("[ dihedrals ]\n")
            f.writelines("%s\n" % l for l in dihedral_)
            f.write("[ dihedrals ]\n")
            f.writelines("%s\n" % l for l in improper_)
    with open(path_output_oeff+moleculename+'_RU_SC_'+str(degree)+'.top', "w") as f:
        f.write("[ defaults ]\n")
        f.write('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
        f.write('1     3     yes     0.5     0.5\n\n')
        f.write('#include "'+moleculename+'_RU_SC_'+str(degree)+'.itp"\n\n')
        f.write('[ system ]\n; Name\n'+moleculename+'_RU_SC_'+str(degree)+'\n\n')
        f.write('[ molecules ]\n'+moleculename+'_RU_SC_'+str(degree)+'  1\n')
    
end = time.time()
with open(path_output+'.log', "a") as f:
    f.write("\n\nRuntime of this step was "+str(round(end - start,1))+" s\n")



for degree in degree_polymer:
    with open(path_output_oeff+'/polymer.tcl', "w") as f:
        f.write('mol new '+moleculename+'_RU_SC_'+degree+'.xyz type {xyz}\n\nmol modstyle 0 0 CPK 1.000000 0.300000 120.000000 120.000000\n\nmol modcolor 0 0 Name\n\ndisplay resetview\n\naxes location Off\n\n')
        f.write("color Display Background white\n\ndisplay projection orthographic\n\n"+"render Tachyon "+moleculename+"_RU_SC_"+degree+" "+VMD_RENDER+" -aasamples 12 %s -format TARGA -o %s.tga -res 1000 1000"+"\nexit\n")
    subprocess.run(['vmd','-dispdev','text','-e', 'polymer.tcl'], cwd = path_output_oeff)
    subprocess.run(['magick', moleculename+'_RU_SC_'+degree+'.tga', moleculename+'_RU_SC_'+degree+'.jpeg'], cwd = path_output_oeff)

    with open(path_output+'.html', "a") as f:
        f.write("\nPolymer models (coordinate and force field files) with DP = "+degree+" were made  <b><br>\n")
        f.write('<aside class="figures">\n')
        f.write('<figure>\n<img src="'+moleculename+'/'+moleculename+'_RU_SC_'+degree+'.jpeg" style="width: 500px; height: 300px; object-fit: cover;" alt="polymer">\n<figcaption>polymer with '+degree+' repeat units</figcaption>\n</figure>')
        f.write('</aside>\n')
with open(path_output+'.html', "a") as f:
    f.write("\n\nRuntime of this step was "+str(round(end - start,1))+" s<br><br>\n")   

