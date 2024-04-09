#!/usr/bin/env python
# coding: utf-8

'''
This is the second code to run in this step.

This code generates the coordinate (.xyz) and force field (.itp and .top according to GROMACS file formats) files
for the repeat unit structure with the sidechains specified in parameters.py.

'''


import collections
import numpy as np
import copy as cp
import subprocess
from matplotlib import pyplot as plt
import os
import time
import re


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
#function to copy all lines between two specific strings
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
#function to find the rotation matrix for aligning two vectors
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

start = time.time()
from parameters import *

parent_dir=cp.deepcopy(PATH_OUTPUT) #the Gaussian log files should be generated in the same folder as oligomer.xyz and .com files
moleculename = OLIGOMERNAME
path_output = os.path.join(parent_dir, moleculename)
path_output_oeff=os.path.join(parent_dir, moleculename+'/') # where the force field files will be coppied

#coordinate of repeat unit without side chain
with open(path_output_oeff+moleculename+'_RU.xyz', 'r') as f:
    coordinates=[]
    coordinates_noatomname=[]
    coordinates_atomname=[]
    coords=f.readlines()[2:]
    for coord in coords:
        coordinates.append(coord.split())
        coordinates_noatomname.append([float(coord.split()[1]), float(coord.split()[2]),float(coord.split()[3])])
        coordinates_atomname.append([coord.split()[0]])
#end of coordinate of repeat unit without side chain

#FF parameters of repeat unit without side chain
with open(path_output_oeff+moleculename+'_RU.itp', 'r') as f:
    lines=f.readlines()
    atomtypes_noside=inbetweenlines('[ atomtypes ]', '[ moleculetype ]')
    atomtypes_noside_copy=cp.deepcopy(atomtypes_noside)
    atomlist_noside=inbetweenlines('[ atoms ]', '[ bonds ]')
    atomlist_noside_copy=cp.deepcopy(atomlist_noside)
    bondlist_noside=inbetweenlines('[ bonds ]', '[ angles ]')
    bondlist_noside_copy=cp.deepcopy(bondlist_noside)
    anglelist_noside=inbetweenlines('[ angles ]', '[ dihedrals ]')
    anglelist_noside_copy=cp.deepcopy(anglelist_noside)
    dihedrallist_noside=inbetweenlines('[ dihedrals ]', '')
    dihedrallist_noside_copy=cp.deepcopy(dihedrallist_noside)
#end of FF parameters of repeat unit without side chain

#finding the carbon atoms (and all connecting hydrogens)to which the side chain should be connected (the teminal carbons wichi have 3 bonds with H)
temp=[]
for each in bondlist_noside:
    temp.append(each[0])
    temp.append(each[1])  
atoms_fourbonds=[item for item, count in collections.Counter(temp).items() if count > 3]
temp=[]
tempnum=[]
for count, each in enumerate(atoms_fourbonds):
    temp.append([int(each)])
    tempnum.append([int(each)])
    for item in bondlist_noside:
        if (each==item[0]):
            temp[count].append(atomlist_noside[int(item[1])-1][4])
            tempnum[count].append(int(item[1]))
        elif (each==item[1]):
            temp[count].append(atomlist_noside[int(item[0])-1][4])
            tempnum[count].append(int(item[0]))

carbon_connect_side=[]
hydrogen_connected_carbon_connect=[]
heavy_connected_carbon_connect=[]
atoms_connected_carbon_connect=[]
for each in temp:
    if([item for item, count in collections.Counter(each).items() if (count == 3 and item=='H')]!=[]):
        carbon_connect_side.append(int(each[0]))
        for number in tempnum:
            if (number[0]==each[0]):
                atoms_connected_carbon_connect.append(number)
                hydrogen=[]
                heavy=[]
                hydrogen.append(int(number[0]))
                for every in number:
                    if(atomlist_noside[int(every)-1][4]=='H'):
                        hydrogen.append(every)
                    else:
                        heavy.append(every)
                hydrogen_connected_carbon_connect.append(hydrogen)
                heavy_connected_carbon_connect.append(heavy)
#end of finding the carbon (and all connecting hydrogens) atoms to which the side chain should be connected (atoms_connect_side contains all of them)

#OPLS values taken from original FF
H_charge=0.06
H_charge_end=0.04
C_charge=-0.120
C_UN_charge=0.000
O_UN_charge=-0.500
CO_UN_charge=0.250
O_UN_sigma=2.96000e-01
O_UN_epsilon=8.78640e-01
C2_UN_sigma=3.90500e-01
C3_UN_sigma=3.96000e-01
C2_UN_epsilon=4.93712e-01
C3_UN_epsilon=6.06680e-01
O_UN_mass=15.99940
C2_UN_mass=14.02700
C3_UN_mass=15.03500
C_mass=12.01100
opls_bondtype=1
opls_pairtype=1
opls_angletype=1
opls_dihedraltype=3
C_H_bondlength=1.0900 
C_C_bondlength=1.5290 
C_O_bondlength=1.4100   #glycol
C_H_bondstiff=284512.0 
C_C_bondstiff=224262.4
C_O_bondstiff=267776.0 #glycol
C_C_C_angle=112.7
C_C_C_anglestiff=488.3 
C_C_H_angle=110.7    
C_C_H_anglestiff=313.8 
O_C_H_angle=109.500     #glycol 
O_C_H_anglestiff=292.880#glycol
H_C_H_angle=107.8
H_C_H_anglestiff=276.1 
C_O_C_angle=109.500     #glycol
C_O_C_anglestiff=502.080#glycol
C_C_O_angle=109.500     #glycol
C_C_O_anglestiff=418.400#glycol
H_C_C_H_dihedralpotential='0.627   1.883   0.000  -2.510   0.000   0.000'
C_C_C_C_dihedralpotential='2.929  -1.464   0.209  -1.674   0.000   0.000'
C_C_C_H_dihedralpotential='0.627   1.883   0.000  -2.510   0.000   0.000'
C_C_C_O_dihedralpotential='2.874   0.582   2.092  -5.548   0.000   0.000'#glycol
C_O_C_C_dihedralpotential='1.715   2.845   1.046  -5.607   0.000   0.000'#glycol
C_O_C_H_dihedralpotential='1.590   4.770   0.000  -6.360   0.000   0.000'#glycol
O_C_C_O_dihedralpotential='-1.151  1.151   0.000   0.000   0.000   0.000'#glycol
C_C_O_H_dihedralpotential='0.979   2.937   0.000  -3.916   0.000   0.000'#glycol


C_connection=[] #save the id of C atoms to which sidechains are attached
H_renamed_C=[] #save the id of H atoms changed after sidechain attachment

addbond=[] #bonds to be added to itp files
addbond_temp=[] #for including the connecting carbon in the side chain
side_type=cp.deepcopy(SIDECHAIN_TYPE) #read the sidechain type infro from the parameter file
side_length=cp.deepcopy(SIDECHAIN_LENGTH) #read the sidechain length info from the parameter file
side_branch_point=cp.deepcopy(SIDECHAIN_BRANCH_POINT) #read the sidechain branching points info from the parameter file
side_branch_length=cp.deepcopy(SIDECHAIN_BRANCH_LENGTH) #read the sidechain branch length info from the parameter file
for count, each in enumerate(heavy_connected_carbon_connect):
    if (side_length[count]!=str(0) and side_type[count]=='a'):
        reference_coord=[coordinates_noatomname[int(each[0])-1][0],coordinates_noatomname[int(each[0])-1][1],coordinates_noatomname[int(each[0])-1][2]]
        coordinates_moved=np.array(coordinates_noatomname)-np.array(reference_coord) # move coordinates to have connecting carbon at 0 0 0
        coordinates_noatomname=cp.deepcopy(coordinates_moved)
        rot_mat = rotation_matrix_from_vectors(np.array([0.0,0.0,1.0]), coordinates_moved[int(each[1])-1]) # find rotation matrix to rotate along heavy-carbon bond
        coordinates_noatomname=np.matmul(coordinates_moved,rot_mat) # rotate along heavy-carbon bond
        coordinates_noatomname[hydrogen_connected_carbon_connect[count][1]-1]= np.array([0.0, -C_H_bondlength, 0.0]) # changing coordinate of first hydrogen
        coordinates_noatomname[hydrogen_connected_carbon_connect[count][2]-1]= np.array([0.0, C_H_bondlength, 0.0]) # changing coordinate of first hydrogen
        coordinates_noatomname[hydrogen_connected_carbon_connect[count][3]-1]= np.array([0.0, 0.0, -C_C_bondlength]) # changing coordinate of first hydrogen (will be change to cabon)
        coordinates_atomname[hydrogen_connected_carbon_connect[count][3]-1]= ['C'] #change atom type in coordinate file
        for item in atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]: #changing atom information in FF [atomtypes]
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][1]='C'
            atomtypes_noside_copy[hydrogen_connected_carbon_connect[count][3]-1][1]='C'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][2]='6'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][3]=str(C2_UN_mass)
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][4]='0.000'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][5]= 'A'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][6]=str(C2_UN_sigma)
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][7]=str(C2_UN_epsilon)
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][4]='C'
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][7]=str(C2_UN_mass)
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        atomlist_noside[hydrogen_connected_carbon_connect[count][0]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        H_renamed_C.append(hydrogen_connected_carbon_connect[count][3])
        C_connection.append(hydrogen_connected_carbon_connect[count][0])

        #redistribute the modified hydrogen charge on the remaining two hydrogens
        atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][6]=str(round(float(atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][6])+float(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][6])/2,4))
        atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][6]=str(round(float(atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][6])+float(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][6])/2,4))
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][6]=str(C_UN_charge) # changinge atom information in FF [atoms]



        for item in bondlist_noside: # correcting bond information for H changed to C and the other two remaining H
            if (int(item[0])==hydrogen_connected_carbon_connect[count][3] or int(item[1])==hydrogen_connected_carbon_connect[count][3]):
                item[3]=str(round(C_C_bondlength*0.1,5)) #angstrom to nm
                item[4]=str(C_C_bondstiff) #opls
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2]):
                item[3]=str(round(C_H_bondlength*0.1,5)) #angstrom to nm
                item[4]=str(C_H_bondstiff) #opls
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][1] or int(item[1])==hydrogen_connected_carbon_connect[count][1]):
                item[3]=str(round(C_H_bondlength*0.1,5)) #angstrom to nm
                item[4]=str(C_H_bondstiff) #opls


        for item in anglelist_noside: # correcting angle information for H changed to C and the other two remaining H
            if (int(item[0])==hydrogen_connected_carbon_connect[count][3] or int(item[1])==hydrogen_connected_carbon_connect[count][3] or int(item[2])==hydrogen_connected_carbon_connect[count][3]):
                if(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_C_angle)
                    item[5]=str(C_C_C_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2] or int(item[2])==hydrogen_connected_carbon_connect[count][2]):
                if(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(H_C_H_angle)
                    item[5]=str(H_C_H_anglestiff)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][1] or int(item[1])==hydrogen_connected_carbon_connect[count][1] or int(item[2])==hydrogen_connected_carbon_connect[count][1]):
                if(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(H_C_H_angle)
                    item[5]=str(H_C_H_anglestiff)
        for dihed_count, item in enumerate(dihedrallist_noside): # deleting dihedral lines for the H changed to C and the other two remaining H (dihedral parameters with side chains will be created later)
            if (int(item[0])==hydrogen_connected_carbon_connect[count][3] or int(item[1])==hydrogen_connected_carbon_connect[count][3] or int(item[2])==hydrogen_connected_carbon_connect[count][3] or int(item[3])==hydrogen_connected_carbon_connect[count][3]):
                dihedrallist_noside.pop(dihed_count)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2] or int(item[2])==hydrogen_connected_carbon_connect[count][2] or int(item[3])==hydrogen_connected_carbon_connect[count][2]):
                dihedrallist_noside.pop(dihed_count)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][1] or int(item[1])==hydrogen_connected_carbon_connect[count][1] or int(item[2])==hydrogen_connected_carbon_connect[count][1] or int(item[3])==hydrogen_connected_carbon_connect[count][1]):
                dihedrallist_noside.pop(dihed_count)
        


        # include the connecting carbon in the side chain FF
   
        addbond_temp.append(str(heavy_connected_carbon_connect[count][0])+'   '+str(heavy_connected_carbon_connect[count][1])+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))
        addbond_temp.append(str(hydrogen_connected_carbon_connect[count][0])+'   '+str(hydrogen_connected_carbon_connect[count][3])+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))


        carbonID=hydrogen_connected_carbon_connect[count][3]#connecting carbon ID
        for i in range(int(side_length[count])):
            coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_C_bondlength*(i+2)]]),axis=0)

            coordinates_atomname.append(['C'])

            atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
            atomlist_noside[-1][0]=str(len(atomlist_noside))
            atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
            atomlist_noside[-1][4]='C'
            atomlist_noside[-1][6]=str(C_UN_charge)

            atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
            atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
            atomtypes_noside[-1][1]='C'
            addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))
            carbonID=len(atomtypes_noside)


            if (i==int(side_length[count])-1):
                atomtypes_noside[-1][3]=str(C3_UN_mass)
                atomtypes_noside[-1][6]=str(C3_UN_sigma)
                atomtypes_noside[-1][7]=str(C3_UN_epsilon)
                atomlist_noside[-1][7]=str(C3_UN_mass)
        #branching
        if (side_branch_point[count]!='0'):
            #if the branching point is not at the first carbon
            if (side_branch_point[count]=='1'):
                reference_coord=[coordinates_noatomname[int(each[0])-1][0],coordinates_noatomname[int(each[0])-1][1],coordinates_noatomname[int(each[0])-1][2]]
                coordinates_moved=np.array(coordinates_noatomname)-np.array(reference_coord) # move coordinates to have connecting carbon at 0 0 0
                coordinates_noatomname=cp.deepcopy(coordinates_moved)
                rot_mat = rotation_matrix_from_vectors(np.array([0.0,1.0,0.0]), coordinates_moved[int(each[1])-1]) # find rotation matrix to rotate along heavy-carbon bond
                coordinates_noatomname=np.matmul(coordinates_moved,rot_mat) # rotate along heavy-carbon bond
                coordinates_noatomname[hydrogen_connected_carbon_connect[count][2]-1]= np.array([0.0, 0.0, -C_C_bondlength]) # fixing the coordinate of the second hydrogen changing to carbon
                coordinates_atomname[hydrogen_connected_carbon_connect[count][2]-1]= ['C'] #change atom type in coordinate file
                for item in atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1]: #changing atom information in FF [atomtypes]
                    atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1][1]='C'
                    atomtypes_noside_copy[hydrogen_connected_carbon_connect[count][2]-1][1]='C'
                    atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1][2]='6'
                    atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1][3]=str(C2_UN_mass)
                    atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1][4]='0.000'
                    atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1][5]= 'A'
                    atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1][6]=str(C2_UN_sigma)
                    atomtypes_noside[hydrogen_connected_carbon_connect[count][2]-1][7]=str(C2_UN_epsilon)
                atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][4]='C'
                atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][7]=str(C2_UN_mass)
                H_renamed_C.append(hydrogen_connected_carbon_connect[count][2])

                for item in bondlist_noside: # correcting bond information for H changed to C and the other two remaining H
                    if (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2]):
                        item[3]=str(round(C_C_bondlength*0.1,5)) #angstrom to nm
                        item[4]=str(C_C_bondstiff) #opls

                for item in anglelist_noside: # correcting angle information for H changed to C and the other two remaining H
                    if (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2] or int(item[2])==hydrogen_connected_carbon_connect[count][2]):
                        if(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                            item[4]=str(C_C_C_angle)
                            item[5]=str(C_C_C_anglestiff)
                        elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                            item[4]=str(C_C_C_angle)
                            item[5]=str(C_C_C_anglestiff)
                        elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                            item[4]=str(C_C_H_angle)
                            item[5]=str(C_C_H_anglestiff)

                #add two hydrogen to the new carbon (the H renamed to C)/ FF params


                # include the connecting carbon in the side chain FF
           
                addbond_temp.append(str(hydrogen_connected_carbon_connect[count][0])+'   '+str(hydrogen_connected_carbon_connect[count][2])+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))

                carbonID=hydrogen_connected_carbon_connect[count][2]#connecting carbon ID
                for i in range(int(side_branch_length[count])):
                    coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_C_bondlength*(i+2)]]),axis=0)

                    coordinates_atomname.append(['C'])

                    atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                    atomlist_noside[-1][0]=str(len(atomlist_noside))
                    atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                    atomlist_noside[-1][4]='C'
                    atomlist_noside[-1][6]=str(C_UN_charge)

                    atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                    atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                    atomtypes_noside[-1][1]='C'
                    addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))
                    carbonID=len(atomtypes_noside)


                    if (i==int(side_length[count])-1):
                        atomtypes_noside[-1][3]=str(C3_UN_mass)
                        atomtypes_noside[-1][6]=str(C3_UN_sigma)
                        atomtypes_noside[-1][7]=str(C3_UN_epsilon)
                        atomlist_noside[-1][7]=str(C3_UN_mass)
            #if the branching point is at the second carbon
            elif (side_branch_point[count]=='2'):
                reference_coord=[coordinates_noatomname[hydrogen_connected_carbon_connect[count][3]-1][0],coordinates_noatomname[hydrogen_connected_carbon_connect[count][3]-1][1],coordinates_noatomname[hydrogen_connected_carbon_connect[count][3]-1][2]]
                coordinates_moved=np.array(coordinates_noatomname)-np.array(reference_coord) # move coordinates to have connecting carbon at 0 0 0
                coordinates_noatomname=cp.deepcopy(coordinates_moved)
                rot_mat = rotation_matrix_from_vectors(np.array([0.0,1.0,0.0]), coordinates_moved[int(each[1])-1]) # find rotation matrix to rotate along heavy-carbon bond
                coordinates_noatomname=np.matmul(coordinates_moved,rot_mat) # rotate along heavy-carbon bond
           
                addbond_temp.append(str(hydrogen_connected_carbon_connect[count][3])+'   '+str(len(atomtypes_noside)+1)+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))

                carbonID=hydrogen_connected_carbon_connect[count][3]#connecting carbon ID
                for i in range(int(side_branch_length[count])):
                    coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_C_bondlength*(i+1)]]),axis=0)

                    coordinates_atomname.append(['C'])

                    atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                    atomlist_noside[-1][0]=str(len(atomlist_noside))
                    atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                    atomlist_noside[-1][4]='C'
                    atomlist_noside[-1][6]=str(C_UN_charge)

                    atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                    atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                    atomtypes_noside[-1][1]='C'
                    addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))
                    carbonID=len(atomtypes_noside)


                    if (i==int(side_length[count])-1):
                        atomtypes_noside[-1][3]=str(C3_UN_mass)
                        atomtypes_noside[-1][6]=str(C3_UN_sigma)
                        atomtypes_noside[-1][7]=str(C3_UN_epsilon)
                        atomlist_noside[-1][7]=str(C3_UN_mass)
            #if the branching point is at the third or further carbon
            else:

                reference_coord=[coordinates_noatomname[carbonID-int(side_length[count])+int(side_branch_point[count])-2-1][0],coordinates_noatomname[carbonID-int(side_length[count])+int(side_branch_point[count])-2-1][1],coordinates_noatomname[carbonID-int(side_length[count])+int(side_branch_point[count])-2-1][2]]
                coordinates_moved=np.array(coordinates_noatomname)-np.array(reference_coord) # move coordinates to have connecting carbon at 0 0 0
                coordinates_noatomname=cp.deepcopy(coordinates_moved)
                rot_mat = rotation_matrix_from_vectors(np.array([0.0,1.0,0.0]), coordinates_moved[int(each[1])-1]) # find rotation matrix to rotate along heavy-carbon bond
                coordinates_noatomname=np.matmul(coordinates_moved,rot_mat) # rotate along heavy-carbon bond
                
                addbond_temp.append(str(carbonID-int(side_length[count])+int(side_branch_point[count])-2)+'   '+str(len(atomtypes_noside)+1)+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))

                carbonID=carbonID-int(side_length[count])+int(side_branch_point[count])-2#connecting carbon ID
                for i in range(int(side_branch_length[count])):
                    coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_C_bondlength*(i+1)]]),axis=0)

                    coordinates_atomname.append(['C'])

                    atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                    atomlist_noside[-1][0]=str(len(atomlist_noside))
                    atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                    atomlist_noside[-1][4]='C'
                    atomlist_noside[-1][6]=str(C_UN_charge)

                    atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                    atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                    atomtypes_noside[-1][1]='C'
                    addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))
                    carbonID=len(atomtypes_noside)


                    if (i==int(side_length[count])-1):
                        atomtypes_noside[-1][3]=str(C3_UN_mass)
                        atomtypes_noside[-1][6]=str(C3_UN_sigma)
                        atomtypes_noside[-1][7]=str(C3_UN_epsilon)
                        atomlist_noside[-1][7]=str(C3_UN_mass)

    
    elif (side_length[count]!=str(0) and side_type[count]=='eg'):
        reference_coord=[coordinates_noatomname[int(each[0])-1][0],coordinates_noatomname[int(each[0])-1][1],coordinates_noatomname[int(each[0])-1][2]]
        coordinates_moved=np.array(coordinates_noatomname)-np.array(reference_coord) # move coordinates to have connecting carbon at 0 0 0
        coordinates_noatomname=cp.deepcopy(coordinates_moved)
        rot_mat = rotation_matrix_from_vectors(np.array([0.0,0.0,1.0]), coordinates_moved[int(each[1])-1]) # find rotation matrix to rotate along heavy-carbon bond
        coordinates_noatomname=np.matmul(coordinates_moved,rot_mat) # rotate along heavy-carbon bond
        coordinates_noatomname[hydrogen_connected_carbon_connect[count][1]-1]= np.array([0.0, -C_H_bondlength, 0.0]) # changing coordinate of first hydrogen
        coordinates_noatomname[hydrogen_connected_carbon_connect[count][2]-1]= np.array([0.0, C_H_bondlength, 0.0]) # changing coordinate of first hydrogen
        coordinates_noatomname[hydrogen_connected_carbon_connect[count][3]-1]= np.array([0.0, 0.0, -C_C_bondlength]) # changing coordinate of first hydrogen (will be change to cabon)
        coordinates_atomname[hydrogen_connected_carbon_connect[count][3]-1]= ['C'] #change atom type in coordinate file
        for item in atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]: #changing atom information in FF [atomtypes]
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][1]='C'
            atomtypes_noside_copy[hydrogen_connected_carbon_connect[count][3]-1][1]='C'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][2]='6'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][3]=str(C2_UN_mass)
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][4]='0.000'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][5]= 'A'
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][6]=str(C2_UN_sigma)
            atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1][7]=str(C2_UN_epsilon)
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][4]='C'
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][7]=str(C2_UN_mass)
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        atomlist_noside[hydrogen_connected_carbon_connect[count][0]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][3]='s'+atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][3] # changinge residue name to side chain (s...) in FF [atoms]
        H_renamed_C.append(hydrogen_connected_carbon_connect[count][3])
        C_connection.append(hydrogen_connected_carbon_connect[count][0])

        #redistribute the modified hydrogen charge on the remaining two hydrogens
        atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][6]=str(round(float(atomlist_noside[hydrogen_connected_carbon_connect[count][1]-1][6])+float(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][6])/2,4))
        atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][6]=str(round(float(atomlist_noside[hydrogen_connected_carbon_connect[count][2]-1][6])+float(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][6])/2,4))
        atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1][6]=str(CO_UN_charge) # changinge atom information in FF [atoms]



        for item in bondlist_noside: # correcting bond information for H changed to C and the other two remaining H
            if (int(item[0])==hydrogen_connected_carbon_connect[count][3] or int(item[1])==hydrogen_connected_carbon_connect[count][3]):
                item[3]=str(round(C_C_bondlength*0.1,5)) #angstrom to nm
                item[4]=str(C_C_bondstiff) #opls
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2]):
                item[3]=str(round(C_H_bondlength*0.1,5)) #angstrom to nm
                item[4]=str(C_H_bondstiff) #opls
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][1] or int(item[1])==hydrogen_connected_carbon_connect[count][1]):
                item[3]=str(round(C_H_bondlength*0.1,5)) #angstrom to nm
                item[4]=str(C_H_bondstiff) #opls


        for item in anglelist_noside: # correcting angle information for H changed to C and the other two remaining H, here is different as compared to alkyl sidechain 
            if (int(item[0])==hydrogen_connected_carbon_connect[count][3] or int(item[1])==hydrogen_connected_carbon_connect[count][3] or int(item[2])==hydrogen_connected_carbon_connect[count][3]):
                if(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_C_angle)
                    item[5]=str(C_C_C_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='O'):
                    item[4]=str(C_C_O_angle)
                    item[5]=str(C_C_O_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='O' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_O_angle)
                    item[5]=str(C_C_O_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='O'):
                    item[4]=str(O_C_H_angle)
                    item[5]=str(O_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='O' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(O_C_H_angle)
                    item[5]=str(O_C_H_anglestiff)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2] or int(item[2])==hydrogen_connected_carbon_connect[count][2]):
                if(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_C_angle)
                    item[5]=str(C_C_C_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='O'):
                    item[4]=str(C_C_O_angle)
                    item[5]=str(C_C_O_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='O' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_O_angle)
                    item[5]=str(C_C_O_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='O'):
                    item[4]=str(O_C_H_angle)
                    item[5]=str(O_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='O' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(O_C_H_angle)
                    item[5]=str(O_C_H_anglestiff)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][1] or int(item[1])==hydrogen_connected_carbon_connect[count][1] or int(item[2])==hydrogen_connected_carbon_connect[count][1]):
                if(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_C_angle)
                    item[5]=str(C_C_C_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='O'):
                    item[4]=str(C_C_O_angle)
                    item[5]=str(C_C_O_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='O' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_O_angle)
                    item[5]=str(C_C_O_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='C'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='C' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(C_C_H_angle)
                    item[5]=str(C_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='H' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='O'):
                    item[4]=str(O_C_H_angle)
                    item[5]=str(O_C_H_anglestiff)
                elif(atomtypes_noside_copy[int(item[0])-1][1]=='O' and atomtypes_noside_copy[int(item[1])-1][1]=='C' and atomtypes_noside_copy[int(item[2])-1][1]=='H'):
                    item[4]=str(O_C_H_angle)
                    item[5]=str(O_C_H_anglestiff)
        for dihed_count, item in enumerate(dihedrallist_noside): # deleting dihedral lines for the H changed to C and the other two remaining H (dihedral parameters with side chains will be created later)
            if (int(item[0])==hydrogen_connected_carbon_connect[count][3] or int(item[1])==hydrogen_connected_carbon_connect[count][3] or int(item[2])==hydrogen_connected_carbon_connect[count][3] or int(item[3])==hydrogen_connected_carbon_connect[count][3]):
                dihedrallist_noside.pop(dihed_count)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][2] or int(item[1])==hydrogen_connected_carbon_connect[count][2] or int(item[2])==hydrogen_connected_carbon_connect[count][2] or int(item[3])==hydrogen_connected_carbon_connect[count][2]):
                dihedrallist_noside.pop(dihed_count)
            elif (int(item[0])==hydrogen_connected_carbon_connect[count][1] or int(item[1])==hydrogen_connected_carbon_connect[count][1] or int(item[2])==hydrogen_connected_carbon_connect[count][1] or int(item[3])==hydrogen_connected_carbon_connect[count][1]):
                dihedrallist_noside.pop(dihed_count)
        #add two hydrogen to the new carbon (the H renamed to C)/ FF params


        # include the connecting carbon in the side chain FF
   
        addbond_temp.append(str(heavy_connected_carbon_connect[count][0])+'   '+str(heavy_connected_carbon_connect[count][1])+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))

        carbonID=hydrogen_connected_carbon_connect[count][3]#connecting carbon ID
        for i in range(int(side_length[count])):
            #first, carbon
            if (i==0):
                coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -(C_O_bondlength+C_C_bondlength)-(2*C_O_bondlength+C_C_bondlength)*(i)]]),axis=0)

                coordinates_atomname.append(['O'])

                atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomlist_noside[-1][0]=str(len(atomlist_noside))
                atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                atomlist_noside[-1][4]='O'
                atomlist_noside[-1][6]=str(O_UN_charge)

                atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                atomtypes_noside[-1][1]='O'
                addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_O_bondlength*0.1,4))+'   '+str(C_O_bondstiff))
                carbonID=len(atomtypes_noside)

                #second, oxygen
                coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_O_bondlength-(C_O_bondlength+C_C_bondlength)-(2*C_O_bondlength+C_C_bondlength)*(i)]]),axis=0)

                coordinates_atomname.append(['C'])

                atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomlist_noside[-1][0]=str(len(atomlist_noside))
                atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                atomlist_noside[-1][4]='C'
                atomlist_noside[-1][6]=str(CO_UN_charge)

                atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                atomtypes_noside[-1][1]='C'
                addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_O_bondlength*0.1,4))+'   '+str(C_O_bondstiff))
                carbonID=len(atomtypes_noside)

            elif(i>0):
                coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_C_bondlength-(2*C_O_bondlength+C_C_bondlength)*(i)]]),axis=0)

                coordinates_atomname.append(['C'])

                atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomlist_noside[-1][0]=str(len(atomlist_noside))
                atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                atomlist_noside[-1][4]='C'
                atomlist_noside[-1][6]=str(CO_UN_charge)

                atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                atomtypes_noside[-1][1]='C'
                addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_C_bondlength*0.1,4))+'   '+str(C_C_bondstiff))
                carbonID=len(atomtypes_noside)

                #second, oxygen
                coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_O_bondlength-C_C_bondlength-(2*C_O_bondlength+C_C_bondlength)*(i)]]),axis=0)

                coordinates_atomname.append(['O'])

                atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomlist_noside[-1][0]=str(len(atomlist_noside))
                atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                atomlist_noside[-1][4]='O'
                atomlist_noside[-1][6]=str(O_UN_charge)

                atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                atomtypes_noside[-1][1]='O'
                addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_O_bondlength*0.1,4))+'   '+str(C_O_bondstiff))
                carbonID=len(atomtypes_noside)

                #third, carbon
                coordinates_noatomname=np.append(coordinates_noatomname,np.array([[0.0, 0.0, -C_O_bondlength-C_C_bondlength-C_O_bondlength-(2*C_O_bondlength+C_C_bondlength)*(i)]]),axis=0)

                coordinates_atomname.append(['C'])

                atomlist_noside.append(cp.deepcopy(atomlist_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomlist_noside[-1][0]=str(len(atomlist_noside))
                atomlist_noside[-1][1]='oeff_'+str(len(atomlist_noside))
                atomlist_noside[-1][4]='C'
                atomlist_noside[-1][6]=str(CO_UN_charge)

                atomtypes_noside.append(cp.deepcopy(atomtypes_noside[hydrogen_connected_carbon_connect[count][3]-1]))
                atomtypes_noside[-1][0]='oeff_'+str(len(atomtypes_noside))
                atomtypes_noside[-1][1]='C'
                addbond.append(str(len(atomtypes_noside))+'   '+str(carbonID)+'   '+str(opls_bondtype)+'   '+str(round(C_O_bondlength*0.1,4))+'   '+str(C_O_bondstiff))
                carbonID=len(atomtypes_noside)
            if (i==int(side_length[count])-1):
                atomtypes_noside[-1][3]=str(C3_UN_mass)
                atomtypes_noside[-1][6]=str(C3_UN_sigma)
                atomtypes_noside[-1][7]=str(C3_UN_epsilon)
                atomlist_noside[-1][7]=str(C3_UN_mass)
        # no branching for glycol chains are considered
        
                   
#add pairs, angles, and dihedrals for side chains
addbond_temp.extend(cp.deepcopy(addbond))
pairs=[s + '{0: <6}'.format(str(opls_pairtype)) for s in angle_dihedral(addbond_temp)[0]]
angles=[s + '{0: <6}'.format(str(opls_angletype)) for s in angle_dihedral(addbond_temp)[1]]
dihedrals=[s + '{0: <6}'.format(str(opls_dihedraltype)) for s in angle_dihedral(addbond_temp)[2]]

for i, each in enumerate(angles):
    if(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C'):
        angle_value=C_C_C_angle
        angle_stiff=C_C_C_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='N'):
        angle_value=C_C_C_angle
        angle_stiff=C_C_C_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='N' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C'):
        angle_value=C_C_C_angle
        angle_stiff=C_C_C_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='H'):
        angle_value=C_C_H_angle
        angle_stiff=C_C_H_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='H' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C'):
        angle_value=C_C_H_angle
        angle_stiff=C_C_H_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='H' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='H'):
        angle_value=H_C_H_angle
        angle_stiff=H_C_H_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='O'):
        angle_value=C_C_O_angle
        angle_stiff=C_C_O_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='O' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C'):
        angle_value=C_C_O_angle
        angle_stiff=C_C_O_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='O' and atomtypes_noside[int(each.split()[2])-1][1]=='C'):
        angle_value=C_O_C_angle
        angle_stiff=C_O_C_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='H' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='O'):
        angle_value=O_C_H_angle
        angle_stiff=O_C_H_anglestiff
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='O' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='H'):
        angle_value=O_C_H_angle
        angle_stiff=O_C_H_anglestiff
    angles[i]=each+'{0: <8}'.format(str(angle_value))+'{0: <8}'.format(str(angle_stiff))
for i, each in enumerate(dihedrals):
    if(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='C'):
        dihedral_value=C_C_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='N' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='C'):
        dihedral_value=C_C_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='N'):
        dihedral_value=C_C_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='H'):
        dihedral_value=C_C_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='H' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='C'):
        dihedral_value=C_C_C_H_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='H'):
        dihedral_value=C_C_C_H_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='H' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='H'):
        dihedral_value=H_C_C_H_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='O' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='C'):
        dihedral_value=O_C_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='O'):
        dihedral_value=O_C_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='O' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='O'):
        dihedral_value=O_C_C_O_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='O' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='C'):
        dihedral_value=C_O_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='C' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='O' and atomtypes_noside[int(each.split()[3])-1][1]=='C'):
        dihedral_value=C_O_C_C_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='O' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='H'):
        dihedral_value=C_C_O_H_dihedralpotential
    elif(atomtypes_noside[int(each.split()[0])-1][1]=='H' and atomtypes_noside[int(each.split()[1])-1][1]=='C' and atomtypes_noside[int(each.split()[2])-1][1]=='C' and atomtypes_noside[int(each.split()[3])-1][1]=='O'):
        dihedral_value=C_C_O_H_dihedralpotential
    dihedrals[i]=each+'{0: <8}'.format(str(dihedral_value))
    
# mass redistribution
redist_flag=cp.deepcopy(MASSDIST)
if (redist_flag=='yes'):
    excess_mass=0.0
    h_number=0
    for each in atomtypes_noside:
        if(each[1]=='C' and float(each[3])>C_mass):
            excess_mass+=float(each[3])-C_mass
        if(each[1]=='H'):
            h_number+=1
    each_h_mass=excess_mass/h_number
    for each in atomtypes_noside:
        if(each[1]=='C' and float(each[3])>C_mass):
            each[3]=str(C_mass)
        if(each[1]=='H'):
            each[3]=str(round(each_h_mass+float(each[3]),5))
    for each in atomlist_noside:
        if(each[4]=='C' and float(each[7])>C_mass):
            each[7]=str(C_mass)
        if(each[4]=='H'):
            each[7]=str(round(each_h_mass+float(each[7]),5))
coordinates_concat=np.append(np.array(coordinates_atomname,dtype=object),coordinates_noatomname,axis=1)

XYZ=[]#writing coordinate info in correct xyz format
for each in coordinates_concat:
    XYZ.append('{0: <12}'.format(each[0])+'{0: <12}'.format(round(each[1],5))+'{0: <12}'.format(round(each[2],5))+'{0: <12}'.format(round(each[3],5)))
atomtype=[]#writing atomtypes info in correct GROMACS itp format
for each in atomtypes_noside:
    atomtypeline ='{0: <12}'.format(each[0])+'{0: <6}'.format(each[1])+'{0: <6}'.format(each[2])+'{0: <12}'.format(each[3])+'{0: <12}'.format(each[4])+'{0: <6}'.format(each[5])+'{0: >15}'.format(each[6])+'{0: >15}'.format(each[7])
    atomtype.append(atomtypeline)
atom=[]#writing atom info in correct GROMACS itp format
for each in atomlist_noside:
    atom_line ='{0: <6}'.format(each[0])+'{0: <12}'.format(each[1])+'{0: <6}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <6}'.format(each[4])+'{0: <6}'.format(each[5])+'{0: >12}'.format(each[6])+'{0: >12}'.format(each[7])
    atom.append(atom_line)
bond=[]#writing bond info in correct GROMACS itp format
for each in addbond:
    bondlist_noside.extend([each.split()])
for each in bondlist_noside:
    bondline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <6}'.format(each[2])+'{0: <10}'.format(each[3])+'{0: <10}'.format(each[4])
    bond.append(bondline)
angle=[]#writing angle info in correct GROMACS itp format
for each in angles:
    anglelist_noside.extend([each.split()])
for each in anglelist_noside:
    angleline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <6}'.format(each[3])+'{0: <10}'.format(each[4])+'{0: <10}'.format(each[5])
    angle.append(angleline)
dihedral=[]#writing dihedral info in correct GROMACS itp format
for each in dihedrals:
    dihedrallist_noside.extend([each.split()])
for each in dihedrallist_noside:
    if (len(each)==8):
        dihedralline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <8}'.format(each[5])+'{0: <8}'.format(each[6])+'{0: <8}'.format(each[7])
        dihedral.append(dihedralline)
    elif (len(each)==7):
        dihedralline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <8}'.format(each[5])+'{0: <8}'.format(each[6])
        dihedral.append(dihedralline)
    elif (len(each)==11):
        dihedralline='{0: <8}'.format(each[0])+'{0: <8}'.format(each[1])+'{0: <8}'.format(each[2])+'{0: <8}'.format(each[3])+'{0: <8}'.format(each[4])+'{0: <8}'.format(each[5])+'{0: <8}'.format(each[6])+'{0: <8}'.format(each[7])+'{0: <8}'.format(each[8])+'{0: <8}'.format(each[9])+'{0: <8}'.format(each[10])
        dihedral.append(dihedralline)


with open(path_output_oeff+moleculename+'_RU_SC.xyz', 'w') as f:
    f.write(str(len(XYZ))+"\n\n")
    f.writelines("%s\n" % l for l in XYZ)
with open(path_output_oeff+moleculename+'_RU_SC.itp', 'w') as f:
    f.write(";\n; GENERATED BY GAMMPS\n; Troisi Lab @ University of Liverpool\n")
    f.write("[ atomtypes ]\n")
    f.writelines("%s\n" % l for l in atomtype)
    f.write("[ moleculetype ]\n;Name     nrexcl\n")
    f.writelines(moleculename+'_RU_SC     3\n')
    f.write("[ atoms ]\n")
    f.writelines("%s\n" % l for l in atom)
    f.write("[ bonds ]\n")
    f.writelines("%s\n" % l for l in bond)
    f.write("[ pairs ]\n")
    f.writelines("%s\n" % l for l in pairs) #only for side chains
    f.write("[ angles ]\n")
    f.writelines("%s\n" % l for l in angle)
    f.write("[ dihedrals ]\n")
    f.writelines("%s\n" % l for l in dihedral)
#    f.write("[ dihedrals ]\n")
#    f.writelines("%s\n" % l for l in improper)
with open(path_output_oeff+moleculename+'_RU_SC.top', "w") as f:
    f.write("[ defaults ]\n")
    f.write('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
    f.write('1     3     yes     0.5     0.5\n\n')
    f.write('#include "'+moleculename+'_RU_SC.itp"\n\n')
    f.write('[ system ]\n; Name\n'+moleculename+'_RU_SC\n\n')
    f.write('[ molecules ]\n'+moleculename+'_RU_SC  1\n')
end = time.time()

with open(path_output_oeff+'/RU_SC.tcl', "w") as f:
    f.write('mol new '+moleculename+'_RU_SC.xyz type {xyz}\n\nmol modstyle 0 0 CPK 1.000000 0.300000 120.000000 120.000000\n\nmol modcolor 0 0 Name\n\ndisplay resetview\n\naxes location Off\n\n')
    f.write("color Display Background white\n\ndisplay projection orthographic\n\n"+"render Tachyon "+moleculename+"_RU_SC "+VMD_RENDER+" -aasamples 12 %s -format TARGA -o %s.tga -res 1000 1000"+"\nexit\n")
subprocess.run(['vmd','-dispdev','text','-e', 'RU_SC.tcl'], cwd = path_output_oeff)
subprocess.run(['magick', moleculename+'_RU_SC.tga', moleculename+'_RU_SC.jpeg'], cwd = path_output_oeff)

with open(path_output+'.html', "a") as f:
    f.write('<p> Atom ID of C atoms to which sidechains are attached: '+str(C_connection)+'</p>\n')
    f.write('<p> Atom ID of H atoms renamed to C due to sidechain attachements are: '+str(H_renamed_C)+'</p>\n')
    f.write('<aside class="figures">\n')
    f.write('<figure>\n<img src="'+moleculename+'/'+moleculename+'_RU_SC.jpeg" style="width: 500px; height: 300px; object-fit: cover;" alt="RU_SC">\n<figcaption>repeat unit with sidechains</figcaption>\n</figure>')
    f.write('</aside>\n')

 
