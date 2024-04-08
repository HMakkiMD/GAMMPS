import collections
import numpy as np
import copy as cp
import subprocess
import os
import time
import re

'''
This is the third code to run in this step.

This code generates the chain and surrounding coordinate (.xyz) files as input for DOS and localization length calculations in Step 4
from the trajectory of MD "soup" simulation as explained in (DOI). The path to the trajectory and the type of trajectory (e.g., .xtc or .trr) 
should be given in the DOSparameters.py.

'''


from DOSparameters import *
moleculename=OLIGOMERNAME
path_to_dosindex=PATH_OUTPUT+moleculename+'/DOS/'
path_to_dosinput=PATH_MD+moleculename+INPUT_TRJ
path_to_dosoutput=path_to_dosinput+'/DOS_'+str(THRESHOLD)+'/'
ndx_file=path_to_dosindex+NDX_FILENAME+'.ndx'
#trr_file=EQ_MODEL+TRJ_TYPE
#tpr_file=EQ_MODEL+'.tpr'
frames=FRAMES
for frame in range(frames):
    '''
    if (SUBMITJOB=='yes'):
        subprocess.run(['mkdir','DOS'],cwd=path_to_dosinput)
        subprocess.run(['gmx','trjconv','-f',trr_file,'-s',tpr_file,'-n', ndx_file,'-o','DOS/pbcmol_conj_'+str(frame)+'.gro','-dump',str(TIME_INTERVAL*(frame+1)),'-pbc','whole'],cwd = path_to_dosinput)
        subprocess.run(['gmx','trjconv','-f',trr_file,'-s',tpr_file,'-n', ndx_file,'-o','DOS/nopbc_conj_'+str(frame)+'.gro','-dump',str(TIME_INTERVAL*(frame+1)),'-pbc','no'],cwd = path_to_dosinput)
    '''
    if (SUBMITJOB=='yes'):
        trr_file=EQ_MODEL+str(frame+1)+TRJ_TYPE
        tpr_file=EQ_MODEL+str(frame+1)+'.tpr'
        subprocess.run(['mkdir','DOS_'+str(THRESHOLD)],cwd=path_to_dosinput)
        subprocess.run(['gmx','trjconv','-f',trr_file,'-s',tpr_file,'-n', ndx_file,'-o','DOS_'+str(THRESHOLD)+'/pbcmol_conj_'+str(frame)+'.gro','-dump',str(TIME_INTERVAL),'-pbc','whole'],cwd = path_to_dosinput)
        subprocess.run(['gmx','trjconv','-f',trr_file,'-s',tpr_file,'-n', ndx_file,'-o','DOS_'+str(THRESHOLD)+'/nopbc_conj_'+str(frame)+'.gro','-dump',str(TIME_INTERVAL),'-pbc','no'],cwd = path_to_dosinput)


    with open(path_to_dosoutput+'pbcmol_conj_'+str(frame)+'.gro', 'r') as f:
        lines=f.readlines()[2:atoms_per_chain+3]
        atomtype_chain=[]
        for each in lines[:-1]:
            atomtype_chain.append(each[14])
    with open(path_to_dosoutput+'pbcmol_conj_'+str(frame)+'.gro', 'r') as f:
        lines=f.readlines()[2:atoms_per_0chargegroup+3]
        atomtype_0chargegroup=[]
        for each in lines[:-1]:
            atomtype_0chargegroup.append(each[14])
    with open(path_to_dosoutput+'nopbc_conj_'+str(frame)+'.gro', 'r') as f:
        lines=f.readlines()[2:]
        nopbc=[]
        for each in lines[:-1]:
            coordinate=each.split()
            nopbc.append([float(coordinate[-3]), float(coordinate[-2]), float(coordinate[-1])])  
        box=lines[-1].split()
        box_dim=[float(box[0]),float(box[1]),float(box[2])]
    with open(path_to_dosoutput+'pbcmol_conj_'+str(frame)+'.gro', 'r') as f:
        lines=f.readlines()[2:atoms_per_chain+3]
        pbcmol=[]
        for each in lines[:-1]:
            coordinate=each.split()
            pbcmol.append([float(coordinate[-3]), float(coordinate[-2]), float(coordinate[-1])])
    chains = [pbcmol[x:x+atoms_per_chain] for x in range(0, len(pbcmol), atoms_per_chain)]
    #chains_nopbc = [nopbc[x:x+atoms_per_chain] for x in range(0, len(nopbc), atoms_per_chain)]
    for index, chain in enumerate(chains):
        images=[]
        for atom in chain:
            images.append(tuple([atom[0]//box_dim[0],atom[1]//box_dim[1],atom[2]//box_dim[2]]))
            images.append(tuple([(atom[0]+threshold_mon)//box_dim[0],(atom[1]+threshold_mon)//box_dim[1],(atom[2]+threshold_mon)//box_dim[2]]))
            images.append(tuple([(atom[0]-threshold_mon)//box_dim[0],(atom[1]-threshold_mon)//box_dim[1],(atom[2]-threshold_mon)//box_dim[2]]))
        images = list(dict.fromkeys(images))
        new_coordinates=[]
        for image in images:
            for atomnr, atom in enumerate(nopbc[atoms_per_chain:]):
                #if(atomnr not in range(index*atoms_per_chain,(index+1)*atoms_per_chain,1)):#excluding the chain of interest
                if (atomnr%(atoms_per_chain_surr) < (atoms_per_chain_surr-atoms_of_endgroups)):
                    new_x=image[0]*box_dim[0]+atom[0]
                    new_y=image[1]*box_dim[1]+atom[1]
                    new_z=image[2]*box_dim[2]+atom[2]
                    new_coordinates.append([new_x,new_y,new_z])
        _0chargegroups = [new_coordinates[x:x+atoms_per_0chargegroup] for x in range(0, len(new_coordinates), atoms_per_0chargegroup)]
        monomers = [chain[x:x+atoms_per_0chargegroup] for x in range(0, len(chain), atoms_per_0chargegroup)]
        monomers.pop()
        center_geom=[]
        for each in monomers:
            center_x = [ sum(row[i] for row in each) for i in range(len(each[0])) ][0] / len(each)
            center_y = [ sum(row[i] for row in each) for i in range(len(each[1])) ][1] / len(each)
            center_z = [ sum(row[i] for row in each) for i in range(len(each[2])) ][2] / len(each)
            center_geom.append([center_x,center_y,center_z])
          
        surronding=[]
        surronding_index=[]
        rightpos=[]
        flag=False
        for each in chain:
        #for each in center_geom:
            for index_0charge, _0charge in enumerate(_0chargegroups):
                for atom in _0charge:
                    if (np.linalg.norm(np.array(atom)-np.array(each))<threshold_mon):
                        flag=True
                        break
                if (flag==True):
                    if (index_0charge not in surronding_index):
                        surronding.append(_0charge)
                        rightpos.append([atom[0],atom[1],atom[2]])
                    surronding_index.append(index_0charge)
                flag=False
        
        surronding_mod=cp.deepcopy(surronding)
        
        for index_surr, each in enumerate(surronding_mod):
            for item in each:
                if (item[0]-rightpos[index_surr][0]>largest_dist):
                    item[0]-=box_dim[0]
                if (item[1]-rightpos[index_surr][1]>largest_dist):
                    item[1]-=box_dim[1]
                if (item[2]-rightpos[index_surr][2]>largest_dist):
                    item[2]-=box_dim[2]
                if (rightpos[index_surr][0]-item[0]>largest_dist):
                    item[0]+=box_dim[0]
                if (rightpos[index_surr][1]-item[1]>largest_dist):
                    item[1]+=box_dim[1]
                if (rightpos[index_surr][2]-item[2]>largest_dist):
                    item[2]+=box_dim[2]        
        to_be_del=[]
        for i in range(len(surronding_mod[:-1])):
            for j in range(i+1,len(surronding_mod)):
                if(np.sqrt(((((np.array(surronding_mod[i]) - np.array(surronding_mod[j]))** 2))*3).mean())<0.0001):
                    to_be_del.append(j)
        surronding_mod_removedrep=[]
        for i in range(len(surronding_mod)):
            if (i not in to_be_del):
                surronding_mod_removedrep.append(surronding_mod[i])

            xyz_chain=[]
            for i, each in enumerate(chain):
                xyz_chain.append('{0: <12}'.format(atomtype_chain[i])+'{0: <12}'.format(round(each[0]*10,5))+'{0: <12}'.format(round(each[1]*10,5))+'{0: <12}'.format(round(each[2]*10,5)))
            with open(path_to_dosoutput+str(frame+1)+'_chain.xyz', 'w') as f:
                f.write(str(len(xyz_chain))+"\n\n")
                f.writelines("%s\n" % l for l in xyz_chain)
            xyz_surronding=[]
            for each in surronding_mod_removedrep:
                for index_item, item in enumerate(each):
                    xyz_surronding.append('{0: <12}'.format(atomtype_0chargegroup[index_item])+'{0: <12}'.format(round(item[0]*10,5))+'{0: <12}'.format(round(item[1]*10,5))+'{0: <12}'.format(round(item[2]*10,5)))
            with open(path_to_dosoutput+str(frame+1)+'_sur.xyz', 'w') as f:
                f.write(str(len(xyz_surronding))+"\n\n")
                f.writelines("%s\n" % l for l in xyz_surronding)
    


