'''
This code generates Gaussian 16 input files from the chain and surrounding coordinate (.xyz)
files located in the 'input_files' directory and the monomer charge file (ending in 'charges.txt') 
located in the current working directory. It submits a single-point calculation for each generated 
input file utilizing 4 CPU cores per calculation (max. 40 CPU cores total in parallel).

The bulk DOS and LL is calculated by averaging over the DOS and LL computed for each 
sample and written to the files '{polymer_name}_DOS_{broadening_1string}.txt' and '{polymer_name}_LL_{broadening_1string}.txt', 
where the variables enclosed in {} should be specified in the 'input_variables.inp' file.

None of the input ('.com'), checkpoint ('.chk'), output ('.log'), or ('.npz') files are retained in the 
the current version of the code.
'''

import numpy as np
import os
import re
import sys
import subprocess # To use bash and external softwares such as g16
import datetime # To get current date and time
import time # To get current date and time
from collections import OrderedDict

#Multiprocess
from multiprocessing import Pool

import cclib

import scipy.sparse as sp
#import matplotlib.pyplot as plt

def checkfile(filename):

    if not os.path.isfile(filename):
        print(" File %s not found!" % filename)
        sys.exit()

def fill_dict(filename):
    '''Fills a dictionary of options.
       Courtesy of Daniele Padula'''

    opts = OrderedDict()
    try:
        checkfile(filename) # Check if the file exists (subroutine)
        with open(filename) as f:
            for line in f:

                # Ignore comments and empty lines
                if line.startswith('#'): #Comments
                    continue

                if not line.strip(): #Empty lines
                    continue

                # AL: Strippo i caratteri inutili
                line = line.replace('=', '')
                line = line.replace('[', '')
                line = line.replace(']', '')
                line = line.replace(',', ' ') # commas with spaces
                line = line.replace('"', ' ') # quotes with spaces
                line = line.replace("'", " ") # single quotes with spaces

                key_dict = line.split()[0] # Keys of dictionary
                data = line.split()[1:] # Values of dictionary

                # Try to understand whether data should be stored as int, float
                # or string
                try:
                    data = np.array(list(map(int, data)))

                # This can occur with both floats and strings
                except ValueError:
                    try:
                        data = np.array(list(map(float, data)))

                    except ValueError:
                        data = np.array(list(map(str, data)))

                if len(data) == 1: # If 1 value only, then extract from list
                    opts[key_dict] = data[0]
                else:
                    opts[key_dict] = data

    except TypeError:
        pass

    return opts

#################################################################################################################################

print('{:>60}'.format(" !-----------------------------------------------------!"))
print('{:>60}'.format(" !                                                     !"))
print('{:>60}'.format(" !  Script to perform electronic structure calculation !"))
print('{:>60}'.format(" !      for an ensemble of SCP chain conformations     !"))
print('{:>60}'.format(" !                                                     !"))
print('{:>60}'.format(" !                Version 1.1: 11/01/24                !"))
print('{:>60}'.format(" !                                                     !"))
print('{:>60}'.format(" !                Written by Colm Burke                !"))
print('{:>60}'.format(" !                                                     !"))
print('{:>60}'.format(" !-----------------------------------------------------!"))

now = datetime.datetime.now()
start = time.time()# To evaluate the elapsed time
print("!-----------------------------------------------------!")
print("    Computations started at:")
print (str(now))
print("!-----------------------------------------------------!")

if len(sys.argv) == 1: 
        print('ERROR: You need to indicate the input file name in the command line')
        print()
        print("   Use of the script: python QC_calculation.py input_filename    ")
        sys.exit(-1)
filename = sys.argv[1]
opts = fill_dict(filename) # Read input file 

# Set variables based on those specified in input file
if "polymer_name" in opts:
    polymer_name = opts["polymer_name"] # name of polymer
if "chain_length" in opts:
    chainlength = opts["chain_length"] # number of monomers in each sample chain
if "nsample_start" in opts:
    nsample_start = opts["nsample_start"] # start index of samples to be used
if "nsample_end" in opts:
    nsample_end = opts["nsample_end"] # end index of samples to be used
if "broadening_1" in opts:
    broadening_1 = opts["broadening_1"] # desired broadening value
if "broadening_2" in opts:
    broadening_2 = opts["broadening_2"] # desired broadening value
if "broadening_3" in opts:
    broadening_3 = opts["broadening_3"] # desired broadening value
if "broadening_4" in opts:
    broadening_4 = opts["broadening_4"] # desired broadening value    
if "E_start" in opts:
    E_start = opts["E_start"] # desired initial energy value
if "E_end" in opts:
    E_end = opts["E_end"] # desired final energy value
    
nsamples = (nsample_end - nsample_start) + 1 # V1.1 i.e., if nsample_end is 30
                                             # and nsample_start is 20 nsamples will be 11

directory = (os.getcwd() + '/input_files')

cwd = os.getcwd() # Get path of charge file 
for file in os.listdir(cwd):
    fullpath = os.path.join(cwd,file)
    if file.endswith('charges.txt'):
        charge_file = fullpath
    else:
        continue

with open(charge_file, 'r') as f: # Get force-field charges on polymer repeat unit (for MM part)
    lines = f.readlines()
    lines.pop(0)
    monomer_charges = [float(line.split()[1]) for line in lines]

# Below loop creates temporary files containing coordinates of surroundings and charges

for filename in os.listdir(directory):
    index = filename.find(".")
    prefix = filename[:index]
    fullpath = os.path.join(directory, filename)
    if filename.endswith('sur.xyz'):
        with open(fullpath, 'r') as f:
            first_line = f.readline().strip('\n')
            no_atoms = float(first_line) # get number of atoms in MM part
            monomer_number = no_atoms / len(monomer_charges) # get number of monomers
            print(monomer_number)
            monomer_number = int(monomer_number)
            print(monomer_number)                           

            charges_list = monomer_number * monomer_charges # extends monomer charges list by factor of the number of monomers

        coords = []
        with open(fullpath, 'r') as f:
            lines = f.readlines()[2:] # skip the first 2 lines
            for line in lines:
                  data = line.split()
                  x, y, z = float(data[1]), float(data[2]), float(data[3])
                  coords.append([x, y, z])

        with open('point_charges_{}.txt'.format(prefix), 'w') as f:
            lines_in_file = [[*sublist, charges_list[i]] for i,sublist in enumerate(coords)] # concatenate coordinates & charge
            f.writelines("   ".join(str(i) for i in sublist) + '\n' for sublist in lines_in_file)

# Below loop creates Gaussian input file containing just coordinates of central chain
# V1.1 - now modified to reformat coords as floats to prevent g16 read fail

for filename in os.listdir(directory):
    if filename.endswith('chain_H.xyz'):
        file_name = os.path.splitext(filename)[0]
        fullpath = os.path.join(directory, filename)
        number = re.findall(r'\d+', file_name)
        number = str(number[0])
        with open(fullpath, 'r') as f1, open(f'{polymer_name}_{number}.com', "w", newline='\n') as output:
            output.write('%nproc=4' + '\n') # specify number of procs per g16 calc (hardcoded now but could be a variable)
            output.write('%mem=20GB' + '\n') # specify total memory for g16 calc (again hardcoded)
            output.write(f'%chk={polymer_name}_{number}.chk' + '\n') # specify chk file
            output.write('#p b3lyp/3-21g* charge nosymm IOp(3/33=1,5/33=1,6/8=2,6/9=2,6/10=2,6/11=2) pop=full' + '\n') # this specifies basis set and full population analysis 
            output.write('\n')
            output.write(f'Gaussian single point energy calculation on {file_name}.xyz' + '\n')
            output.write('\n')
            output.write('0,1' + '\n')
            i = 0
            for line in f1:
               if i < 2:
                   i += 1
                   continue               
               data = line.split()
               x, y, z = float(data[1]), float(data[2]), float(data[3])
               output.write(str(data[0])+"   "+str(x)+"   "+str(y)+"   "+str(z)+'\n')
            output.write('\n')

# Below loop appends point charges to the corresponding Gaussian input file (by matching the chain number) then removes point charge files

for filename in os.listdir(os.curdir):
    if filename.startswith('point_charges'):
        number = re.findall(r'\d+', filename)
        number = str(number[0])
        with open(filename, 'r') as f2, open(f'{polymer_name}_{number}.com', "a", newline='\n') as output:
            for line in f2:
                output.write(line)
            output.write('\n')
        os.remove(filename)

now = datetime.datetime.now()
print("!-----------------------------------------------------!")
print("    g16 input generated at:")
print (str(now))
print("!-----------------------------------------------------!")
########################################################################################
################################ LAUNCHING G16 JOBS ###################################

def run_g16_parallel_list(filename):
    '''
    Used to launch Gaussian calculations in parallel
    'filename' = name of Gaussian input file 

    '''
    # Check the files are not launched more than 1 times
    print(filename, flush = True)

    return filename, subprocess.call(["g16", filename], cwd="./")

list_calcg16 = []

for i in range(nsample_start,(nsample_end+1)): # V1.1 only the samples specified will be launched
    list_calcg16.append(f'{polymer_name}_{i}.com')
        
with Pool(10) as p: # 10 concurrent worker processes (utilising 4 cores per calculation)
    results = []
    r = p.map_async(run_g16_parallel_list, list_calcg16,
                    callback=results.append)
    r.get()
    results = results[0]
    
now = datetime.datetime.now()
print("!-----------------------------------------------------!")
print("    g16 calculations finished at:")
print (str(now))
print("!-----------------------------------------------------!")
#########################################################################################
######################### EXTRACTING NECESSARY DATA FROM LOG ########################## 

def save_to_npz(file): # V1.1: Now a try-except block to catch when g16 job has failed
    
    try:
        prefix,extension = os.path.splitext(file)
        
        data = cclib.io.ccread(file)

        overlap = data.aooverlaps
        mocoeffs = data.mocoeffs
        moenergies = data.moenergies
        coords = data.atomcoords
        masses = data.atommasses
        nbasis = data.nbasis
        natoms = data.natom
        mocoeffs = mocoeffs[0]

        np.savez(f'{prefix}.npz', coords=coords,masses=masses,overlap=overlap,mocoeffs=mocoeffs,moenergies=moenergies,nbasis=nbasis,natoms=natoms)
        return None
    
    except Exception as e:
        print(f'Error processing file {file}: {e}')
        return file
                
log_files = []
for file in os.listdir(os.curdir):
    if file.endswith('.log'):
        log_files.append(file)

if nsamples >= 40:
    ntasks_savelog = 40
else: 
    ntasks_savelog = nsamples

with Pool(ntasks_savelog) as p: # 40 concurrent worker processes (1 core per process)
    results = []
    r = p.map_async(save_to_npz, log_files,
                    callback=results.append)
    r.get()
    
num_failures = sum(result is not None for result in results)
nsamples = nsamples - num_failures # V1.1: updates number of samples to account for failures

for filename in os.listdir(os.curdir): #V1.1: remove uneccessary chk and log files
    if filename.endswith('.chk'):
        os.remove(filename)
    elif filename.endswith('.log'):
        os.remove(filename)
    elif filename.endswith('.com'):
        os.remove(filename)
        
now = datetime.datetime.now()
print("!-----------------------------------------------------!")
print("    Necessary data extracted from log files at:")
print (str(now))
print("!-----------------------------------------------------!")
#########################################################################################
############## CALCULATING + PLOTTING DOS & LOCALISATION LENGTH #######################

def calculate_centroid(points): # From input of points (x,y,z) returns centroid of points
    
    if not points:
        return None  # Return None for an empty list of points

    # Initialize the sums of x, y, and z coordinates to zero
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # Iterate through each point and add its coordinates to the sums
    for point in points:
        x, y, z = point
        sum_x += x
        sum_y += y
        sum_z += z

    # Calculate the average coordinates to get the centroid
    num_points = len(points)
    centroid_x = sum_x / num_points
    centroid_y = sum_y / num_points
    centroid_z = sum_z / num_points

    return [centroid_x, centroid_y, centroid_z]

def center_of_mass(points, masses): # From input of points (x,y,z) and corresponding masses, returns center of mass of points 

    if len(points) != len(masses):
        raise ValueError("The number of points and masses must be the same.")

    total_mass = sum(masses)
    if total_mass == 0:
        raise ValueError("Total mass must be greater than zero.")

    weighted_sum = [0, 0, 0]  # Initialize a list to store the weighted sum of coordinates

    for i in range(len(points)):
        for j in range(3):  # Iterate over x, y, and z coordinates
            weighted_sum[j] += points[i][j] * masses[i]

    center_x = weighted_sum[0] / total_mass
    center_y = weighted_sum[1] / total_mass
    center_z = weighted_sum[2] / total_mass

    return (center_x, center_y, center_z)

def calc_LL(W, E, sigma, Emin, Emax, dE, LL_m):
       
    NP = int((Emax-Emin)/dE)
    X = np.linspace(Emin, Emax, NP)
    C = -1/(2*sigma**2)
    values_allE = np.zeros(NP)
    values_MO = np.zeros(len(E))
    values_MO_IPR = np.zeros(len(E))
    sum_MO = np.sum(values_MO)
    constant = 0.39894228 / sigma / 10
    
    for i in range(len(X)):
        values_MO[:] = W[:] * np.exp(C*(X[i]-E[:])**2)
        values_MO *= constant
        values_MO_IPR[:] = values_MO[:] * LL_m[:]
        sum_MO = np.sum(values_MO)
        values_allE[i] = np.sum(values_MO_IPR) / sum_MO if sum_MO != 0 else 0
        
    return X, values_allE

def calc_LL_m(mocoeffs, nbasis, Pkm_ranges, overlap, coords, monomer_ranges, masses, chainlength):
    
    # This function calculates the localisation length for each molecular orbital in the chain 
    
    Pkms = np.zeros((chainlength, nbasis))
    coms = []

    overlap_sparse = sp.csr_matrix(overlap)
    
    # Below calculates weight of each molecular orbital on each monomer
    
    for k, (start, stop, special_j_values) in enumerate(Pkm_ranges): 
        j_range = np.arange(start, stop)
        if special_j_values:
            j_range = np.concatenate([j_range, special_j_values])
        j_indices = np.hstack([j_range])
        j_values = mocoeffs[:, j_indices] * overlap_sparse[j_indices, :] * mocoeffs
        Pkms[k] = np.sum(j_values, axis=1)

    # Below calculates center of mass of each monomer
 
    for k, (start, stop, special_j_values) in enumerate(monomer_ranges):
        monomer_coords = []
        monomer_masses = []
        j_range = np.arange(start, stop)
        if special_j_values:
            j_range = np.concatenate([j_range, special_j_values])
        for i in j_range:
            monomer_coords.append(coords[i])
            monomer_masses.append(masses[i])

        coms.append(center_of_mass(monomer_coords, monomer_masses))

    MO_centroids = []
    for m in range(nbasis):
        R_m = np.zeros(3)  # Initialize R^{(m)} as a 3D array to store x, y, z coordinates

    # Iterate over each monomer (k)
        for k in range(chainlength):
        # Get the weight of molecular orbital m on monomer k
            P_km = Pkms[k, m]

        # Get the centroid coordinates of monomer k
            com_k = np.array(coms[k])

        # Calculate the contribution of monomer k to R^{(m)} and add it
            R_m += com_k * P_km

    # Append the calculated R^{(m)} to the list
        MO_centroids.append(R_m)

    LL_m = np.zeros(nbasis)
    for m in range(nbasis):

        loclength_m = 0

        for k in range(chainlength):
            P_km = Pkms[k,m]
            if P_km < 0:
                P_km = abs(P_km)
            R_m = np.array(MO_centroids[m])
            r_k = np.array(coms[k])

            LL = abs(np.linalg.norm(r_k - R_m))
            LL = LL ** 2
            LL = LL * P_km
            loclength_m += LL

        loclength_m = 2 * np.sqrt(loclength_m)
        LL_m[m] = loclength_m

    return Pkms, LL_m

def broaden(W, E, sigma, Emin, Emax, dE, chainlength, nchains):
    
    NP = int((Emax-Emin)/dE) # Number of points at which E-dep DOS is evaluated
    X = np.linspace(Emin, Emax, NP) # X axis of DOS
    C = -1/(2*sigma**2) # constant term 
    
    f = np.dot(W, np.exp(C * (np.outer(X, np.ones(len(E))) - np.outer(np.ones(NP), E))**2).T) # calculate DOS 
    
    # Normalise the DOS by chain length  
    f = 0.39894228 * f / sigma / chainlength
    f = f / nchains 
    
    return (X,f)

def process_directories(broaden_func, Pkm_ranges,Monomer_ranges,chainlength,nsamples,E_start,E_end,broadening):
    
    LL_lists = [] # initiate empty list which contains the energy-dependent LL for each chain
    f_lists = [] # initiate empty list which contains the enrgy-dependent DOS for each chain

    cwd = os.getcwd()

    for i, filename in enumerate(os.listdir(cwd)): # Looping over all '.npz' files in working directory
        if filename.endswith('.npz'):
            file_path = os.path.join(cwd, filename)

            test = np.load(file_path) # Read necessary data from file
            overlap = test['overlap']
            mocoeffs = test['mocoeffs']
            MO_energies = test['moenergies']
            masses = test['masses']
            coords = test['coords']
            coords = coords[0]
            nbasis = len(mocoeffs)

            MO_energies = [item for sublist in MO_energies for item in sublist]

            Pkms, LL_m = calc_LL_m(mocoeffs, nbasis, Pkm_ranges, overlap, coords, Monomer_ranges, masses, chainlength) # Calculate LL for each MO 

            weights = [1 for i in range(len(MO_energies))] # REDUNDANT IN THIS CODE

            X, values_allE = calc_LL(weights, MO_energies, broadening, E_start, E_end, 0.01, LL_m) # Calculate energy dependent LL for current chain

            LL_lists.append(values_allE) # Append energy-dependent LL to list
            
            (X, f) = broaden_func(weights, MO_energies, broadening, E_start, E_end, 0.01, chainlength, nsamples) # Calculate energy dependent DOS for current chain
            f_lists.append(f) # Append energy-dependent DOS to list
            
    bulk_DOS = np.sum(f_lists,axis=0) # Average to get bulk DOS

    bulk_LL = np.sum(LL_lists, axis=0) # Average to get bulk LL
    bulk_LL = [x / nsamples for x in bulk_LL] # Divide by nchains to get averaged LL (done for DOS within broaden function)
    
    return X, bulk_LL, bulk_DOS

def LL_DOS_1():
    return process_directories(broaden, coeff_ranges, monomer_ranges, chainlength, nsamples,E_start,E_end,broadening_1)

def LL_DOS_2():
    return process_directories(broaden, coeff_ranges, monomer_ranges, chainlength, nsamples,E_start,E_end,broadening_2)

def LL_DOS_3():
    return process_directories(broaden, coeff_ranges, monomer_ranges, chainlength, nsamples,E_start,E_end,broadening_3)

def LL_DOS_4():
    return process_directories(broaden, coeff_ranges, monomer_ranges, chainlength, nsamples,E_start,E_end,broadening_4)

cwd = os.getcwd()
for file in os.listdir(cwd): # Reads one .npz file to get nbasis and natoms values
    if file.endswith('.npz'):
        file_path = os.path.join(cwd, file)
        test = np.load(file_path)
        mocoeffs = test['mocoeffs']
        natoms = test['natoms']
        nbasis = len(mocoeffs)
        break

def calculate_coeff_ranges(nbasis, num_monomers): #V1.1: coefficient ranges now generated based on specified chain length
                                                  
    coeff_ranges = []

    for i in range(num_monomers):
        start = i * ((nbasis - 4) // num_monomers)
        end = (i + 1) * ((nbasis - 4) // num_monomers) if i != num_monomers - 1 else nbasis - 4

        if i == 0:
            start_range = (0, end, [(nbasis - 4), (nbasis - 3)])
        elif i == num_monomers - 1:
            start_range = (start, end, [(nbasis - 2), (nbasis - 1)])
        else:
            start_range = (start, end, [])

        coeff_ranges.append(start_range)

    return coeff_ranges

def calculate_monomer_ranges(natoms, num_monomers): #V1.1: monomer ranges now generated based on specified chain length
    monomer_ranges = []

    for i in range(num_monomers):
        start = i * ((natoms - 2) // num_monomers)
        end = (i + 1) * ((natoms - 2) // num_monomers) if i != num_monomers - 1 else natoms - 2

        if i == 0:
            start_range = (0, end, [(natoms - 2)])
        elif i == num_monomers - 1:
            start_range = (start, end, [(natoms - 1)])
        else:
            start_range = (start, end, [])

        monomer_ranges.append(start_range)

    return monomer_ranges

coeff_ranges = (calculate_coeff_ranges(nbasis, chainlength))
monomer_ranges = (calculate_monomer_ranges(natoms,chainlength))

if 'broadening_1' in locals() and locals()['broadening_1'] is not None:
    X, bulk_LL_1, bulk_DOS_1 = LL_DOS_1() # call function to calculate DOS and LL
    
    broadening_1string = str(broadening_1).replace('.', '')
    
    with open(f'{polymer_name}_DOS_{broadening_1string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_DOS_1):
            file.write(f"{item1}\t{item2}\n")
            
    with open(f'{polymer_name}_LL_{broadening_1string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_LL_1):
            file.write(f"{item1}\t{item2}\n")

if 'broadening_2' in locals() and locals()['broadening_2'] is not None:
    X, bulk_LL_2, bulk_DOS_2 = LL_DOS_2() # call function to calculate DOS and LL
    
    broadening_2string = str(broadening_2).replace('.', '')
    
    with open(f'{polymer_name}_DOS_{broadening_2string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_DOS_2):
            file.write(f"{item1}\t{item2}\n")
            
    with open(f'{polymer_name}_LL_{broadening_2string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_LL_2):
            file.write(f"{item1}\t{item2}\n")
    
if 'broadening_3' in locals() and locals()['broadening_3'] is not None:
    X, bulk_LL_3, bulk_DOS_3 = LL_DOS_3() # call function to calculate DOS and LL
    
    broadening_3string = str(broadening_3).replace('.', '')
    
    with open(f'{polymer_name}_DOS_{broadening_3string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_DOS_3):
            file.write(f"{item1}\t{item2}\n")
            
    with open(f'{polymer_name}_LL_{broadening_3string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_LL_3):
            file.write(f"{item1}\t{item2}\n")
    
if 'broadening_4' in locals() and locals()['broadening_4'] is not None:
    X, bulk_LL_4, bulk_DOS_4 = LL_DOS_4() # call function to calculate DOS and LL
    
    broadening_4string = str(broadening_4).replace('.', '')

    with open(f'{polymer_name}_DOS_{broadening_4string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_DOS_4):
            file.write(f"{item1}\t{item2}\n")
            
    with open(f'{polymer_name}_LL_{broadening_4string}.txt', 'w') as file:
        # Iterate through both lists simultaneously
        for item1, item2 in zip(X, bulk_LL_4):
            file.write(f"{item1}\t{item2}\n")


for filename in os.listdir(os.curdir):
    if filename.endswith('.npz'):
        os.remove(filename)
        
####################### END OF LOOP #########################
        
now = datetime.datetime.now()
print("!-----------------------------------------------------!")
print("    Computations finished at:")
print (str(now))
print("!-----------------------------------------------------!")
end = time.time()# To evaluate the elapsed time

print("!-----------------------------------------------------!")
print("    Elapsed time (s):")
print(str(end - start))
print("!-----------------------------------------------------!")
