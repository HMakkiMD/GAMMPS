
PATH_OUTPUT='/users/hmakki/sharedscratch/GAMMPS/Polymers/'
PATH_MD='/users/hmakki/sharedscratch/MD/'

'''bulk
EQ_MODEL='eq_1'
TRJ_TYPE='.xtc'
OLIGOMERNAME='IDTBO'
NDX_FILENAME='IDTBOcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=492 # total number of atoms in each chain
atoms_per_chain_surr=492 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=49 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup mono
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='IDTBO'
NDX_FILENAME='IDTBOcharge_10_1'
INPUT_TRJ='/soup_mono/900/DOS/300'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=492 # total number of atoms in each chain
atoms_per_chain_surr=51 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=49 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup di
EQ_MODEL='900'
TRJ_TYPE='.trr'
OLIGOMERNAME='IDTBO'
NDX_FILENAME='IDTBOcharge_10_2'
INPUT_TRJ='/soup_di/900'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=492 # total number of atoms in each chain
atoms_per_chain_surr=100 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=49 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''
EQ_MODEL='eq_1'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPTT'
NDX_FILENAME='DPPTTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=382 # total number of atoms in each chain
atoms_per_chain_surr=382 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=38 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup mono
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPTT'
NDX_FILENAME='DPPTTcharge_10_1'
INPUT_TRJ='/soup_mono/900/DOS/300'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=382 # total number of atoms in each chain
atoms_per_chain_surr=40 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=38 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup di
EQ_MODEL='900'
TRJ_TYPE='.trr'
OLIGOMERNAME='DPPTT'
NDX_FILENAME='DPPTTcharge_10_2'
INPUT_TRJ='/soup_di/900'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=382 # total number of atoms in each chain
atoms_per_chain_surr=78 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=38 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''
EQ_MODEL='eq_1'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPBTZ'
NDX_FILENAME='DPPBTZcharge_10'
INPUT_TRJ='/bulk/'
FRAMES=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=472 # total number of atoms in each chain
atoms_per_chain_surr=472 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=47 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup mono
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPBTZ'
NDX_FILENAME='DPPBTZcharge_10_1'
INPUT_TRJ='/soup_mono/900/DOS/300'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=472 # total number of atoms in each chain
atoms_per_chain_surr=49 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=47 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup di
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPBTZ'
NDX_FILENAME='DPPBTZcharge_10_2'
INPUT_TRJ='/soup_di/900/DOS/300'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=472 # total number of atoms in each chain
atoms_per_chain_surr=96 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=47 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''
EQ_MODEL='eq_1'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPDTT'
NDX_FILENAME='DPPDTTcharge_10'
INPUT_TRJ='/bulk/'
FRAMES=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=422 # total number of atoms in each chain
atoms_per_chain_surr=422 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=42 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup mono
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPDTT'
NDX_FILENAME='DPPDTTcharge_10_1'
INPUT_TRJ='/soup_mono/900/DOS/300'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=422 # total number of atoms in each chain
atoms_per_chain_surr=44 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=42 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup di
EQ_MODEL='900'
TRJ_TYPE='.trr'
OLIGOMERNAME='DPPDTT'
NDX_FILENAME='DPPDTTcharge_10_2'
INPUT_TRJ='/soup_di/900'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=422 # total number of atoms in each chain
atoms_per_chain_surr=86 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=42 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''bulk
EQ_MODEL='eq_1'
TRJ_TYPE='.xtc'
OLIGOMERNAME='IDTBT'
NDX_FILENAME='IDTBTcharge_10'
INPUT_TRJ='/bulk/'
FRAMES=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=492 # total number of atoms in each chain
atoms_per_chain_surr=492 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=49 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup mono
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='IDTBT'
NDX_FILENAME='IDTBTcharge_10_1'
INPUT_TRJ='/soup_mono/900/DOS/300'
FRAMES=250
TIME_INTERVAL=1000 #ps between each frame
atoms_per_chain=492 # total number of atoms in each chain
atoms_per_chain_surr=51 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=49 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup di
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='IDTBT'
NDX_FILENAME='IDTBTcharge_10_2'
INPUT_TRJ='/soup_di/900/300'
FRAMES=250
TIME_INTERVAL=1000 #ps between each frame
atoms_per_chain=492 # total number of atoms in each chain
atoms_per_chain_surr=100 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=49 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='CPDTBFF'
NDX_FILENAME='CPDTBFFcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=312 # total number of atoms in each chain
atoms_per_chain_surr=312 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=31 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=2.5 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup mono
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='CPDTBFF'
NDX_FILENAME='CPDTBFFcharge_10_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=4000 #ps between each frame
atoms_per_chain=312 # total number of atoms in each chain
atoms_per_chain_surr=33 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=31 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=2.5 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup di
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='CPDTBFF'
NDX_FILENAME='CPDTBFFcharge_10_2'
INPUT_TRJ='/soup_di/900/300'
FRAMES=250
TIME_INTERVAL=1000 #ps between each frame
atoms_per_chain=312 # total number of atoms in each chain
atoms_per_chain_surr=64 # total number of atoms in each chain of the surronding molecules 
atoms_per_0chargegroup=31 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=2.5 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='TIFBT'
NDX_FILENAME='TIFBTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=612 # total number of atoms in each chain
atoms_per_chain_surr=612 # total number of atoms in each chain of the surronding molecules
atoms_per_0chargegroup=61 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side ch$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup mono
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='TIFBT'
NDX_FILENAME='TIFBTcharge_10_1'
INPUT_TRJ='/soup_mono/900/DOS/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
atoms_per_chain=612 # total number of atoms in each chain
atoms_per_chain_surr=63 # total number of atoms in each chain of the surronding molecules
atoms_per_0chargegroup=61 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side ch$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''soup di
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='TIFBT'
NDX_FILENAME='TIFBTcharge_10_2'
INPUT_TRJ='/soup_di/900/DOS/300'
FRAMES=250
TIME_INTERVAL=1000 #ps between each frame
atoms_per_chain=612 # total number of atoms in each chain
atoms_per_chain_surr=124 # total number of atoms in each chain of the surronding molecules
atoms_per_0chargegroup=61 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side ch$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = 2.0 # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

''' 1chain
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPBTZ'
NDX_FILENAME='DPPBTZcharge_10'
INPUT_TRJ='/bulk/'
FRAMES=5
THRESHOLD=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=472 # total number of atoms in each chain
atoms_per_chain_surr=472 # total number of atoms in each chain of the surronding molecules
atoms_per_0chargegroup=47 # total number of atoms in each monomer to be considered for surronding charges
                    # in case of IDTBT, I consider conjogated part plus the first methyl group of the side chains
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''bulk 1chain
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='IDTBT'
NDX_FILENAME='IDTBTcharge_10'
INPUT_TRJ='/bulk/'
FRAMES=5
THRESHOLD=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=492 # total number of atoms in each chain
atoms_per_chain_surr=492 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=49 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='TIFBT'
NDX_FILENAME='TIFBTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=5
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=612 # total number of atoms in each chain
atoms_per_chain_surr=612 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=61 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup newtifbt
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='TIFBT'
NDX_FILENAME='TIFBTcharge_10_1'
INPUT_TRJ='/soup_mono/soup_400/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=3
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=612 # total number of atoms in each chain
atoms_per_chain_surr=63 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=61 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='AIGBT'
NDX_FILENAME='AIGBTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=432 # total number of atoms in each chain
atoms_per_chain_surr=432 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=43 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup 
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='AIGBT'
NDX_FILENAME='AIGBTcharge_10_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=432 # total number of atoms in each chain
atoms_per_chain_surr=45 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=43 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''


'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='CPDTBT'
NDX_FILENAME='CPDTBTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=322 # total number of atoms in each chain
atoms_per_chain_surr=322 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=32 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup 
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='CPDTBT'
NDX_FILENAME='CPDTBTcharge_10_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=322 # total number of atoms in each chain
atoms_per_chain_surr=34 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=32 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='F8BT'
NDX_FILENAME='F8BTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=382 # total number of atoms in each chain
atoms_per_chain_surr=382 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=38 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup 
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='F8BT'
NDX_FILENAME='F8BTcharge_10_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=382 # total number of atoms in each chain
atoms_per_chain_surr=40 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=38 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='P3HT'
NDX_FILENAME='P3HTcharge_20'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=202 # total number of atoms in each chain
atoms_per_chain_surr=202 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=10 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=20 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup *
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='P3HT'
NDX_FILENAME='P3HTcharge_20_1'
INPUT_TRJ='/soup_mono/600/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=202 # total number of atoms in each chain
atoms_per_chain_surr=12 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=10 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=20 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''


'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='PBTTT'
NDX_FILENAME='PBTTTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=302 # total number of atoms in each chain
atoms_per_chain_surr=302 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=30 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup 
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='PBTTT'
NDX_FILENAME='PBTTTcharge_10_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=302 # total number of atoms in each chain
atoms_per_chain_surr=32 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=30 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''


'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='PFO'
NDX_FILENAME='PFOcharge_20'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=542 # total number of atoms in each chain
atoms_per_chain_surr=542 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=27 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=20 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='PFO'
NDX_FILENAME='PFOcharge_20_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=542 # total number of atoms in each chain
atoms_per_chain_surr=29 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=27 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=20 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''


'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='PNBT'
NDX_FILENAME='PNBTcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=552 # total number of atoms in each chain
atoms_per_chain_surr=552 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=55 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''soup
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='PNBT'
NDX_FILENAME='PNBTcharge_10_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=552 # total number of atoms in each chain
atoms_per_chain_surr=57 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=55 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''

'''bulk
EQ_MODEL='eq_3'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPDPP'
NDX_FILENAME='DPPDPPcharge_10'
INPUT_TRJ='/bulk'
FRAMES=5
THRESHOLD=2
TIME_INTERVAL=25000 #ps between each frame
atoms_per_chain=642 # total number of atoms in each chain
atoms_per_chain_surr=642 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=64 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
'''
#'''soup
EQ_MODEL='300_'
TRJ_TYPE='.xtc'
OLIGOMERNAME='DPPDPP'
NDX_FILENAME='DPPDPPcharge_10_1'
INPUT_TRJ='/soup_mono/900/300'
FRAMES=250
TIME_INTERVAL=5000 #ps between each frame
FRAMES=250
THRESHOLD=2
atoms_per_chain=642 # total number of atoms in each chain
atoms_per_chain_surr=66 # total number of atoms in each chain of the surrondin$
atoms_per_0chargegroup=64 # total number of atoms in each monomer to be conside$
                    # in case of IDTBT, I consider conjogated part plus the fir$
repeatunits_per_chain=10 # 5mer=5, 10mer=10, 20mer=20, ...
atoms_of_endgroups=2 # the end gropu of chains are two C united atoms
threshold_mon = float(THRESHOLD) # the threshold for finding monomers around chains in [nm]
largest_dist=3.0 #nm largest distance of atoms in one monomer (0charge group)
SUBMITJOB='yes'
#'''
