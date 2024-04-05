PATH_GAMMPS='PATH_to_Scripts'
PATH_FRAGMENTS='Path_to_MonomerInputFiles'
PATH_OUTPUT='Path_to_ForceField_Output'
PATH_OUTPUT_MD='Path_to_MD_Files'
COLOR=['olive', 'navy', 'crimson','orchid','darkorange','tan','darkviolet','pink']
ALLFRAGMENTS=['T','DPP','BTZ'] #all fragment
ALLFRAGMENTS_MASS=[82.12, 134.09, 117.11]
FIVERINGS=['T']
SIXRINGCONNECT=['BTZ']
THEORY=['b3lyp/3-21g*','b3lyp/6-31g*']
EMTOL='10.0'
RESTRAINT_STRENGTH_HIGH=95000#110000
RESTRAINT_STRENGTH_MID=55000#55000
RESTRAINT_STRENGTH_LOW=25000#25000
VMD_RENDER="'/mnt/data1/users/software/VMD/vmd-1.9.4a38/lib/tachyon/tachyon_LINUXAMD64'"
nochainsinbox = {'1' : 300, '2':150, '5':100, '10':50,'20':25}
CPUNODE='Name_of_CPUnode'
GPUNODE='Name_of_GPUnode'
PolymerLength=10



OLIGOMERNAME = 'DPPBTZ' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','BTZ']#DPPBTZ
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
							   # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['16','16','7'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','1']
SIDECHAIN_BRANCH_LENGTH=['0','0','8']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'


HALF_FITTING='no'
MASSDIST='yes'
SUBMITJOB='yes'
