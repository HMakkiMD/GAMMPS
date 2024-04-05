PATH_GAMMPS='/users/hmakki/sharedscratch/GAMMPS/'
PATH_FRAGMENTS='/users/hmakki/sharedscratch/GAMMPS/fragments/'
PATH_OUTPUT='/users/hmakki/sharedscratch/GAMMPS/Polymers/'
PATH_OUTPUT_MD='/users/hmakki/sharedscratch/MD/'
COLOR=['olive', 'navy', 'crimson','orchid','darkorange','tan','darkviolet','pink']
ALLFRAGMENTS=['IDT','T','BTT','DPP','BT','BPD','BTZ','BTZFF','BPT','BO',
                'TO','TCB','TCE','TCN','TCNCN','DPPCC','IG','AIG','TIG', 'BFF',
                    'BFFFF', 'CPDT','TIF', 'NDI', 'Y6','BTPCL','FL','TF','TFF','BTF',
                        'BTFF','BPDF','BPDFF','BTTF','BTTFF','TVT','TVTFF'] #all fragment
ALLFRAGMENTS_MASS=[320.47, 82.12, 138.21, 134.09, 134.16,161.16, 117.11,169.13, 161.16, 120.11, 
                     114.17,98.17, 98.17,109.15,134.16,164.16, 290.32, 292.29, 302.37, 114.09, 
                        150.07, 206.33, 422.60, 294.26, 975.03,1012.79,194.27,102.13,120.12,154.16,
                            172.16, 179.15, 197.14,158.22,176.21,192.3,264.26]
FIVERINGS=['T','TO','TCB','TCE', 'TF', 'TFF','TCN','TCNCN']
SIXRINGCONNECT=['BTZ','IG','AIG','BFF','FL','BT','BTF','BTFF','BPD','BPDF','BPDFF']
THEORY=['b3lyp/3-21g*','b3lyp/6-31g*']
EMTOL='10.0'
RESTRAINT_STRENGTH_HIGH=95000#110000
RESTRAINT_STRENGTH_MID=75000#55000
RESTRAINT_STRENGTH_LOW=22000#25000
VMD_RENDER="'/mnt/data1/users/software/VMD/vmd-1.9.4a38/lib/tachyon/tachyon_LINUXAMD64'"
nochainsinbox = {'1' : 300, '2':150, '5':100, '10':50,'20':25}
CPUNODE='troisi'
GPUNODE='gputroisi'
PolymerLength=10
'''
OLIGOMERNAME = 'DPPTT' # no number is allowed in oligomer name
FRAGMENTLIST = ['BTT','DPP','BTT']#DPPTT

SIDECHAIN_TYPE=['a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
							   # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
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
'''
'''
OLIGOMERNAME = 'DPPDTT' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','BTT']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
							   # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'IDTBO' # no number is allowed in oligomer name
FRAGMENTLIST = ['IDT','BO']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
							   # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'IDTBT' # no number is allowed in oligomer name
FRAGMENTLIST = ['IDT','BT']
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'P3MEEET' # no number is allowed in oligomer name
FRAGMENTLIST = ['TC','TC']
SIDECHAIN_TYPE=['eg','eg'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['3','3'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''
'''
OLIGOMERNAME = 'gT-TT' # no number is allowed in oligomer name
FRAGMENTLIST = ['BTT','TO','TO']
SIDECHAIN_TYPE=['eg','eg'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['3','3'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='yes'
'''
'''
OLIGOMERNAME = 'gDPP-TT' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPPCC','T','T','T']
SIDECHAIN_TYPE=['eg','eg'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['3','3'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''
'''
OLIGOMERNAME = 'gt-gt' # no number is allowed in oligomer name
FRAGMENTLIST = ['TO','TO','TO','TO']
SIDECHAIN_TYPE=['eg','eg','eg','eg'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['2','2','4','4'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'TIGT' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIG','T']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''
'''
OLIGOMERNAME = 'TIGBT' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIG','BT']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'AIGBT' # no number is allowed in oligomer name
FRAGMENTLIST = ['AIG','BT']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''

'''
OLIGOMERNAME = 'CPDTBFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['CPDT','BFF']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['4','4'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['2','2']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
PolymerLength=20
'''
'''
OLIGOMERNAME = 'TIFBT' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIF','BT']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'Y6' # no number is allowed in oligomer name
FRAGMENTLIST = ['Y6']#Y6
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['4','4','9','9'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2','0','0']
SIDECHAIN_BRANCH_LENGTH=['2','2','0','0']
DP=['1']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'BTPCL' # no number is allowed in oligomer name
FRAGMENTLIST = ['BTPCL']#BTPCL
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'DPPDPP' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','T','DPP','T']#BTPCL
SIDECHAIN_TYPE=['eg','eg','a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['3','3','10','10'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['0','0','2','2']
SIDECHAIN_BRANCH_LENGTH=['0','0','8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'PNBT' # no number is allowed in oligomer name
FRAGMENTLIST = ['NDI','T','BT','T']#BTPCL
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['12','12'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['10','10']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'NDITT' # no number is allowed in oligomer name
FRAGMENTLIST = ['T', 'NDI','T']#BTPCL
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'CPDTBT' # no number is allowed in oligomer name
FRAGMENTLIST = ['CPDT','BT']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['14','14'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''
'''
OLIGOMERNAME = 'P3HT' # no number is allowed in oligomer name
FRAGMENTLIST = ['TCB']
SIDECHAIN_TYPE=['a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['4'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'PFO' # no number is allowed in oligomer name
FRAGMENTLIST = ['FL']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['6','6'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'F8BT' # no number is allowed in oligomer name
FRAGMENTLIST = ['FL','BT']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['6','6'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='yes'
'''
'''
OLIGOMERNAME = 'PBTTT' # no number is allowed in oligomer name
FRAGMENTLIST = ['TCB','BTT','TCE']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can b$
                                                           # 'pg' as propylene $
SIDECHAIN_LENGTH=['12','12'] # for alkyl the final length of each chain will be$
SIDECHAIN_BRANCH_POINT=['0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''

'''
OLIGOMERNAME = 'IDTBTF' # no number is allowed in oligomer name
FRAGMENTLIST = ['IDT','BTF']
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'IDTBTFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['IDT','BTFF']
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'IDTBPDF' # no number is allowed in oligomer name
FRAGMENTLIST = ['IDT','BPDF']
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'IDTBPDFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['IDT','BPDFF']
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'TIFBTF' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIF','BTF']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'TIFBTFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIF','BTFF']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'TIFBPD' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIF','BPD']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'TIFBPDF' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIF','BPDF']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'TIFBPDFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIF','BPDFF']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'TIFBO' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIF','BO']#IDTBO
SIDECHAIN_TYPE=['a', 'a', 'a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycol$
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['14','14','14','14'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','0','0']
SIDECHAIN_BRANCH_LENGTH=['0','0','0','0']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'DPPDTTF' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','BTTF']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'DPPDTTFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','BTTFF']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'DPPTVT' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','TVT']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'DPPTVTFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','TVTFF']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''

'''•
OLIGOMERNAME = 'DPPTCN' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','TCN']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''•
OLIGOMERNAME = 'DPPTCNCN' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','TCNCN','T']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'DPPBFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','BFF']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''
#''' step1 not completed from here down
OLIGOMERNAME = 'DPPBFFFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','BFFFF']#DPPDTT
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
#'''
'''•
OLIGOMERNAME = 'DPPBTZFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['T','DPP','T','BTZFF']#DPPBTZ
SIDECHAIN_TYPE=['a', 'a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['16','16','7'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['0','0','1']
SIDECHAIN_BRANCH_LENGTH=['0','0','8']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'DPPTTF' # no number is allowed in oligomer name
FRAGMENTLIST = ['BTTF','DPP','BTTF']#DPPTT

SIDECHAIN_TYPE=['a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'DPPTTFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['BTTFF','DPP','BTTFF']#DPPTT

SIDECHAIN_TYPE=['a', 'a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                               # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='no'
'''

'''
OLIGOMERNAME = 'TIGBTF' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIG','BTF']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'TIGBTFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['TIG','BTF']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['syn']
ROTATEFIRSTMONOMER='no'
'''
'''
OLIGOMERNAME = 'AIGBTF' # no number is allowed in oligomer name
FRAGMENTLIST = ['AIG','BTF']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''
'''
OLIGOMERNAME = 'AIGBTFF' # no number is allowed in oligomer name
FRAGMENTLIST = ['AIG','BTFF']
SIDECHAIN_TYPE=['a','a'] # type of each sidechain should be mentioned. It can be 'a' as alkyl, 'eg' as ethylene glycole,
                                                           # 'pg' as propylene glycol, and 'bg' as butylene glycole
SIDECHAIN_LENGTH=['10','10'] # for alkyl the final length of each chain will be these numbers + 2
SIDECHAIN_BRANCH_POINT=['2','2']
SIDECHAIN_BRANCH_LENGTH=['8','8']
DP=['1','2','5','10','20']
TACTICITY=['iso']
ROTATEFIRSTMONOMER='yes'
'''
HALF_FITTING='no'
MASSDIST='yes'
SUBMITJOB='yes'
