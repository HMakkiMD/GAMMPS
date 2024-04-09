%nproc=10
%mem=9GB
%chk=TIFBT-chelpg.chk
#p opt b3lyp/6-31g* pop=chelpg

Gaussian optimization calculation on oligomer4

0,1
C        -17.12705       -0.26280        0.00161
C        -17.90370       -1.47697        0.00172
C        -17.84746        0.92121        0.00196
N        -17.43289       -2.74053        0.00143
N        -19.94931       -2.63365        0.00219
S        -18.73976       -3.75584        0.00176
H        -19.74913        1.94194        0.00262
H        -17.31114        1.86294        0.00192
C        -19.27025        0.97014        0.00238
C        -19.36075       -1.41641        0.00211
C        -20.03350       -0.16806        0.00243
H        -21.11408       -0.14458        0.00271
C         -0.40641        0.91391       -0.00173
C         -1.24510        2.00211       -0.00128
S         -1.32379       -0.60995       -0.00217
C         -7.35068        1.52826       -0.00062
C         -6.65160        0.29613       -0.00117
C         -7.33845       -0.92506       -0.00137
C         -8.72658       -0.89490       -0.00101
C         -9.42559        0.33722       -0.00046
C         -8.73880        1.55839       -0.00027
C         -9.70473       -2.07703       -0.00106
C        -11.05776       -1.35304       -0.00048
C        -10.86746        0.05980       -0.00014
C         -6.37255        2.71040       -0.00048
C         -5.01949        1.98643       -0.00091
C         -5.20997        0.57347       -0.00136
C        -12.32444       -1.89814       -0.00026
C        -13.43997       -1.03042        0.00032
C        -13.23116        0.37925        0.00064
C        -11.94976        0.93527        0.00041
C         -3.75298        2.53150       -0.00088
C         -2.63700        1.66375       -0.00130
C         -2.84631        0.25377       -0.00177
C         -4.12749       -0.30206       -0.00182
C         -9.53548       -2.93560       -1.27912
C         -9.53478       -2.93623        1.27649
C         -6.54235        3.56952       -1.27811
C         -6.54193        3.56899        1.27752
H         -9.27860        2.49940        0.00016
H         -6.79869       -1.86609       -0.00179
H         -3.59783        3.60498       -0.00052
H         -4.27246       -1.37598       -0.00216
H        -12.47966       -2.97160       -0.00053
H        -11.80483        2.00920        0.00063
H         -8.53899       -3.39126       -1.30006
H        -10.28342       -3.73628       -1.29876
H         -9.65908       -2.31473       -2.17090
H         -9.65789       -2.31579        2.16865
H        -10.28271       -3.73692        1.29614
H         -8.53828       -3.39189        1.29667
H         -7.53846        4.02456        1.29838
H         -6.41838        2.94812        2.16932
H         -5.79404        4.36972        1.29723
H         -7.53889        4.02510       -1.29846
H         -5.79445        4.37026       -1.29772
H         -6.41908        2.94901       -2.17019
H        -15.21660       -2.37737        0.00054
H         -0.86287        3.01105       -0.00096
S        -14.75397        1.24269        0.00123
C        -14.83298       -1.36865        0.00066
C        -15.66994       -0.28108        0.00124
C          3.99104        0.77730       -0.00165
C          3.29720        2.03734       -0.00109
C          3.17901       -0.34795       -0.00234
N          3.86668        3.25825       -0.00057
N          1.36240        3.35539       -0.00058
S          2.65623        4.38037       -0.00007
H          1.23338       -1.23535       -0.00294
H          3.63988       -1.32874       -0.00280
C          1.76883       -0.29318       -0.00240
C          1.83564        2.09404       -0.00120
C          1.04643        0.89155       -0.00181
C         20.54884       -1.67529        0.00320
C         19.66059       -2.70499        0.00192
S         19.78010       -0.10283        0.00403
C         13.59685       -1.76041       -0.00055
C         14.38902       -0.58639        0.00080
C         13.79889        0.68383        0.00153
C         12.41224        0.76134        0.00091
C         11.62003       -0.41265       -0.00041
C         12.21045       -1.68319       -0.00117
C         11.52865        2.01579        0.00147
C         10.12349        1.39891        0.00038
C         10.20397       -0.02457       -0.00070
C         14.48066       -3.01487       -0.00113
C         15.88595       -2.39800        0.00013
C         15.80641       -0.97539        0.00122
C          8.90307        2.04048        0.00044
C          7.72314        1.26181       -0.00057
C          7.82260       -0.16009       -0.00161
C          9.05695       -0.81353       -0.00171
C         17.10672       -3.04113        0.00028
C         18.28567       -2.26449        0.00152
C         18.18897       -0.84314        0.00259
C         16.95309       -0.18863        0.00245
C         11.76456        2.86008       -1.27553
C         11.76344        2.85816        1.27995
C         14.24575       -3.85703       -1.27970
C         14.24441       -3.85927        1.27569
H         11.59952       -2.57966       -0.00219
H         14.40986        1.58024        0.00255
H         17.17489       -4.12381       -0.00054
H         16.89353        0.89351        0.00328
H          8.83175        3.12277        0.00126
H          9.11825       -1.89545       -0.00255
H         12.79323        3.23753       -1.29535
H         11.08060        3.71612       -1.29483
H         11.59414        2.25148       -2.16806
H         11.59222        2.24823        2.17141
H         11.07948        3.71418        1.29990
H         12.79211        3.23557        1.30126
H         13.21577       -4.23691        1.29545
H         14.41470       -3.25066        2.16823
H         14.92844       -4.71530        1.29509
H         13.21712       -4.23461       -1.30124
H         14.92978       -4.71303       -1.29985
H         14.41702       -3.24686       -2.17098
H          6.05880        2.74246       -0.00002
H         19.95056       -3.74678        0.00126
S          6.23775       -0.90337       -0.00277
C          6.36167        1.70693       -0.00063
C          5.44111        0.68701       -0.00156
H         21.62589       -1.74272        0.00373

