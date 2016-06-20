import numpy as np
nar = np.array
def get_real_data():
    M1x= [1.217665, 1.431674, 1.548269, 1.756108, 2.066012, 2.310685, 2.478404, 2.478714, 2.814167, 2.845678]
    A1x= [1.837977, 1.069593, 1.061287, 0.795001, 0.845782, 0.716979, 0.773236, 0.450656, 0.546341, 0.305136]
    m2x= [1.218471,   1.434830,   1.551471,   1.756267,   2.066233,   2.308028,   2.476005,   2.817269,   2.842953,   ]
    A2x= [0.999267, 0.744210, 0.688218, 0.629504, 0.615768, 0.523427, 0.310401, 0.277059, 0.181711]
    M1 = 10**nar(M1x)
    A1 = 10**nar(A1x)
    return M1,A1

    LogR_12co=[ -0.305786,  -0.275184,  -0.200224,  -0.167759,  -0.030060,  0.067204,   0.052935,   0.178727,   0.330684,   0.348847]   
    LogSigma_12co=[ 3.144261, 2.854906, 2.875071, 2.814946, 2.944083, 2.944222, 3.061680, 2.841230, 2.867236, 2.970413]
    R_12co=10**nar(LogR_12co)
    Sigma_12co=10**nar(LogSigma_12co)

    R_13co=10**(nar([-0.305819,  -0.275338,  -0.202412,  -0.169898,  -0.032177,  0.050769,   0.062983,   0.178708,   0.344724,   0.330593]))
    Sigma_13co=10**(nar([ 3.19010, 3.06980, 3.10142, 2.97253, 3.07015, 3.25651, 3.17917, 2.86701, 3.06782, 2.99331]))

    LadaCoresList = na.array([[ 1,  0.38  , 1.85 ,  0.73],
    [ 2,  0.46  , 2.04 ,  0.66],
    [ 3,  0.49  , 2.12 ,  0.63],
    [ 4,  0.36  , 1.85 ,  0.69],
    [ 5,  0.23  , 1.50 ,  0.82],
    [ 6,  3.14  , 3.57 ,  0.84],
    [ 7,  4.69  , 4.45 ,  0.65],
    [ 8,  3.26  , 3.76 ,  0.75],
    [ 9,  0.56  , 2.27 ,  0.59],
    [10,  0.51  , 2.09 ,  0.68],
    [11,  3.37  , 3.93 ,  0.68],
    [12,  20.37 , 7.06 ,  0.71],
    [13,  0.54  , 2.04 ,  0.78],
    [14,  9.73  , 5.19 ,  0.85],
    [15,  2.64  , 3.58 ,  0.70],
    [16,  3.29  , 4.66 ,  0.40],
    [17,  0.69  , 2.27 ,  0.72],
    [18,  0.35  , 1.82 ,  0.70],
    [19,  0.34  , 1.73 ,  0.79],
    [20,  2.28  , 3.53 ,  0.64],
    [21,  2.66  , 4.29 ,  0.41],
    [22,  1.01  , 2.42 ,  0.86],
    [23,  1.87  , 3.27 ,  0.65],
    [24,  0.70  , 2.29 ,  0.71],
    [25,  1.10  , 2.78 ,  0.63],
    [26,  0.37  , 1.85 ,  0.72],
    [27,  3.09  , 4.37 ,  0.45],
    [28,  0.32  , 1.54 ,  1.07],
    [29,  0.43  , 2.07 ,  0.59],
    [30,  0.41  , 1.99 ,  0.64],
    [31,  1.95  , 3.61 ,  0.50],
    [32,  0.45  , 1.97 ,  0.72],
    [33,  4.27  , 4.37 ,  0.62],
    [34,  2.66  , 3.85 ,  0.57],
    [35,  0.52  , 2.19 ,  0.60],
    [36,  1.69  , 3.16 ,  0.65],
    [37,  1.97  , 2.70 ,  1.22],
    [38,  1.10  , 2.75 ,  0.64],
    [39,  1.07  , 2.51 ,  0.83],
    [40,  9.23  , 5.77 ,  0.59],
    [41,  1.08  , 2.34 ,  1.03],
    [42,  2.79  , 2.96 ,  1.31],
    [43,  0.85  , 2.74 ,  0.51],
    [44,  0.50  , 2.17 ,  0.60],
    [45,  0.64  , 2.47 ,  0.52],
    [46,  0.28  , 1.64 ,  0.79],
    [47,  1.41  , 2.93 ,  0.68],
    [48,  4.18  , 4.80 ,  0.46],
    [49,  0.85  , 2.68 ,  0.54],
    [50,  0.40  , 1.94 ,  0.67],
    [51,  1.20  , 2.87 ,  0.62],
    [52,  0.24  , 1.54 ,  0.82],
    [53,  2.13  , 4.09 ,  0.38],
    [54,  1.40  , 3.03 ,  0.61],
    [55,  0.30  , 1.67 ,  0.80],
    [56,  5.18  , 4.74 ,  0.59],
    [57,  0.31  , 1.60 ,  0.92],
    [58,  0.79  , 2.74 ,  0.47],
    [59,  0.37  , 1.85 ,  0.72],
    [60,  0.43  , 2.02 ,  0.63],
    [61,  2.60  , 3.72 ,  0.62],
    [62,  2.30  , 3.72 ,  0.55],
    [63,  0.36  , 1.70 ,  0.89],
    [64,  0.41  , 1.50 ,  1.47],
    [65,  0.72  , 1.99 ,  1.12],
    [66,  0.98  , 2.34 ,  0.94],
    [67,  2.84  , 4.56 ,  0.37],
    [68,  0.38  , 1.82 ,  0.76],
    [69,  1.98  , 3.43 ,  0.60],
    [70,  1.14  , 2.31 ,  1.12],
    [71,  0.42  , 1.79 ,  0.89],
    [72,  0.72  , 2.17 ,  0.86],
    [73,  0.68  , 2.47 ,  0.55],
    [74,  2.96  , 4.06 ,  0.54],
    [75,  0.25  , 1.50 ,  0.90],
    [76,  0.52  , 2.19 ,  0.60],
    [77,  0.40  , 1.91 ,  0.71],
    [78,  0.34  , 1.85 ,  0.65],
    [79,  1.49  , 3.10 ,  0.61],
    [80,  3.21  , 4.57 ,  0.41],
    [81,  0.43  , 1.88 ,  0.78],
    [82,  0.44  , 1.99 ,  0.68],
    [83,  0.76  , 2.38 ,  0.68],
    [84,  0.49  , 2.02 ,  0.73],
    [85,  0.66  , 2.49 ,  0.52],
    [86,  1.12  , 2.42 ,  0.96],
    [87,  10.29 , 5.24 ,  0.87],
    [88,  2.25  , 3.63 ,  0.57],
    [89,  1.36  , 2.62 ,  0.93],
    [90,  0.49  , 2.02 ,  0.72],
    [91,  1.09  , 2.19 ,  1.26],
    [92,  1.61  , 2.74 ,  0.95],
    [93,  3.55  , 3.66 ,  0.88],
    [94,  1.06  , 2.51 ,  0.82],
    [95,  0.70  , 1.97 ,  1.12],
    [96,  1.11  , 2.45 ,  0.92],
    [97,  5.86  , 5.68 ,  0.39],
    [98,  1.34  , 2.72 ,  0.81],
    [99,  2.22  , 3.10 ,  0.91],
    [100,   0.61 ,  2.31,   0.60],
    [101,   1.87 ,  2.58,   1.34],
    [102,   6.71 ,  5.79,   0.42],
    [103,   0.27 ,  1.54,   0.89],
    [104,   0.53 ,  2.09,   0.71],
    [105,   1.64 ,  2.89,   0.83],
    [106,   0.83 ,  2.02,   1.24],
    [107,   0.46 ,  1.94,   0.77],
    [108,   0.78 ,  2.49,   0.62],
    [109,   3.63 ,  3.72,   0.86],
    [110,   0.37 ,  1.76,   0.82],
    [111,   0.22 ,  1.50,   0.79],
    [112,   1.59 ,  2.38,   1.43],
    [113,   2.39 ,  3.07,   1.01],
    [114,   1.14 ,  2.95,   0.55],
    [115,   0.89 ,  2.58,   0.64],
    [116,   1.20 ,  2.70,   0.75],
    [117,   0.58 ,  1.97,   0.93],
    [118,  0.62  , 2.07 ,  0.86],
    [119,  0.88  , 2.56 ,  0.64],
    [120,  0.42  , 1.88 ,  0.77],
    [121,  2.15  , 4.25 ,  0.34],
    [122,  1.34  , 2.89 ,  0.68],
    [123,  1.55  , 2.96 ,  0.73],
    [124,  0.34  , 1.82 ,  0.69],
    [125,  0.26  , 1.64 ,  0.72],
    [126,  1.50  , 3.35 ,  0.48],
    [127,  1.49  , 2.89 ,  0.75],
    [128,  0.27  , 1.50 ,  0.96],
    [129,  0.36  , 1.88 ,  0.66],
    [130,  0.75  , 2.22 ,  0.84],
    [131,  2.91  , 4.07 ,  0.53],
    [132,  4.67  , 4.76 ,  0.53],
    [133,  1.94  , 3.66 ,  0.48],
    [134,  2.19  , 4.24 ,  0.35],
    [135,  0.44  , 1.97 ,  0.71],
    [136,  1.79  , 3.63 ,  0.46],
    [137,  0.30  , 1.60 ,  0.88],
    [138,  0.26  , 1.64 ,  0.73],
    [139,  1.68  , 3.53 ,  0.47],
    [140,  1.07  , 2.87 ,  0.55],
    [141,  0.75  , 2.42 ,  0.64],
    [142,  0.39  , 1.79 ,  0.82],
    [143,  0.35  , 1.76 ,  0.78],
    [144,  0.31  , 1.57 ,  0.98],
    [145,  0.56  , 2.09 ,  0.74],
    [146,  0.71  , 2.31 ,  0.69],
    [147,  0.37  , 1.79 ,  0.78],
    [148,  0.45  , 1.91 ,  0.79],
    [149,  0.39  , 1.76 ,  0.86],
    [150,  0.81  , 2.51 ,  0.62],
    [151,  1.35  , 2.91 ,  0.67],
    [152,  0.46  , 1.94 ,  0.77],
    [153,  0.83  , 2.66 ,  0.54],
    [154,  0.60  , 2.09 ,  0.80],
    [155,  2.22  , 3.94 ,  0.44],
    [156,  0.27  , 1.57 ,  0.85],
    [157,  0.76  , 2.22 ,  0.84],
    [158,  1.76  , 3.35 ,  0.57],
    [159,  0.84  , 2.75 ,  0.49]])
    LadaCores = {'number':LadaCoresList[:,0],'MassMsun':LadaCoresList[:,1],'RadiusCMe17':LadaCoresList[:,2],'DensityPerCCe4':LadaCoresList[:,3]}
