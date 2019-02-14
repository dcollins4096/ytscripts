import unit_tests_EB as ut
reload(ut)
import numpy as np
import copy
this_spiral = ut.SpiralQU(Nzones=64)

#qswap = ut.Spiral['Q'].swapaxes(0,1)
#uswap = ut.Spiral['U'].swapaxes(0,1)

#qswap = this_spiral['Q'].swapaxes(0,1)
#uswap = this_spiral['U'].swapaxes(0,1)
def dostuff(array):
    b = np.zeros_like(array)
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            b[j,i]=array[i,j]
    b= array.swapaxes(0,1) #np.roll(copy.copy(array.swapaxes(0,1)), 0, axis=0)
    return b
Q1 = this_spiral['Q']
U1 = this_spiral['U']
qswap = dostuff(Q1)
uswap = dostuff(U1)
ebswap=ut.QU_TO_EB(qswap,uswap)
ut.PlotQUEB(qswap,uswap,ebswap['E'],ebswap['B'],outname='Spiral_swap.png')
