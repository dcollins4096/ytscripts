import numpy as np

sfrm = "%12s"
if 'thresh' not in dir():
    thresh = 1e-11
def nz(arr):
    return  np.where(np.abs(arr) > thresh)
def nonzero(arr):
    return arr[ nz(arr)]
def pnz(arr):
    print(arr[ nz(arr)])
    print(nz(arr))
cf = '({0.real:5.2e} + {0.imag:5.2e}i)'
def ampl_str(arr):
    out = ""
    for v in arr:
        out += cf.format(v)
    return out
