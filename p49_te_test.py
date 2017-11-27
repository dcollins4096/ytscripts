import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cmbtools
import pylab 
from scipy.ndimage.filters import gaussian_filter
import numpy as np


N = np.array([128,128],dtype = np.int32)
xsize = 1 * np.pi / 180
size2d = np.array([xsize,xsize])

Delta = size2d/N

print("N = ",N)
print("Delta = ",Delta)

Deltal = cmbtools.Delta2l(Delta,N)

print("Deltal = ",Deltal)

border = int(N[0]/10)
mask = np.zeros(N)
mask[border:-border,border:-border] = 1.0
mask = np.ones(N)
fmask = np.sum(mask)/np.prod(N)

lmax = 10000
lbins = np.linspace(0,lmax,50)
lcent = lbins[:-1] + np.diff(lbins)/2.

# Make some E/B

N1 = pylab.normal(size=N)
T = gaussian_filter(N1,2,mode='wrap')
E = gaussian_filter(N1,2,mode='wrap')
B = np.zeros(N)

Eharm = cmbtools.map2harm(E,Delta)
Bharm = cmbtools.map2harm(B,Delta)

ClEE = cmbtools.harm2cl(Eharm,Deltal,lbins)
ClBB = cmbtools.harm2cl(Bharm,Deltal,lbins)

Qharm, Uharm = cmbtools.EB2QU(Eharm,Bharm,Deltal)

Q = cmbtools.harm2map(Qharm,Delta) * mask
U = cmbtools.harm2map(Uharm,Delta) * mask


Tharm = cmbtools.map2harm(T,Delta)
Qharm = cmbtools.map2harm(Q,Delta) 
Uharm = cmbtools.map2harm(U,Delta) 

Eharm, Bharm = cmbtools.QU2EB(Qharm,Uharm,Deltal)

clTTout = cmbtools.harm2cl(Tharm,Deltal,lbins)
clEEout = cmbtools.harm2cl(Eharm,Deltal,lbins)
clBBout = cmbtools.harm2cl(Bharm,Deltal,lbins)
clTEout = cmbtools.harm2clcross_samegrid(Tharm,Eharm, Deltal,lbins)

