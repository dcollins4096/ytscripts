import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/home/huffenbe/projects/cmbtools')

import cmbtools
#from pylab import *

#close('all')

N = np.array([128,128],dtype = np.int32)
xsize = 1 * np.pi / 180
size2d = np.array([xsize,xsize])

Delta = size2d/N

print("N = ",N)
print("Delta = ",Delta)

Deltal = cmbtools.Delta2l(Delta,N)

print("Deltal = ",Deltal)


E = np.zeros(N)
#E = array([[cos(3*2*pi*(i+2*j)/N[1]) for i in range(N[1]) ] for j in range(N[0])])
# + sin(2*pi*i/N[1])

B = np.array([[np.cos(3*2*np.pi*(i+j)/N[1]) for i in range(N[1]) ] for j in range(N[0])]);

Eharm = cmbtools.map2harm(E,Delta)
Bharm = cmbtools.map2harm(B,Delta)

Qharm, Uharm = cmbtools.EB2QU(Eharm,Bharm,Deltal)

Q = cmbtools.harm2map(Qharm,Delta)
U = cmbtools.harm2map(Uharm,Delta)

Eharm2, Bharm2 = cmbtools.QU2EB(Qharm,Uharm,Deltal)

E2 = cmbtools.harm2map(Eharm2,Delta)
B2 = cmbtools.harm2map(Bharm2,Delta)


#figure()
plt.clf()
plt.imshow(E,interpolation='nearest')
plt.title('E')
plt.colorbar()
plt.savefig('p49_E.png')

plt.clf()
plt.imshow(B,interpolation='nearest')
plt.title('B')
plt.colorbar()
plt.savefig('p49_B.png')

plt.clf()
plt.imshow(Q,interpolation='nearest')
plt.title('Q')
plt.colorbar()
plt.savefig('p49_Q.png')

plt.clf()
plt.imshow(U,interpolation='nearest')
plt.title('U')
plt.colorbar()
plt.savefig('p49_U.png')

plt.clf()
plt.imshow(E2,interpolation='nearest')
plt.title('E2')
plt.colorbar()
plt.savefig('p49_E2.png')

plt.clf()
plt.imshow(B2,interpolation='nearest')
plt.title('B2')
plt.colorbar()
plt.savefig('p49_B2.png')

