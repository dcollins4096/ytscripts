if 'ef' not in dir():
    execfile('go')
import os
import sys
import pdb
#sys.path.append('cmbtools')
sys.path.append('../cmbtools_nofits')
if os.path.isdir('code') :
    sys.path.append('code')
import p49_fields
import cmbtools
#import pyfits
import astropy.io.fits as pyfits
from pylab import *
import re
import glob
if 'frame' not in dir():
    frame = 0
if 1:
    car = taxi.taxi('p06')
    car.derived_fields['QU']=p49_fields.add_QU
    car.axis=['z']
    car.fields=['Qz']
    car.frames=[frame]
    car.plot()
    Q  = car.the_plot.frb['Qz']
    U  = car.the_plot.frb['Uz']
    car.callbacks=['magnetic_field']
    car.fields=['Qz','density','Uz', 'y-velocity','pressure']
    car.plot()

#Q = array(pyfits.open(Qfile)[0].data,dtype=double)
#U = array(pyfits.open(Ufile)[0].data,dtype=double)
if 0:
    A = 1.0 # 0.00
    B = 0.0 # -1.
    C = 0.5
    offset=0.2
    sigma= 1e-5 #0.03
    noise_amp = 0
    prefix='hor_line_noise_0.1'
    Nx = Ny = 256
    dx = 1./Nx
    dy = 1./Ny
    y,x = np.mgrid[0:1:dx,0:1:dy]
    if 0:
        d = np.abs( A*x + B*y + C)/np.sqrt(A**2+B**2)
        image = np.exp( -(d)**2/(2*sigma**2)) + offset
        image = image/image.max()

        Q = image
        U = image*0
    if 1:
        Q = np.zeros_like(x)
        U = np.zeros_like(x)
        ok = x > 0.4
        ok = np.logical_and(ok, x < 0.6)
        ok = np.logical_and(ok, y > 0.2)
        ok = np.logical_and(ok, y < 0.8 )
        Q[ok] = 1

N = array(shape(Q),dtype = int32)
xsize = 5 * pi / 180
size2d = array([xsize,xsize])
Delta = size2d/N

print("N = ",N)
print("Delta = ",Delta)

Deltal = cmbtools.Delta2l(Delta,N)

Qharm = cmbtools.map2harm(Q,Delta)
Uharm = cmbtools.map2harm(U,Delta)

Eharm, Bharm = cmbtools.QU2EB(Qharm,Uharm,Deltal)

E = cmbtools.harm2map(Eharm,Delta)
B = cmbtools.harm2map(Bharm,Delta)


lmax = Deltal[0]*N[0]
lbins = linspace(0,lmax,100)
lcent = lbins[:-1] + diff(lbins)/2.

ClEE = cmbtools.harm2cl(Eharm,Deltal,lbins)
ClBB = cmbtools.harm2cl(Bharm,Deltal,lbins)

plt.clf()
fig, axes = plt.subplots(2,2,sharex=True,sharey=True)
args={'interpolation':'nearest','origin':'lower'}
axes[0][0].imshow(Q,**args)
axes[0][0].set_title('Q')
axes[1][0].imshow(U,**args)
axes[1][0].set_title("U")
oot=axes[0][1].imshow(E,**args)
fig.colorbar(oot)
axes[0][1].set_title("E")
axes[1][1].imshow(B,**args)
axes[1][1].set_title("B")
outname = 'P49_play_%s_%04d.png'%(car.outname,frame)
fig.savefig(outname)
print(outname )
