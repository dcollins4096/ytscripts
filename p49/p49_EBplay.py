if 'ef' not in dir():
    execfile('go')
import os
import sys
import pdb
#sys.path.append('cmbtools')
#sys.path.append('../cmbtools_nofits')
sys.path.append('/Users/dcollins/local-other-2018-01-05/cmbtools_nofits')
if os.path.isdir('code') :
    sys.path.append('code')
import p49_fields
#import cmbtools
import cmbtools_handler as cmbtools
#import pyfits
import astropy.io.fits as pyfits
from pylab import *
import re
import glob
if 'frame' not in dir():
    frame = 0
EB_KLUDGE = -1
if 1:
    car = taxi.taxi('p12')
    #car.outname = 'p12_field'
    simname = car.outname
    car.derived_fields['QU']=p49_fields.add_QU
    car.axis=['z']
    car.fields=['Qz']
    car.frames=[frame]
    car.plot()
    Q  = car.the_plot.frb['Qz'].v
    U  = car.the_plot.frb['Uz'].v
    car.callbacks=['magnetic_field']
    car.fields=['Qz','density','Uz', 'y-velocity','pressure']
    car.plot()

#Q = array(pyfits.open(Qfile)[0].data,dtype=double)
#U = array(pyfits.open(Ufile)[0].data,dtype=double)
if 0:
    simname = 'spiral'
    xmax = 5 *np.pi/180
    ymax = 5 *np.pi/180
    a = 1.0 * np.pi/180
    b = 0.15
    thetamin = -np.pi/2
    thetamax = 3*np.pi/2
    Nsteps = 300
    r = np.linspace(a*np.exp(b*thetamin),a*np.exp(b*thetamax),Nsteps)
    theta = np.log(r/a)/b if (b!=0.0) else np.linspace(thetamin,thetamax,Nsteps)
    angtan = np.pi/2 if (b==0.0) else np.arctan(1/b) # angle between tangent line and radial direction
    angpol = theta + angtan - np.pi/2
    xc = 0.5*xmax
    yc = 0.5*ymax
    x = r*np.cos(theta) + xc
    y = r*np.sin(theta) + yc
    r0 = 0.1 * pi/180
    u = r0*np.cos(angpol)
    v = r0*np.sin(angpol)
    psi = angpol - np.pi/2
    N = np.array([512,512],dtype = int32)
    size2d = np.array([xmax,ymax])
    Delta = size2d/N
    Deltal = cmbtools.Delta2l(Delta,N)
    lmax = 10000
    lbins = np.linspace(0,lmax,50)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    Q =np.zeros(N)
    U =np.zeros(N)
# identify coordinates of places to fill in
    s = 1.0
    sres = 20
    for i in range(Nsteps) :
        crossx = u[i]*np.linspace(-s,s,sres)
        crossy = v[i]*np.linspace(-s,s,sres)
        cox = crossy/5
        coy = -crossx/5

        regx = (x[i] + [[ x0+x1 for x0 in crossx] for x1 in cox]).flatten()
        regy = (y[i] + [[ y0+y1 for y0 in crossy] for y1 in coy]).flatten()

        pixx = np.array(regx/Delta[1],dtype=int)
        pixy = np.array(regy/Delta[0],dtype=int)

        for p in zip(pixy,pixx):
            Q[p] =np.cos(2*psi[i])
            U[p] =np.sin(2*psi[i])

    from scipy.ndimage.filters import gaussian_filter
    Q = gaussian_filter(Q,1.0/60.*pi/180/Delta[0],mode='wrap')
    U = gaussian_filter(U,1.0/60.*pi/180/Delta[0],mode='wrap')


if 0:
    A = 1.0 # 0.00
    B = -1 # -1.
    C = 0.5
    offset=0.2
    sigma= 1e-5 #0.03
    noise_amp = 0

    Nx = Ny = 32
    dx = 1./Nx
    dy = 1./Ny
    y,x = np.mgrid[0:1:dx,0:1:dy]
    Q = np.zeros_like(x)
    U = np.zeros_like(x)
    if 1:
        pass
    if 0:
        simname = 'h2'
        d = np.abs( A*x + B*y + C)/np.sqrt(A**2+B**2)
        image = np.exp( -(d)**2/(2*sigma**2)) + offset
        image = image/image.max()

        Q = image
        U = image*0
    if 0:
        ok = x > 0.4
        ok = np.logical_and(ok, x < 0.6)
        ok = np.logical_and(ok, y > 0.2)
        ok = np.logical_and(ok, y < 0.8 )
        Q[ok] = 1
    if 0:
        simname = 'point_Q'
        Q[int(Nx/2),int(Ny/2)] = 1
    if 0:
        simname = 'line_horiz_Q'
        Q[int(Nx/2),int(Ny/4):int(3.*Ny/4)] = 1
    if 0:
        simname = 'line_vert_Q'
        Q[int(Nx/4):int(3.*Nx/4),int(Ny/2)] = 1
    if 0:
        simname = 'line_vert_Qm'
        Q[int(Nx/4):int(3.*Nx/4),int(Ny/2)] = -1


N = array(shape(Q),dtype = int32)
xsize = 5 * pi / 180
size2d = array([xsize,xsize])
Delta = size2d/N

print("N = ",N)
print("Delta = ",Delta)

Deltal = cmbtools.Delta2l(Delta,N)
if not Q.flags['C_CONTIGUOUS']:
    Q = np.ascontiguousarray(Q)
if not U.flags['C_CONTIGUOUS']:
    U = np.ascontiguousarray(U)

Qharm = cmbtools.map2harm(Q,Delta)
Uharm = cmbtools.map2harm(U,Delta)

Eharm, Bharm = cmbtools.QU2EB(Qharm,EB_KLUDGE*Uharm,Deltal)

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
outname = '%s_%04d_P49_play.png'%(simname,frame)
fig.savefig(outname)
print(outname )
plt.clf()
plt.imshow(E,**args)
plt.savefig("%s_%04d_E.png"%(simname,frame))
plt.clf()
plt.imshow(B,**args)
plt.savefig("%s_%04d_B.png"%(simname,frame))
