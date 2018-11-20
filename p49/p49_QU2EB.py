import os
import sys
import pdb
sys.path.append('cmbtools')
if os.path.isdir('code') :
    sys.path.append('code')

import cmbtools
#import pyfits
import astropy.io.fits as pyfits
from pylab import *
import re
import glob

frbname = "OH NO"
#if len(sys.argv) > 1:
#    rootdir = sys.argv[1]
#    frame = int(sys.argv[2])
#    Qlist = glob.glob(rootdir+'/DD%04d_Q[xyz]*.fits'%frame)
#    print Qlist
#else:
#    rootdir = './output'
#    Qlist = glob.glob(rootdir+'/*/DD*_Q[xyz]*.fits')
def QU2EB(rootdir,frame):
    Qlist = glob.glob(rootdir+'/DD%04d_Q[xyz]*.fits'%frame)
    Ulist = []
    for Qfile in Qlist:
        mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
        Ufile = mo.group(1)+'U'+mo.group(2)+'.fits'
        Ulist.append(Ufile)

    QUlist = zip(Qlist,Ulist)

    for Qfile, Ufile in QUlist :
        print( "WORKING!", Qfile)

        Q = array(pyfits.open(Qfile)[0].data,dtype=double)
        U = array(pyfits.open(Ufile)[0].data,dtype=double)

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


        # write the output near the input
        mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
        outroot = mo.group(1)
        outsuf = mo.group(2)

        Efile = outroot+'E'+outsuf+'.fits'
        Bfile = outroot+'B'+outsuf+'.fits'
        Clfile = outroot+'Cl'+outsuf+'.dat'


        print( Efile, Bfile, Clfile)

        hdu = pyfits.PrimaryHDU(E)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Efile,clobber=True)

        hdu = pyfits.PrimaryHDU(B)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Bfile,clobber=True)

        savetxt(Clfile, zip(lcent,ClEE,ClBB))

'''
# Plot Q/U
figure(figsize=[20,6])
subplot(1,4,1)
title('Q')
imshow(Q,interpolation='nearest')
colorbar()
subplot(1,4,2)
title('U')
imshow(U,interpolation='nearest')
colorbar()

# Plot E/B
subplot(1,4,3)
title('E')
imshow(E,interpolation='nearest')
colorbar()
subplot(1,4,4)
imshow(B,interpolation='nearest')
title('B')
colorbar()

#plot Cl
figure()
plot(lcent,ClEE,label='EE')
plot(lcent,ClBB,label='BB')
loglog()
legend()

show()
'''
'''
N = array([128,128],dtype = int32)
xsize = 1 * pi / 180
size2d = array([xsize,xsize])

Delta = size2d/N

print("N = ",N)
print("Delta = ",Delta)

Deltal = cmbtools.Delta2l(Delta,N)

print("Deltal = ",Deltal)

border = int(N[0]/10)
mask = zeros(N)
mask[border:-border,border:-border] = 1.0
#mask = ones(N)
fmask = sum(mask)/prod(N)

lmax = 10000
lbins = linspace(0,lmax,50)
lcent = lbins[:-1] + diff(lbins)/2.

# Make some E/B

E = gaussian_filter(normal(size=N),2,mode='wrap')
#B = zeros(N)
B = normal(size=N)
B = gaussian_filter(B,4,mode='wrap') - gaussian_filter(B,2,mode='wrap')


#E = zeros(N)#array([[cos(3*2*pi*(i+j)/N[1]) for i in range(N[1]) ] for j in range(N[0])])
#E = array([[cos(3*2*pi*(i+2*j)/N[1]) for i in range(N[1]) ] for j in range(N[0])])
# + sin(2*pi*i/N[1])

#B = array([[cos(3*2*pi*(i+j)/N[1]) for i in range(N[1]) ] for j in range(N[0])]);


# Now do the analysis

Eharm = cmbtools.map2harm(E,Delta)
Bharm = cmbtools.map2harm(B,Delta)

ClEE = cmbtools.harm2cl(Eharm,Deltal,lbins)
ClBB = cmbtools.harm2cl(Bharm,Deltal,lbins)

Qharm, Uharm = cmbtools.EB2QU(Eharm,Bharm,Deltal)

Q = cmbtools.harm2map(Qharm,Delta) * mask
U = cmbtools.harm2map(Uharm,Delta) * mask

Qharm2 = cmbtools.map2harm(Q,Delta)
Uharm2 = cmbtools.map2harm(U,Delta)

Eharm2, Bharm2 = cmbtools.QU2EB(Qharm2,Uharm2,Deltal)

#ClB = cmbtools.harm2cl(Bharm2)
#ClE = cmbtools.harm2cl(Eharm2)
ClEEout = cmbtools.harm2cl(Eharm2,Deltal,lbins)
ClBBout = cmbtools.harm2cl(Bharm2,Deltal,lbins)


E2 = cmbtools.harm2map(Eharm2,Delta)
B2 = cmbtools.harm2map(Bharm2,Delta)


# Plot E/B
figure()
subplot(1,2,1)
title('E')
imshow(E,interpolation='nearest')
colorbar()
subplot(1,2,2)
imshow(B,interpolation='nearest')
title('B')
colorbar()

# Plot Q/U
figure()
subplot(1,2,1)
title('Q')
imshow(Q,interpolation='nearest')
colorbar()
subplot(1,2,2)
title('U')
imshow(U,interpolation='nearest')
colorbar()

# Plot Q/U harmonics
figure()
subplot(2,2,1)
imshow(real(Qharm2), interpolation='nearest')
colorbar()
title('Qharm2')
subplot(2,2,2)
title('Uharm2')
imshow(real(Uharm2), interpolation='nearest')
colorbar()
subplot(2,2,3)
imshow(imag(Qharm2), interpolation='nearest')
colorbar()
subplot(2,2,4)
imshow(imag(Uharm2), interpolation='nearest')
colorbar()


# Plot output E/B harmonics
figure()
subplot(2,2,1)
imshow(real(Eharm2), interpolation='nearest')
colorbar()
title('Eharm2')
subplot(2,2,2)
title('Bharm2')
imshow(real(Bharm2), interpolation='nearest')
colorbar()
subplot(2,2,3)
imshow(imag(Eharm2), interpolation='nearest')
colorbar()
subplot(2,2,4)
imshow(imag(Bharm2), interpolation='nearest')
colorbar()

# Plot output E/B
figure()
subplot(1,2,1)
title('E2')
imshow(E2,interpolation='nearest')
colorbar()
subplot(1,2,2)
title('B2')
imshow(B2,interpolation='nearest')
colorbar()


# Plot input/output spectrum
figure()
plot(lcent,ClEE,label='EE')
plot(lcent,ClEEout/fmask,':',lw=6,label='EEout')
plot(lcent,ClBB,label='BB')
plot(lcent,ClBBout/fmask,':',lw=6,label='BBout')
legend()

show()
'''
