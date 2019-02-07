from pylab import *
import os
import sys
if os.path.isdir('code') :
    sys.path.append('code')
sys.path.append('/Users/dcollins/local-other-2018-01-05/cmbtools_nofits')

from scipy.ndimage.filters import gaussian_filter
from matplotlib.gridspec import *
    
#import cmbtools
import cmbtools_handler as cmbtools

close('all')


def hide_axes():
    gca().get_xaxis().set_visible(False)
    gca().get_yaxis().set_visible(False)

    return()

xmax = 5 * pi/180
ymax = 5 * pi/180


# define r,theta for core of filament

#a = 0.5
#b = 0.2

#a = 0.3 * pi/180
#b = 0.2

a = 1.0 * pi/180
b = 0.15


#thetamin = pi/2
#thetamax = 3*pi
thetamin = -pi/2
thetamax = 3*pi/2
Nsteps = 300

# logarithmic spiral

#theta = thetamin*exp((log(thetamax)-log(thetamin))*linspace(0.0,1.0,Nsteps))
#r = a*exp(b*theta) # logarithmic spiral

r = linspace(a*exp(b*thetamin),a*exp(b*thetamax),Nsteps)
theta = log(r/a)/b if (b!=0.0) else linspace(thetamin,thetamax,Nsteps)



angtan = pi/2 if (b==0.0) else arctan(1/b) # angle between tangent line and radial direction

angpol = theta + angtan - pi/2

xc = 0.5*xmax
yc = 0.5*ymax

x = r*cos(theta) + xc
y = r*sin(theta) + yc

'''
# straight line
x = xc * ones(Nsteps)
y = yc + 1.0*pi/180 * linspace(-1.0,1.0,Nsteps)
angpol = 0.0 * ones(Nsteps)
'''


r0 = 0.1 * pi/180
u = r0*cos(angpol)
v = r0*sin(angpol)


psi = angpol - pi/2


# Set up space
N = array([512,512],dtype = int32)
size2d = array([xmax,ymax])

Delta = size2d/N
Deltal = cmbtools.Delta2l(Delta,N)
lmax = 10000
lbins = linspace(0,lmax,50)
lcent = lbins[:-1] + diff(lbins)/2.

Q = zeros(N)
U = zeros(N)




# identify coordinates of places to fill in
s = 1.0
sres = 20
for i in range(Nsteps) :
    crossx = u[i]*linspace(-s,s,sres)
    crossy = v[i]*linspace(-s,s,sres)
    cox = crossy/5
    coy = -crossx/5

    regx = (x[i] + [[ x0+x1 for x0 in crossx] for x1 in cox]).flatten()
    regy = (y[i] + [[ y0+y1 for y0 in crossy] for y1 in coy]).flatten()

    pixx = array(regx/Delta[1],dtype=int)
    pixy = array(regy/Delta[0],dtype=int)

    for p in zip(pixy,pixx):
        Q[p] = cos(2*psi[i])
        U[p] = sin(2*psi[i])


Q = gaussian_filter(Q,1.0/60.*pi/180/Delta[0],mode='wrap')
U = gaussian_filter(U,1.0/60.*pi/180/Delta[0],mode='wrap')
    
if not Q.flags['C_CONTIGUOUS']:
    Q = np.ascontiguousarray(Q)
if not U.flags['C_CONTIGUOUS']:
    U = np.ascontiguousarray(U)

Qharm = cmbtools.map2harm(Q,Delta)
Uharm = cmbtools.map2harm(U,Delta)

Eharm, Bharm = cmbtools.QU2EB(Qharm,-Uharm,Deltal)

ClEEout = cmbtools.harm2cl(Eharm,Deltal,lbins)
ClBBout = cmbtools.harm2cl(Bharm,Deltal,lbins)


E = cmbtools.harm2map(Eharm,Delta)
B = cmbtools.harm2map(Bharm,Delta)


        
figure(1,figsize=(6,6))

gs= GridSpec(2,2)
gs.update(wspace=0.0,hspace=0.0,left=0.0,right=1.0,top=1.0,bottom=0.0)

subplot(gs[0,0],aspect='equal')

titlex = 0.05
titley = 0.9


jet()
subplot(gs[0])
imshow(Q,origin='lower',extent=[0,xmax,0,ymax])
text(titlex,titley,'Q',transform=gca().transAxes)
#colorbar()
clim(-1,1)
hide_axes()

subplot(gs[1])
text(titlex,titley,'U',transform=gca().transAxes)
imshow(U,origin='lower',extent=[0,xmax,0,ymax])
#colorbar()
clim(-1,1)
hide_axes()



subplot(gs[2])
text(titlex,titley,'E',transform=gca().transAxes)
imshow(E,origin='lower',extent=[0,xmax,0,ymax])
#colorbar()
clim(-1,1)
hide_axes()

subplot(gs[3])
text(titlex,titley,'B',transform=gca().transAxes)
imshow(B,origin='lower',extent=[0,xmax,0,ymax])
#colorbar()
clim(-1,1)
hide_axes()


for i in range(4):
    subplot(gs[i])
    #plot(x,y)
    quivskip = 3
    quiver(x[::quivskip],y[::quivskip],u[::quivskip],v[::quivskip],pivot='middle',headwidth=1.0,headlength=0.0,units='x',scale=1.0,scale_units='x')
    #text(titlex,titley,'Polarization',transform=gca().transAxes)
    
    xlim(0,xmax)
    ylim(0,ymax)
    hide_axes()

savefig('spiral2.pdf',bbox_inches='tight')

show()
