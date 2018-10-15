"""
Make the relative orientation field, 
\phi = arctan( cross(B, grad(N))/dot(B, grad(N)))
using gaussian derivatives.
This uses scipy, so if you didn't build yt with 
INST_SCIPY=1
or otherwise install it, you can do

% source /home/dcollins4096/local-yt-2015-05-14/bin/activate

which I *believe* will set everything up properly for you.

The gaussian and blur_image functions are unused.

"""

import pyfits
from scipy import signal
def gaussian(x,p):
    #p = {mean, width, norm}
    return p[2]/np.sqrt(p[1]*2*np.pi)*np.exp(-( (x-p[0])/p[1])**2 )
def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    y,x = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()
def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)

def dx_gauss(size, sizey=None):
    """ dG/dx where G is a gaussian with a width *size*"""
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    y,x = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    norm = 1/g.sum()
    return -2*x*norm*g

def dx_image(im, n, ny=None) :
    """ gaussian derivative in x of *im* with a gaussian width of *n*"""
    g = dx_gauss(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)
def dy_gauss(size, sizey=None):
    """ dG/dy where G is a gaussian with a width *size*"""
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    y,x = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    norm = 1/g.sum()
    return -2*y*norm*g

def dy_image(im, n, ny=None) :
    """ gaussian derivative in y of *im* with a gaussian width of *n*"""
    g = dy_gauss(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)



def phi(density_frb, magnetic_frb,n_smooth=5, axis=2):
    """Compute the field of 
    phi = arctan( B \cross \nabla N/ B \dot \nabla N)
    Needs: to be generalized a bit, presently assumes projection along 
    \hat{z}"""
    density =  density_frb['density'].v
    B1 = [None,None,'Bx'][axis]
    B2 = [None,None,'By'][axis]
    bx_full  = magnetic_frb[B1].v
    by_full  = magnetic_frb[B2].v
    sl=slice(n_smooth,-n_smooth)
    bx = bx_full[sl,sl]
    by = by_full[sl,sl]

    gradN_x = dx_image(density,n_smooth)
    gradN_y = dy_image(density,n_smooth)
    phi = np.zeros_like(density)
    tan_phi = (bx*gradN_y-by*gradN_x)/(bx*gradN_x+by*gradN_y)
    phi_work = np.arctan( tan_phi)
    phi[sl,sl] = phi_work
    return phi/np.pi*180


ef('p20_stokes.py')
def make_plotses(frb_density, frb_mag, prefix, n_smooth=5,axis=2):
    sl=slice(n_smooth,-n_smooth)
    nplots=0
    density = frb_density['density'].v
    nx,ny = density.shape
    B1 = 'Bx'; B2 = 'By'
    #mag = frb['magnetic_energy'].v
    #plt.imshow(dx_image(density,n_smooth),origin='lower',interpolation='nearest')
    #plt.colorbar()
    #plt.savefig('%s_dx.png'%prefix); nplots+=1
    #plt.clf()
    #plt.imshow(dy_image(density,n_smooth),origin='lower',interpolation='nearest')
    #plt.colorbar()
    #plt.savefig('%s_dy.png'%prefix); nplots+=1
    #plt.clf()
    #plt.imshow(density[sl,sl],origin='lower',interpolation='nearest')
    #plt.savefig('%s_d.png'%prefix); nplots += 1
    #plt.clf()
    #plt.imshow(mag[sl,sl],origin='lower',interpolation='nearest')
    #plt.savefig('%s_b2.png'%prefix); nplots += 1
    if 0:
        plt.clf()
        plt.imshow(phi(frb_density, frb_mag,3,axis), origin='lower', interpolation='nearest',cmap='hsv')
        plt.colorbar()
        plt.savefig('%s_phi.png'%prefix); nplots += 1
        plt.clf()
        plt.imshow(density,origin='lower',interpolation='nearest',cmap='gray')
        plt.colorbar()
    Y,X = np.mgrid[0:nx:1, 0:ny:1]
    #plt.streamplot(X,Y, frb_mag[B1].v, frb_mag[B2].v,color='b')

    if 1:
        plt.clf()
        Q = frb_density['Q']
        U = frb_density['U']
        theta = 0.5*np.arctan(U/Q)*180/np.pi
        plt.imshow(theta, origin='lower', interpolation='nearest',cmap='hsv')
        plt.colorbar()
        plt.savefig('%s_stokes.png'%prefix);nplots+=1


    if 0:
        dx_work = dx_image(density,n_smooth)
        dy_work = dy_image(density,n_smooth)
        dx_full = np.zeros_like(density)
        dx_full[sl,sl]=dx_work
        dy_full = np.zeros_like(density)
        dy_full[sl,sl]=dy_work
        plt.clf()
    if 0:
        theta = np.arctan(dx_full/dy_full)/np.pi*180
    if 1:
        #theta = np.arcsin(dy_full/np.sqrt(dx_full**2+dy_full**2))/np.pi*180:
        theta = dy_full/np.sqrt(dx_full**2+dy_full**2)
        theta = np.sqrt(dx_full**2+dy_full**2)


    if 0:
        plt.imshow(theta, 
                   origin='lower', interpolation='nearest',cmap='hsv')
        plt.colorbar()
        plt.savefig("%s_theta_grad.png"%prefix); nplots+=1
        plt.clf()
        plt.imshow(np.arctan(frb_mag[B1].v/frb_mag[B2].v)/np.pi*180, 
                   origin='lower', interpolation='nearest',cmap='hsv')
        plt.colorbar()
        plt.savefig("%s_theta_field.png"%prefix); nplots+=1


    if 0:
        plt.streamplot(X,Y,dx_full,dy_full, color='y')
        plt.xlim(0,nx)
        plt.ylim(0,ny)
        plt.savefig('%s_dq.png'%prefix); nplots+=1

    if 0:
        plt.clf()
        plt.imshow(dx_full, origin='lower',interpolation='nearest',cmap='gray')
        plt.colorbar()
        plt.savefig('%s_dx.png'%prefix); nplots+=1
        plt.clf()
        plt.imshow(dy_full, origin='lower',interpolation='nearest',cmap='gray')
        plt.colorbar()
        plt.savefig('%s_dy.png'%prefix); nplots+=1

    print "Plotted %d plots %s"%(nplots, prefix)
        

if 0:
    basedir ='/scratch1/dcollins/Paper20/a02_sphere_nofield'
    basedir ='/scratch1/dcollins/Paper20/a03_sphere_b2'
    basedir ='/scratch1/dcollins/Paper20/a04_sphere_strongb'
    basedir = '/scratch1/dcollins/Paper20/a05_OT'
    if 'frame' not in dir():
        frame = 100
    prefix = 'a05p_%04d'%frame
    setname = '%s/DD%04d/data%04d'%(basedir,frame,frame)
    line_of_sight = 2
    ds = yt.load(setname)
    proj= yt.ProjectionPlot(ds,line_of_sight,'density')
    proj.annotate_streamlines('Bx','By')
    print proj.save(prefix)
    #proj= yt.ProjectionPlot(ds,2,'magnetic_energy')
    shape = 128
    proj_den = ds.proj('density',line_of_sight) 
    frb_den  = proj_den.to_frb(1,[128,128])
    proj_mag = ds.proj('Bx',line_of_sight, weight_field = 'cell_mass') 
    frb_mag  = proj_mag.to_frb(1,[128,128])
    make_plotses(frb_den,frb_mag,prefix,nx=5,axis= line_of_sight)
    plt.clf()

if 0:
    basedir ='/scratch1/dcollins/Paper20/a02_sphere_nofield'
    basedir ='/scratch1/dcollins/Paper20/a03_sphere_b2'
    shape = 128
    prefix = 'a03s_%04d'%frame
    setname = '%s/DD%04d/data%04d'%(basedir,frame,frame)
    ds = yt.load(setname)
    for z in [0.5]:
        sp = yt.SlicePlot(ds,2,'density', center = [0.5,0.5,z])
        sp.annotate_streamlines('Bx','By')
        slice_prefix="%s_z%0.2f"%(prefix,z)
        sp.save(slice_prefix)
        frb=ds.slice('z',z).to_frb(1,[shape]*2)
        make_plotses(frb,slice_prefix,nx=5)


