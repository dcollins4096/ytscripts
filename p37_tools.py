
from matplotlib.backends.backend_agg import FigureCanvasAgg
from scipy.optimize import leastsq
from scipy import signal
import pyfits
#ef('p14_alpha.py')
def gaussian(x,p):
    #p = {mean, width, norm}
    return p[2]/np.sqrt(p[1]*2*np.pi)*np.exp(-( (x-p[0])/p[1])**2 )
def noisy_gaussian(x,p,amplitude):
    noise = np.random.normal(size=x.size)
    return gaussian(x,p)+noise*amplitude

def residuals(p, y, x):
        err = y-gaussian(x,p)
        return err
def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
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


if 0:
  distance = 140. #pc
  resolution = 40.#arcsec
  smooth_pc = distance*resolution/206264. #in pc
  smooth_px = int(smooth_pc/(4.6/Nzones))
  ef('p37_tools.py')
  blur = blur_image(den,smooth_px)
  hdu = pyfits.PrimaryHDU(blur)
  hdulist = pyfits.HDUList([hdu])
  hdulist.writeto('b02_512_smoothed_%04d_density.fits'%(Nzones))
  hdulist.close()
