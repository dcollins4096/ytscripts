ef=execfile
execfile('go_lite')
import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
import p44_random_power as p44rp
reload(p44rp)
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
def gauss_kern_1d(size):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    x = np.mgrid[-size:size+1]
    g = np.exp(-(x**2/float(size)))
    return g / g.sum()
def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)

def fwhm(arr, x=None):
    if x is None:
        x = np.arange(arr.size)
    half_max = 0.5*(arr.max())
    over_half = arr > half_max
    max_x = x[over_half].max()
    min_x = x[over_half].min()
    max_pos = x[arr==arr.max()]
    return half_max, max_pos, 1.0*(max_x-min_x)
    
def gauss(x, off,a, b):
    return off * np.exp(-(((x-a)**2)/(2*(b**2))))
def gauss_fit(centers,this_spectra):
    norm = this_spectra.sum()
    vbar= (this_spectra*centers).sum()/norm
    sigma2 = (this_spectra*(centers-vbar)**2).sum()/norm
    sigma = (sigma2)**0.5
    halfmax=this_spectra.max()*0.5
    #plt1.plot([vbar-0.5*sigma,vbar+0.5*sigma],[halfmax,halfmax])
    #ydata = gauss(centers, vbar, sigma)
    popt, pcov = curve_fit(gauss, centers, this_spectra, p0=[this_spectra.max(),vbar,sigma] ) #,  bounds=(0, [(1.1*np.max(this_spectra)), 0.15, np.max(this_spectra)]), method='trf')
    return {'vbar':vbar,'norm':norm,'sigma':sigma,'fit_norm':popt[0],'fit_center':popt[1],'fit_width':popt[2]}

#sig = np.repeat([0., 1., 0.], 100)
nzones=300
nzones_window = 500
width_signal = 10.
width_beam = 20
Norm_Signal = 1
Norm_Powerlaw = 0
Offset_Powerlaw = 0.0
k_min=2
k_max = nzones
slope = -1


these_x = np.arange(nzones)
window_x = np.arange(nzones_window)
sig_gauss = Norm_Signal*signal.gaussian(nzones,width_signal) #/np.sqrt(2*np.pi*width_signal**2)
if 'actual_powerlaw' not in dir() or new_random:
    actual_powerlaw = p44rp.make_random_power(nzones,k_min,k_max,slope).real
    Offset_Powerlaw = -1*min(actual_powerlaw)+Offset_Powerlaw
    
sig_powerlaw = Norm_Powerlaw*(actual_powerlaw) + Offset_Powerlaw
#sig_powerlaw = (p44rp.make_random_power(nzones,2,5,-1))
sig_powerlaw = 0.2

if 0:
    #test powerlaw properties
    k_space = np.fft.fft(sig_powerlaw)
    sig_power = p44rp.spectral_slope(np.abs(k_space[:nzones/2])**2,these_x[:nzones/2],'sig_power')
    #slope of this should be equal to the slope above 
    print sig_power.slope

sig =sig_powerlaw.real+sig_gauss
#win = signal.hann(50)
win = signal.gaussian(nzones_window,width_beam) #/np.sqrt(2*np.pi*width_beam**2)
filtered = signal.convolve(sig, win, mode='same') / sum(win)

half_max_sig, max_pos_sig, full_width_sig = fwhm(sig)
fit_sig = gauss_fit(these_x, sig)
fit_win = gauss_fit(window_x, win)
half_max_win, max_pos_win, full_width_win = fwhm(win)

half_max_filt, max_pos_filt, full_width_filt = fwhm(filtered)
fit_filt = gauss_fit(these_x, filtered)
deconv = np.sqrt(fit_filt['fit_width']**2-fit_win['fit_width']**2)

import matplotlib.pyplot as plt
fig, (ax_orig, ax_win, ax_filt) = plt.subplots(3, 1, sharex=True)


ax_orig.plot(sig_gauss, c='g') #, marker='*')
ax_orig.plot(sig_powerlaw.real, c='b') #, marker='*')
ax_orig.plot(sig,c='k')
#ax_orig.plot( [max_pos_sig-0.5*full_width_sig, max_pos_sig+0.5*full_width_sig], [half_max_sig,half_max_sig], c='b')
#ax_orig.set_title('Original pulse fwhm %0.2f fit %0.2f'%(full_width_sig, fit_sig['fit_width']))
ax_orig.set_title('Gauss width %0.2f Combined Fit %0.2f'%(width_signal, fit_sig['fit_width']))
#ax_orig.set_ylim(0.0,1.0)
ax_orig.margins(0, 0.1)

ax_win.plot(win, marker='*')
ax_win.set_title('Beam width fit fit %0.2f'%( fit_win['fit_width']))
ax_win.margins(0, 0.1)
ax_win.set_ylim(0.0,1.0)

ax_filt.plot(filtered)
ax_filt.plot( [max_pos_filt-0.5*full_width_filt, max_pos_filt+0.5*full_width_filt], [half_max_filt,half_max_filt], c='b')
ax_filt.set_title('Filtered signal, fit %0.2f deconv %0.2f'%(fit_filt['fit_width'], deconv))
ax_filt.margins(0, 0.1)
#ax_filt.set_ylim(0.0,1.0)

nfigs = len(glob.glob('p44_smooth*pdf'))+1
nfigs=1
outname = 'p44_smooth_%d.pdf'%nfigs
fig.tight_layout()
fig.savefig(outname)
print outname

if 1:
    print "fit_sig %0.2f fit_win %0.2f fit_filt %0.2f"%tuple([a['fit_width'] for a in [fit_sig, fit_win, fit_filt]])
    print "width squares",fit_filt['fit_width']**2 , fit_sig['fit_width']**2 , fit_win['fit_width']**2
    print "width error",fit_filt['fit_width']**2 - fit_sig['fit_width']**2 - fit_win['fit_width']**2
    print "width frac",(fit_filt['fit_width']**2 - fit_sig['fit_width']**2 - fit_win['fit_width']**2)/( fit_sig['fit_width']**2 + fit_win['fit_width']**2)

if 0:
    sab = full_width_filt
    sa = full_width_sig
    sb = full_width_win
    print "quadrature ab^2,b^2+a^2", full_width_filt**2, full_width_sig**2+full_width_win**2
    print "error    1 - (b^2+a^2)/ab^2", 1.-(full_width_sig**2+full_width_win**2)/full_width_filt**2, 
