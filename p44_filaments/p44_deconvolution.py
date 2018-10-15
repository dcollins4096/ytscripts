ef=execfile
#execfile('go_lite')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
from p44_random_tools import *    

#sig = np.repeat([0., 1., 0.], 100)
beam_list= [1,2,5,8]
c_list = ['r','g','b','c']
nzones=70
nzones_window = 50
width_signal = 10.
width_beam = 8
Norm_Signal = 1

if 'n_plot' not in dir():
    #n_plot dictates which powerlaw to take.
    n_plot=1
powerlaws=[]
Mach = 9.
#var_linear=(0.5 Mach)**2
powerlaws.append(powerlaw_parameters(var=2.0,mean=1.,k_min=2,k_max=nzones,slope=-1))
powerlaws.append(powerlaw_parameters(var=2.0,mean=1.,k_min=2,k_max=nzones,slope=-2))
powerlaws.append(powerlaw_parameters(var=2.0,mean=1.,k_min=2,k_max=nzones,slope=-0.5))
powerlaws.append(powerlaw_parameters(var=1.0,mean=1.,k_min=2,k_max=nzones,slope=-0.5))

this_powerlaw = powerlaws[n_plot]
print this_powerlaw
signal_x = np.arange(nzones)
window_x = np.arange(nzones_window)

if 'deconvolved_widths' not in dir():
    deconvolved_widths ={}


if 'new_random' not in dir():
    new_random=True
if 'phases' not in dir() or new_random:
    phases=None

phases, sig_powerlaw = make_random_power(nzones,this_powerlaw.k_min,this_powerlaw.k_max,this_powerlaw.slope,phases=phases, 
        var = this_powerlaw.var, mean=this_powerlaw.mean)

plt.clf()
plt.hist(sig_powerlaw,histtype='step')

plt.savefig('p44_hist.pdf')
check_powerlaw( sig_powerlaw )


sig_gauss = Norm_Signal*signal.gaussian(nzones,width_signal) 
sig =sig_powerlaw.real+sig_gauss
fit_sig = gauss_fit(signal_x, sig)


fig, (ax_orig, ax_win, ax_filt) = plt.subplots(3, 1, sharex=True)


ax_orig.plot(sig_gauss, c='g') #, marker='*')
ax_orig.plot(sig_powerlaw.real, c='b') #, marker='*')
ax_orig.plot(signal_x,gauss_me(signal_x,fit_sig['fit_norm'],fit_sig['fit_center'],fit_sig['fit_width']),c='r')
ax_orig.plot(sig,c='k')
ax_orig.set_title('Gauss width %0.2f Combined Fit %0.2f'%(width_signal, fit_sig['fit_width']))
ax_orig.margins(0, 0.1)
ax_orig.text(1,1,powerlaws[n_plot])
#powerlaws.append(powerlaw_parameters(norm=2.0,offset=0.1,kmin=2,kmax=nzones,slope=-0.5))

deconvolved_widths[n_plot]=[]

for n_window,width_beam in enumerate(beam_list):
    win = signal.gaussian(nzones_window,width_beam) #/np.sqrt(2*np.pi*width_beam**2)
    filtered = signal.convolve(sig, win, mode='same') / sum(win)

    fit_win = gauss_fit(window_x, win, fit_fwhm_only=False)
    fit_filt = gauss_fit(signal_x, filtered)

    deconvolved_widths[n_plot].append( np.sqrt(fit_filt['fit_width']**2-fit_win['fit_width']**2))

    ax_win.plot(win, c=c_list[n_window])
    #ax_win.set_title('Beam width fit fit %0.2f'%( fit_win['fit_width']))
    ax_win.margins(0, 0.1)
    ax_win.set_ylim(0.0,1.0)

    ax_filt.plot(filtered, c=c_list[n_window])
    ax_filt.plot(signal_x,gauss_me(signal_x,fit_filt['fit_norm'],fit_filt['fit_center'],fit_filt['fit_width']),
            c=c_list[n_window], linestyle='--')
    #ax_filt.set_title('Filtered signal, fit %0.2f deconv %0.2f'%(fit_filt['fit_width'], deconv))
    ax_filt.margins(0, 0.1)

outname = 'p44_smooth_%d.pdf'%n_plot
fig.tight_layout()
fig.savefig(outname)
print outname
plt.close(fig)

fig_deconv, ax_deconv = plt.subplots(1, 1)
ax_deconv.plot(beam_list, width_signal*np.ones(len(beam_list),dtype='float'),label=r'$\sigma_{\rm{signal}}$')
ax_deconv.plot(beam_list, fit_sig['fit_width']*np.ones(len(beam_list),dtype='float'),label=r'$\sigma_{\rm{fit}}$')

deconv_range = None
for n in deconvolved_widths.keys():
    ax_deconv.plot(beam_list, deconvolved_widths[n],label=str(powerlaws[n]))
    deconv_range = collect_extrema(deconvolved_widths[n], deconv_range)
    print deconv_range
ax_deconv.set_xlabel('beam width')
ax_deconv.legend(loc=0)
ax_deconv.set_ylabel('deconvolved width')
ax_deconv.set_ylim(0,np.ceil(deconv_range.max()))
outname = 'p44_beam_vs_deconv.pdf'
fig_deconv.savefig(outname)
plt.close(fig_deconv)
print outname


