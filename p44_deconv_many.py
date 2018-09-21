
nzones_window = 500
width_beam = 2

nzones=300
width_signal = 10.
Norm_Signal = 1

Norm_Powerlaw = 0.2 #0.5
Offset_Powerlaw = 0.2
k_min=2
k_max = nzones
slope = -1

these_x = np.arange(nzones)
window_x = np.arange(nzones_window)
sig_gauss = Norm_Signal*signal.gaussian(nzones,width_signal) #/np.sqrt(2*np.pi*width_signal**2)
sig_powerlaw = Norm_Powerlaw*(p44rp.make_random_power(nzones,k_min,k_max,slope)+Offset_Powerlaw)

sig =sig_powerlaw.real+sig_gauss
win = signal.gaussian(nzones_window,width_beam) #/np.sqrt(2*np.pi*width_beam**2)
filtered = signal.convolve(sig, win, mode='same') / sum(win)


fit_sig = gauss_fit(these_x, sig)
fit_win = gauss_fit(window_x, win)
fit_filt = gauss_fit(these_x, filtered)

deconv = np.sqrt(fit_filt['fit_width']**2-fit_win['fit_width']**2)

