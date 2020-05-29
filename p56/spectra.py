
import matplotlib.patches as patches
spec={}
spec['dspec']=dpy("./p56/power_density.h5", ['k','power'])

if 1:
    fig, ax1 = plt.subplots(1,1)
    the_k=spec['dspec'][0]
    ax1.plot( the_k, np.abs(spec['dspec'][1]))
    k0 = the_k[5]
    #P0 = spec['dspec'][1][5]
    #Kolmogorov_ish = P0*(the_k[5:10]/k0)**(-5./3)
    #ax0.plot( the_k[5:10], Kolmogorov_ish)
    axbonk(ax1,xscale='log',yscale='log',xlabel='$k$',ylabel=r'$P_\rho(k)$')


    AC = np.fft.irfft(spec['dspec'][1])*128/(np.pi*2)
    the_k=spec['dspec'][0]
    AC = AC[:the_k.size]
    dL = the_k[1:]-the_k[:-1]
    ACcen = np.abs(0.5*(AC[1:]+AC[:-1]))
    Length = (dL*ACcen).sum()/ACcen[0]
    print("LENGTH",Length)

    rect = patches.Rectangle((0,0),Length,ACcen[0],linewidth=1,edgecolor='r',facecolor='r')
    

    # Add the patch to the Axes
    ax1.add_patch(rect)
    the_x = np.linspace(0,1,AC.size)
    this_dumb_line = np.abs(AC)
    ax1.plot(the_x, this_dumb_line,'g:')
    #ax1.plot(the_x, AC.real,'g:')
    #ax1.plot( the_k, AC.imag,'r')

if 1:
    from scipy.optimize import curve_fit
    def density_ac(r, A, r0, gamma):
        return A/(1+(r/r0)**gamma)
    fit_x = the_x[1:]
    fit_ac = AC[1:]
    fit_ac = 1/(1+(fit_x/0.2)**2)+0.2*np.random.random(fit_ac.size)
    popt, pcov = curve_fit(density_ac, fit_x, fit_ac, p0=[1,0.2,2])
    po2 = AC[0]/2, 10, 1
    the_x = np.linspace(0,0.25,100)
    ax1.plot( the_x, density_ac(the_x, *popt),c='k')
    ax1.plot( fit_x, fit_ac,c='g')
    #ax1.plot(the_x, this_y)
axbonk(ax1,xscale='linear',yscale='linear',xlabel='$r$',ylabel=r'$\langle \rho(x)\rho(x+r)\rangle$')#,xlim=[the_k[1],the_x[-1]],ylim=[0,1])
fig.savefig('p56_corr.png')

if 0:
    fig, ax1 = plt.subplots(1,1)
    the_k=spec['dspec'][0]
    ax1.plot(the_k, 1/(1+(the_k/0.2)**2))
    fig.savefig('p56_corr.png')

plt.close('all')
