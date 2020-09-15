import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
def make_k_freqs(nk,real=False, d=1):

    ny = nk
    nz = nk
    if real:
        nz = nk//2+1

    k_freq = np.zeros([3,nk,ny,nz])
    k1=np.fft.fftfreq(nk,d=d)
    #kx, ky, kz = np.meshgrid(k1,k1,k1)
    #k_freq[0,...]=kx
    #k_freq[1,...]=ky
    #k_freq[2,...]=kz
    x = np.repeat(k1,nk*nk)
    x.shape = (nk,nk,nk)
    y = np.repeat(k1,nk*nk)
    y.shape = (nk,nk,nk)
    y=y.swapaxes(0,1)
    z = np.repeat(k1,nk*nk)
    z.shape = (nk,nk,nk)
    z=z.swapaxes(0,2)
    if real:
        x = x[:,:,:nz]
        y = y[:,:,:nz]
        z = z[:,:,:nz]
    k_freq[0,...]=x
    k_freq[1,...]=y
    k_freq[2,...]=z
    return k_freq

if 0:
    fig,axes=plt.subplots(1,2)
    ax, ax1=axes
    ax.plot(a,'k')
    ax.plot(cr,'r')

    #ax1.plot(cr_hat_1,'k')
    #ax1.plot(cr_hat_1*cr_hat_1.conj(),'r')

fig,axes=plt.subplots(1,2)
ax0, ax1=axes

import scipy.signal
if 0:
    """2d circle"""
    dx=0.1
    x,y = np.mgrid[-1:1:dx,-1:1:dx] #,0:1:dx]
    rho1 = np.zeros_like(x)
    rho2 = np.zeros_like(x)
    x1=0
    y1=0
    r1=np.sqrt( (x-x1)**2+(y-y1)**2)
    rho1[ r1<0.25] = 1
    x2=0.4
    y2=0
    r2=np.sqrt( (x-x2)**2+(y-y2)**2)
    rho2[ r2<0.25] = 2
    #rho1[1,1]=1
    #rho1[1,2]=1

    fig2,axes2=plt.subplots(2,2)
    a20,a21=axes2[0]
    a23,a24=axes2[1]
    a20.imshow( rho2 )
    AC2d=scipy.signal.correlate(rho1,rho2,mode='full',method='direct')
    a21.imshow(AC2d.real)

    rhohat = np.fft.ifftn(rho1)
    rhohatdag = rhohat.conj()
    AC2dft = np.fft.fftn(rhohat*rhohatdag)

    FreqCompRows = np.fft.fftfreq(AC2dft.shape[0])#,d=2)
    FreqCompCols = np.fft.fftfreq(AC2dft.shape[1])#,d=2)
    Kx = np.r_[ FreqCompRows.size*[FreqCompCols]]
    Ky = np.r_[ FreqCompCols.size*[FreqCompRows]].transpose()
    Kr = np.sqrt( Kx**2 + Ky**2)

    a23.imshow(AC2dft.real)
    a24.scatter( Kr, AC2dft.real)

    fig2.savefig('../PigPen/slice1.png')
    print('saved')

if 0:
    """3d sphere"""
    axis=2
    import tools_turb.multi_imshow
    dx=0.05
    axis = 0
    x,y,z = np.mgrid[-1:1:dx,-1:1:dx, -1:1:dx] 

    rho1 = np.zeros_like(x)
    rho2 = np.zeros_like(x)
    x1=0
    y1=0
    z1=0
    r1=np.sqrt( (x-x1)**2+(y-y1)**2 + (z-z1)**2)
    rho1[ r1<0.25] = 1
    x2=0
    y2=0
    z2=0
    r2=np.sqrt( (x-x2)**2+(y-y2)**2 + (z-z2)**2)
    rho2[ r2<0.25] = 1

if 1:
    if 'u05' not in dir():
        u05=taxi.load('u05')
        ds=u05.load(0)
        u05.make_cg(0)

    rho1 = u05.cg['density'].v -1 #[:40,:40,:40]
    rho2=rho1

if 1:
    fig2,axes2=plt.subplots(2,2,figsize=figsize)
    a20,a21=axes2[0]
    a23,a24_dep=axes2[1]
    a20.imshow( rho2.sum(axis=axis), cmap='Greys' )

    fig3,a24=plt.subplots(1,1,figsize=figsize)
if 1:
    if 'AC3d' not in dir() or True:
        print('Correlate')
        AC3d=scipy.signal.correlate(rho1,rho2,mode='same',method='fft')
        AC3d=np.roll(AC3d, AC3d.shape[0]//2,axis=0)
        AC3d=np.roll(AC3d, AC3d.shape[1]//2,axis=1)
        AC3d=np.roll(AC3d, AC3d.shape[2]//2,axis=2)
        print('rolled')
    a21.imshow( (AC3d.real).sum(axis=axis) , cmap='Greys')
              

    print( 'fft')
    rhohat = np.fft.ifftn(rho1)
    rhohatdag = rhohat.conj()
    rho_prod = rhohat*rhohatdag
    print( 'fft2')
    AC3dft = np.fft.fftn(rho_prod)

    k_array = make_k_freqs( AC3d.shape[0],real=False)
    kmag = np.sqrt(k_array[0,...]**2 + k_array[1,...]**2 + k_array[2,...]**2)


if 1:    
    proj = AC3dft.real.sum(axis=axis)
    #proj     = np.roll(proj, proj.shape[0]//2,axis=0)
    #rollproj = np.roll(proj, proj.shape[1]//2,axis=1)
    a23.imshow(proj,cmap='Greys')



    ktt=kmag.flatten()
    R_sphere = 0.25

    hh = R_sphere-ktt
    ok = hh>0
    hh = hh[ok]
    AC_from_caps = (np.pi/3*hh**2*(3*R_sphere-hh))


    import radial_binner as rb
    reload(rb)

    #a24.plot( ktt[ok], AC_from_caps,'k',label='Cap')
    #a24.plot( kmag.flatten(), AC3dft.real.flatten(),label='ACfft',c='b')
    #a24.plot( kmag.flatten(), AC3d.real.flatten(),label='ACfft',c='b')
if 1:    

    bins = np.fft.fftfreq(AC3dft.shape[0])
    ok = bins > 0
    bins = np.r_[0,bins[ok]]
    db = bins[1:]-bins[:-1]
    binned=rb.rb( k_array, AC3dft.real,bins=bins)
    a24.plot( binned[0],binned[1],c='r',label='binned ACfft')

    binned2=rb.rb( k_array, AC3d.real/AC3d.size,bins=bins)
    print('butts')
    a24.plot( binned2[0],binned2[1],c='g',label='binned AC')

    AC = binned2[1]
    L = np.sum(AC*db)/AC[0]
    rect=patches.Rectangle((0,0),L,AC[0],facecolor=[0.5]*3)
    a24.add_patch(rect)
    #a24.set_yscale('log')
    #x0a24.set_xscale('log')
    

if 1:
    a24.legend(loc=0)

    fig2.savefig('../PigPen/slice_3d.png')
    print('saved')


    

if 0:
    ax0.clear()
    ax1.clear()
    ax0.plot(X,'k')
    ax0.plot(H,'b')
    new_X = np.hstack( (X[1:],X))
    AC = scipy.correlate(new_X, H)
    ax0.plot(AC,c='r')
    
    x = np.fft.ifft(X)
    h = np.fft.ifft(H)
    hdag = h.conj()
    XH = np.fft.fft(x*hdag)
    ax1.plot(x,'k')
    ax1.plot(h,'b')
    ax1.plot(x*h,'r')
    ax1.plot(XH,'k:')
    AC_from_ft = XH*Nx
    renorm=(XH*Nx).max()/AC.max()
    renorm=0.5
    AC_from_ft = AC_from_ft/renorm
    ax0.plot(AC_from_ft,'k:')
    



fig2.savefig('../PigPen/p56_convolve_4.png')
fig3.savefig('../PigPen/p56_convolve_5.png')

