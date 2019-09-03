from go import *
import cmbtools
import fourier_tools_py3.fourier_filter as Filter
reload(Filter)
import davetools as dt
reload(dt)

def shell_average(power,debug=1,filename='PowerSpectrum.h5'):
    ff = Filter.FourierFilter(power)
    print("POWER MAX",np.abs(power).max())
    #power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    stuff=[]
    if 1:
        power_1d=np.zeros(ff.nx,dtype='complex128')
        empty = np.zeros_like(power)
        for bbb in range(ff.nx):
            this_flag = ff.get_shell(bbb)
            this_shell=power[this_flag]
            power_1d[bbb] = this_shell.sum(dtype='complex128')
            empty[this_flag] += 1
            if np.abs(this_shell).max() > 1e-6:
                stuff.append(this_shell)
                print("bin %2d shellmax %0.2e power %0.2e"%(bbb,np.abs(this_shell).max(),power_1d[bbb]))
    fptr = h5py.File(filename,'w')

    fptr.create_dataset('power',power_1d.shape,data=power_1d,dtype=power_1d.dtype)
    kspace=ff.get_shell_k()
    fptr.create_dataset('k',kspace.shape,data=kspace)
    fptr.close()
    return ff,stuff

def Make_other(arr, BoxSize=1):
    N = np.array(arr.shape,dtype = np.int32)
    arr = np.ascontiguousarray(arr,dtype=np.double)
    xsize = N.max() #5 * np.pi / 180*BoxSize
    size2d = np.array([xsize,xsize])
    Delta = size2d/N
    Deltal = cmbtools.Delta2l(Delta,N)
    harm = cmbtools.map2harm(arr,Delta)
    lmax = Deltal[0]*N[0]/2
    lbins = np.arange(N[0]//2+1)*Deltal[0] #np.linspace(0,lmax,N//2)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    ClBB = cmbtools.harm2cl(harm,Deltal,lbins)
    output={'Delta':Delta,'Deltal':Deltal,'harm':harm,'lbins':lbins,'ClBB':ClBB,'lmax':lmax}
    return output

dx=1./16
x,y=np.mgrid[0:1:dx,0:1:dx]
r=np.sqrt( (x)**2+(y)**2)
kx=5;ky=2
rho = np.cos(kx*2*np.pi*x)*np.sin(ky*2*np.pi*y)# np.ones(x.shape,dtype=np.double)*0.1
rho = np.ascontiguousarray(rho, dtype=np.double)
rhohat1 = np.fft.fftn(rho)
rhohat1 = rhohat1*rhohat1.conj()

rs=Make_other(rho)
#rho[ r<0.5] = 1.0
#rhohat=cmbtools.map2harm(rho, Delta)
#rhx,y,ohat = np.zeros_like(rho)+1.0
#rhx,y,ohat[r<0.5] = 10
#rho = np.ascontiguousarray(rho, dtype=np.double)
#rhohat = np.ascontiguousarray(rhohat, dtype=np.double)
plt.close('all')
if 1:
    #image the two harmonic spaces
    fig,ax=plt.subplots(1,2)
    ax0=ax[0];ax1=ax[1]
    ppp=ax0.pcolormesh(np.abs(rhohat1))
    ppp=ax1.pcolormesh(np.abs(rs['harm']))
    fig.savefig('hats.png')

if 1:
    #plot the shell averaged spaces
    fig,ax=plt.subplots(1,1)
    ff,empty=shell_average(rhohat1, filename="./spectra_hat.h5")
    k, p = dt.dpy("./spectra_hat.h5",['k','power'])
    p=np.abs(p)
    ok=slice(None)#k>0
    kok=k[ok]*np.pi*2
    
    
    
    ax.plot( kok,p[ok]/p[ok].max(), c='k',marker='o',label='Spec_tool' )
    the_y=rs['ClBB']
    the_x = rs['lbins'][:-1]
    ax.plot(the_x , the_y/the_y.max(),c='b',marker='o',label='cmbtool')
    ax.legend(loc=0)
    dt.axbonk(ax,xlabel='k',ylabel='rho(||k||)',yscale='log')
    ax.set_xscale('symlog',linthreshx=1)
    ax.set_xlim(0,max([the_x.max(),kok.max()]))

    fig.savefig('shells.png')



"""

shell_average(rhohat, filename="./spectra_hat.h5")
#plave(rhohat,"rhohat.png")
#plt.scatter(x.flatten(),y.flatten(),c='r')
if 0:
    ff = Filter.FourierFilter(rhohat)
    rmap = dt.rainbow_map(ff.nx)
    for bbb in range(ff.nx):
        these = ff.get_shell(bbb)
        ax.scatter( x[these], y[these],c=rmap(bbb,these.sum()))
        print(bbb)
    #power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
fig.savefig("rhohat.png")
plt.close(fig)
#MakeDensitySpectra(rho)

plt.clf()
plt.plot( k[ok],p[ok]/kok, c='k',marker='o',label='Spec_tool' )
#plt.plot( kok,kok*256*2*np.pi, c='k',marker='o',label='Spec_tool' )
if 0:
    BoxSize=1
    N=np.array(np.shape(rho), dtype=np.int32)
    xsize = 5 * np.pi / 180*BoxSize
    xsize = 1.0
    size2d = np.array([xsize,xsize])
    Delta = size2d/N
#Deltal = cmbtools.Delta2l(Delta,N)
    Deltal = 0.5*nar([dx,dx])
#harm = cmbtools.map2harm(arr,Delta)
    lmax = Deltal[0]*N[0]
    lbins = np.linspace(0,lmax,0.5/dx)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    ClBB = cmbtools.harm2cl(rhohat,Deltal,lbins)

if 0:
    scale = 1 #p.max()/ClBB.max()
    plt.plot( lbins[:-1], scale*ClBB,marker='o', c='r')

#plt.xscale('log');plt.yscale('log')
plt.savefig('fitting.png')

"""

if 0:
    """I think i need this"""
    if filename is None:
        filename = "%s/power_%s.h5"%(oober.product_dir(frame), field)
