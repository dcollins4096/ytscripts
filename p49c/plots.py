from go import *
import turb_quan
def QUangle(F,Q,U,outname):
    fig,ax=plt.subplots()
    size=F.shape
    dx = 1./nar(size)
    #Healpix: +q is south, +U is south-east, psi is counter-clockwise from +Q
    #However, don't forget that imshow has "Y" on the vertical axis.
    #ax.imshow(F, origin='lower',interpolation='nearest')
    #
    #ax.streamplot(X,Y,Vx,Vy)
    fig.savefig(outname)

    plt.close(fig)
def field_with_streams(F,Vx,Vy,outname):
    fig,ax=plt.subplots()
    size=F.shape
    dx = 1./nar(size)
    Y,X = np.mgrid[0:size[1]:1,0:size[0]:1]
    ax.imshow(F, origin='lower',interpolation='nearest')
    Vxbar = np.mean(Vx)
    varx = Vx-Vxbar
    Vybar = np.mean(Vy)
    vary = Vy-Vybar
    stat(varx)
    ax.streamplot(X,Y,10*varx+0.1*Vxbar,10*vary+0.1*Vybar)
    fig.savefig(outname)
    plt.close(fig)
    turb_quan.plotter2([Vx-np.mean(Vx),Vy-np.mean(Vy)],
                       'p49c_plots/means.png',norm='ind',
                       #axis_labels=stuff['axis_labels'],
                       labs=['x-xbar','y-ybar'])


import matplotlib.colors as colors
def moreplots(stuff):
    flat = {}
    for f in ['n','H2','Hh','Hv']:
        flat[f] = np.sum(stuff[f],axis=0)
        flat[f].shape = (stuff[f].shape[0],stuff[f].shape[1])
    field_with_streams(flat['n'],flat['Hh'],flat['Hv'],stuff,'p49c_plots/with_streams.png')
    fig, axes = plt.subplots(1,1,sharex=True,figsize=(8,8))
    psi=0.5*np.arctan2(stuff['U'],stuff['Q'])
    imax=axes
    b = imax.imshow(flat['n'],origin='lower',interpolation='nearest')#,norm=colors.Normalize(vmin=0,vmax=np.pi))
    psi=0.5*np.arctan2(stuff['U'],stuff['Q'])
    imax.quiver(2*np.cos(psi),2*np.sin(psi))
    psi2=1e3*(psi-np.mean(psi))+np.mean(psi)
    print('arrow',np.std(psi2))
    
    imax.quiver(2*np.cos(psi2),2*np.sin(psi2),color='r',headlength=0,headwidth=0)
    #b = axes[0].imshow(psi-np.mean(psi))#,norm=colors.Normalize(vmin=0,vmax=np.pi))
    #axes[1].hist(psi,histtype='step')
    fig.savefig('p49c_plots/psi.png')


#moreplots(stuff)
