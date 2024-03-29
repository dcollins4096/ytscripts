from go import *
import scipy.signal
import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
import p56.radial_binner as rb
reload(rb)
from scipy.optimize import curve_fit
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
def powerlaw2(r, r0,alpha):
    rhosquared=1
    return alpha*np.log10(r/r0) + np.log10(rhosquared)
axis=0
labelmap = {'p52_441':'D22','p52_432':'D26','p52_433':'D27','p52_434':'D28'}
sim_list = ['p52_441','p52_432','p52_433','p52_434']


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


def add_fields(obj):
    def bx_hat(field,data):
        return data['magnetic_field_x']/data['magnetic_field_strength']
    def by_hat(field,data):
        return data['magnetic_field_y']/data['magnetic_field_strength']
    def bz_hat(field,data):
        return data['magnetic_field_z']/data['magnetic_field_strength']
    obj.add_field("bx_hat",bx_hat,units='dimensionless',sampling_type='cell')
    obj.add_field("by_hat",by_hat,units='dimensionless',sampling_type='cell')
    obj.add_field("bz_hat",bz_hat,units='dimensionless',sampling_type='cell')

    def vx_hat(field,data):
        vmag =data['velocity_magnitude']
        out= data['velocity_x']/vmag
        out[ vmag < 1e-16 ] = 0
        return out
    def vy_hat(field,data):
        vmag =data['velocity_magnitude']
        out= data['velocity_y']/vmag
        out[ vmag < 1e-16 ] = 0
        return out
    def vz_hat(field,data):
        vmag =data['velocity_magnitude']
        out= data['velocity_z']/vmag
        out[ vmag < 1e-16 ] = 0
        return out
    obj.add_field("vx_hat",vx_hat,units='dimensionless',sampling_type='cell')
    obj.add_field("vy_hat",vy_hat,units='dimensionless',sampling_type='cell')
    obj.add_field("vz_hat",vz_hat,units='dimensionless',sampling_type='cell')

class make_ac():
    def __init__(self,carname,frame=0,prefix='luge',use='yt',do_mean=False,c=None,l=None):

        self.prefix=prefix
        self.frame=frame
        self.plot_args={'c':c,'linestyle':l}
        if use == 'yt':
            car=taxi.load(carname)
            ds=car.load(frame)
            add_fields(ds)

            car.derived_fields['hats']=add_fields
            car.make_cg(frame)
            cg = car.cg 

            self.bx_hat=cg['bx_hat'].v
            self.by_hat=cg['by_hat'].v
            self.bz_hat=cg['bz_hat'].v
            if do_mean:
                print("DOES NOT WORK to take the mean off the yt reader. Fix it.")
                raise
        if use == 'ytvel':
            car=taxi.load(carname)
            ds=car.load(frame)
            add_fields(ds)

            car.derived_fields['hats']=add_fields
            car.make_cg(frame)
            cg = car.cg 
            self.cg=cg

            self.bx_hat=cg['vx_hat'].v
            self.by_hat=cg['vy_hat'].v
            self.bz_hat=cg['vz_hat'].v

            if do_mean:
                print("DOES NOT WORK to take the mean off the yt reader. Fix it.")
                raise
        elif use == 'fits':
            bxname = "%s%s"%(carname,"magnetic_field_x.fits")
            byname = "%s%s"%(carname,"magnetic_field_y.fits")
            bzname = "%s%s"%(carname,"magnetic_field_z.fits")
            Bx = pyfits.open(bxname)[0].data
            By = pyfits.open(byname)[0].data
            Bz = pyfits.open(bzname)[0].data
            if do_mean:
                Bx_mean = np.mean(Bx)
                By_mean = np.mean(By)
                Bz_mean = np.mean(Bz)
                Bx = Bx - Bx_mean
                By = By - By_mean
                Bz = Bz - Bz_mean
            Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
            self.bx_hat= Bx/Bmag
            self.by_hat= By/Bmag
            self.bz_hat= Bz/Bmag
        elif use == 'test':
            self.bx_hat = np.ones([32,32,32])
            self.by_hat = np.zeros([32,32,32])
            self.bz_hat = np.zeros([32,32,32])


        print('do fft')
        self.bxhathat = np.fft.ifftn(self.bx_hat)
        self.byhathat = np.fft.ifftn(self.by_hat)
        self.bzhathat = np.fft.ifftn(self.bz_hat)
        self.bxhatdag = self.bxhathat.conj()
        self.byhatdag = self.byhathat.conj()
        self.bzhatdag = self.bzhathat.conj()
        size=self.bxhathat.size
        self.ac_bxhat = np.fft.fftn( self.bxhathat*self.bxhatdag)*size
        self.ac_byhat = np.fft.fftn( self.byhathat*self.byhatdag)*size
        self.ac_bzhat = np.fft.fftn( self.bzhathat*self.bzhatdag)*size
        self.AC3dft = self.ac_bxhat+self.ac_byhat+self.ac_bzhat

        bins = np.fft.fftfreq(self.AC3dft.shape[0])
        ok = bins > 0
        bins = np.r_[0,bins[ok]]
        k_array = make_k_freqs( self.AC3dft.shape[0],real=False)
        self.kmag = np.sqrt(k_array[0,...]**2 + k_array[1,...]**2 + k_array[2,...]**2)
        print('do binning')
        self.binned_ac=rb.rb( k_array, self.AC3dft.real/self.AC3dft.size,bins=bins)

    def PUT_COMPUTE_AC_HERE(self):
        pass
    def plot_bxhat(self):
        fig2,axes2=plt.subplots(2,2,figsize=figsize)
        a20,a21=axes2[0]
        a22,a23=axes2[1]

        a20.imshow(self.bx_hat.sum(axis=axis), cmap='Greys' )
        a21.imshow(self.by_hat.sum(axis=axis), cmap='Greys' )
        a22.imshow(self.bz_hat.sum(axis=axis), cmap='Greys' )
        frame=60
        fig2.savefig('p52_%s_bxhats_n%04d.png'%(self.prefix,frame))
        plt.close(fig2)


prefix='434'
if 'aclist' not in dir():
    aclist={}
if 'vel_cor' not in dir():
    vel_cor={}
    

linemap = {'p52_434':'-','p52_441':'--', 'p52_432':':', 'p52_433':'-.','turb':'-'}
colordict={10:'k',50:'r',60:'r'}
clobber=False
for frame in [0]:
    if frame not in aclist:
        aclist[frame]={}
    for carname in ['p52_434']:
        if carname not in aclist[frame] or clobber:
            aclist[frame][carname] =  make_ac(carname=carname,frame=0,prefix=r'$%s\ \mathbf{b}\ t=0s$'%(labelmap[carname]),c='k'
                    ,l=linemap[carname])

for frame in [0,20,30,40,50,60,70,90,100]:
    if frame not in aclist:
        aclist[frame]={}
    if frame not in vel_cor:
        vel_cor[frame]={}

    time = frame/100
    for carname in ['p52_434','p52_441']:#,'p52_432','p52_433']:
        if carname not in aclist[frame] or clobber:
            prefix = r'$%s\ \mathbf{b}\ t=%0.1fs$'%(labelmap[carname],time)
            c=colordict.get(frame,'g')
            if sim == 'turb':
                c='b'
            aclist[frame][carname] =  make_ac(carname=carname,frame=frame,prefix=prefix,c=c,l=linemap[carname])

        if carname not in vel_cor[frame] or clobber:
            prefix = r'$%s\ \mathbf{v}\ t=%0.1fs$'%(labelmap[carname],time)
            vel_cor[frame][carname] =  make_ac(carname=carname,frame=frame,prefix=prefix,c=c,l=linemap[carname],use='ytvel')
#
#    aclist.append( make_ac(carname='p52_441',frame=frame,use='ytvel',prefix=r'$D22\ \mathbf{v}\ t=%0.1fs$'%time,c='r',l='--'))
#    aclist.append( make_ac(carname='p52_432',frame=frame,use='ytvel',prefix=r'$D26\ \mathbf{v}\ t=%0.1fs$'%time,c='r',l=':'))
#    aclist.append( make_ac(carname='p52_433',frame=frame,use='ytvel',prefix=r'$D27\ \mathbf{v}\ t=%0.1fs$'%time,c='r',l='-.'))
#    aclist.append( make_ac(carname='p52_434',frame=frame,use='ytvel',prefix=r'$D28\ \mathbf{v}\ t=%0.1fs$'%time,c='r',l='-'))

if 'turb' not in aclist[0]:
    if 'turb' not in aclist[0] or clobber:
        aclist[0]['turb']=make_ac(carname="/data/astrofs_2/davidc/P58_synchrotron/mshalf_ma2/cube_256_mshalf_ma2_0011_",
                          use='fits',prefix=r'$turb$',frame=11,do_mean=True,c='b',l='-')

#aclist.append(make_ac(carname="/data/astrofs_2/davidc/P58_synchrotron/mshalf_ma2/cube_256_mshalf_ma2_0011_",
#                      use='fits',prefix='mshalf_ma2',frame=11))
#aclist.append(make_ac(carname="/data/astrofs_2/davidc/P58_synchrotron/ms1_ma2/cube_256_ms1_ma2_0011_",
#                      use='fits',prefix='ms1_ma2',frame=11))
#aclist.append(make_ac(carname="/data/astrofs_2/davidc/P58_synchrotron/ms1_mahalf/cube_256_ms1_mahalf_0011_",
#                      use='fits',prefix='ms1_mahalf',frame=11))
#aclist.append(make_ac(carname="/data/astrofs_2/davidc/P58_synchrotron/ms1_mahalf/cube_256_ms1_mahalf_0011_",
#                      use='fits',prefix='ms1_mahalf_rel',frame=11,do_mean=True))
#aclist.append(make_ac(carname=None, use='test',prefix='unit',frame=11))
#aclist[-1] =make_ac(carname=None, use='test',prefix='unit',frame=11)

#aclist.append( make_ac(carname='p52_432',frame=0,prefix='432'))
#aclist.append( make_ac(carname='p52_433',frame=0,prefix='433'))
#aclist.append( make_ac(carname='p52_434',frame=0,prefix='434'))
#aclist.append( make_ac(carname='p52_441',frame=0,prefix='441'))
#aclist.append( make_ac(carname='p52_432',frame=60, prefix='432'))
#aclist.append( make_ac(carname='p52_433',frame=60, prefix='433'))
#aclist.append( make_ac(carname='p52_434',frame=60, prefix='434'))
#aclist.append( make_ac(carname='p52_441',frame=60, prefix='441'))

#aclist.append( make_ac(carname='p52_441',frame=100,prefix='441-vel',use='ytvel'))
#aclist.append( make_ac(carname='p52_434',frame=100,prefix='434-vel',use='ytvel'))
x_units = 400000000
def compute_Lac(self):
    bins = np.fft.fftfreq(self.AC3dft.shape[0])
    ok = bins > 0
    bins = np.r_[0,bins[ok]]
    db = bins[1:]-bins[:-1]
    AC = self.binned_ac[1]
    Lac = np.sum(AC*db)/AC[0]
    Lac = Lac*x_units
    return Lac
def plot_ac_rect(self,ax,**kwargs):
    Lac = compute_Lac(self) #change later
    rect=patches.Rectangle((0,0),Lac,self.binned_ac[1][0],**kwargs)
    ax.add_patch(rect)


if 1:

    fig,ax = plt.subplots(1,1)
    frame_list = sorted(list(aclist.keys()))
    for frame in frame_list:
        sims = sorted(list(aclist[frame].keys()))
        for sim in sims:
            ac = aclist[frame][sim]
            c=colordict.get(frame,'g')
            if sim == 'turb':
                c='b'
            ax.plot( ac.binned_ac[0]*x_units, ac.binned_ac[1], label= '%s'%(ac.prefix),
                    c=c, linestyle=linemap.get(sim,'-'))

    acr = aclist[0]['p52_434']
    #plot_ac_rect(acr,ax,facecolor=[1,0,0,0.5])
    #plot_ac_rect(aclist[0]['turb'],ax,facecolor=[0,1,0,0.5])
    ax.set_xlabel(r'$\Delta r\ [cm]$')
    ax.set_ylabel(r'$AC_\theta$')
    prefix='4xx'
    ax.legend(loc=0)
    outname='p52n_%s_ac_theta_n%04d.pdf'%(prefix,frame)
    fig.savefig(outname)
    plt.close(fig)
    print(outname)

if 1:
    my_Lacs={}
    frame_list_all = frame_list.copy()
    frame_list_all.pop(0)
    fig,ax=plt.subplots(1,1)
    for sim in ['p52_441','p52_434']:
        my_Lacs[sim]=[]
        for frame in frame_list_all:
            ac = aclist[frame][sim]
            my_Lacs[sim].append( compute_Lac(ac))
        ax.plot(nar(frame_list_all)/100,nar(my_Lacs[sim])/2e8,label=labelmap[sim])
    axbonk(ax,xlabel=r'$t\ [\rm{s}]$',ylabel=r'$L_{\rm{AC}}/R_{\rm{WD}}$')
    ax.legend(loc=0)
    outname = 'p52_ac_time.pdf'
    fig.savefig(outname)
    print(outname)


if 1:
    my_Lacs_v={}
    frame_list_all = sorted(list(vel_cor.keys()))
    fig,ax=plt.subplots(1,1)
    for sim in ['p52_441','p52_434']:
        my_Lacs_v[sim]=[]
        for frame in frame_list_all:
            ac = vel_cor[frame][sim]
            my_Lacs_v[sim].append( compute_Lac(ac))
        ax.plot(frame_list_all,my_Lacs_v[sim],label=sim)
    axbonk(ax,xlabel=r'$t\ [\rm{s}]$',ylabel=r'$L_{\rm{AC}}$')
    ax.legend(loc=0)
    outname = 'p52_ac_v_time.pdf'
    fig.savefig(outname)
    print(outname)










if 0:
    fig3,a24=plt.subplots(1,1,figsize=figsize)
    if 1:
        if 'AC3d' not in dir():
            print('Correlate')
            AC3d=scipy.signal.correlate(rho1,rho2,mode='same',method='fft')
            AC3d=np.roll(AC3d, AC3d.shape[0]//2,axis=0)
            AC3d=np.roll(AC3d, AC3d.shape[1]//2,axis=1)
            AC3d=np.roll(AC3d, AC3d.shape[2]//2,axis=2)
            print('rolled')
                  

        print( 'fft')
        rhohat = np.fft.ifftn(rho1)
        rhohatdag = rhohat.conj()
        rho_prod = rhohat*rhohatdag
        print( 'fft2')
        AC3dft = np.fft.fftn(rho_prod)

        k_array = make_k_freqs( AC3dft.shape[0],real=False)
        kmag = np.sqrt(k_array[0,...]**2 + k_array[1,...]**2 + k_array[2,...]**2)
        ktt=kmag.flatten()
        R_sphere = 0.25


        a21.imshow( (AC3d.real).sum(axis=axis) , cmap='Greys')
        proj = AC3dft.real.sum(axis=axis)
        a23.imshow(proj,cmap='Greys')

    if 1:    

        bins = np.fft.fftfreq(AC3dft.shape[0])
        ok = bins > 0
        bins = np.r_[0,bins[ok]]
        db = bins[1:]-bins[:-1]
        print('binning')
        binned=rb.rb( k_array, AC3dft.real,bins=bins)
        binned2=rb.rb( k_array, AC3d.real/AC3d.size,bins=bins)

    if 1:
        #a24.plot( binned[0],binned[1],c='r',label='binned ACfft')
        a24.plot( binned2[0],binned2[1],c='g',label=r'$AC_\rho(r)$')

        ok = binned2[1]>0
        #fits = np.polyfit(np.log10(binned2[0][ok]), np.log10(binned2[1][ok]),2)
        #a24.plot( binned2[0][ok], binned2[0][ok]**fits[-1])

    if 0:
        """fit to power laws.  Clearly this is not a great thing to do."""
        popt, pcov = curve_fit(powerlaw, binned2[0][ok][:20],np.log10(binned2[1][ok][:20]), p0=[1,1,-2])
        fit_rho0, fit_r0, fit_alpha = popt

        popt2, pcov2 = curve_fit(powerlaw2, binned2[0][ok][:20],np.log10(binned2[1][ok][:20]), p0=[1,-2])
        fit_r02, fit_alpha2 = popt2

        rrr = binned2[0]
        acac = binned2[1]
        a24.plot( rrr, 10**powerlaw(rrr, fit_rho0, fit_r0, fit_alpha),label='powerlaw')
        a24.plot( rrr, 10**powerlaw2(rrr, fit_r02, fit_alpha2),label='tweak2')
        a24.plot( rrr, 10**powerlaw2(rrr, L, fit_alpha2),label='tweak2')

    if 1:
        """actual correlation length"""

        AC = binned2[1]
        L = np.sum(AC*db)/AC[0]

    if 1:
        """do we want the rectangle?"""
        rect=patches.Rectangle((0,0),L,AC[0],facecolor=[0.2]*3)
        a24.add_patch(rect)

    if 1:
        axbonk(a24,xlabel='$r$', ylabel=r'$\rm{AC}_\rho(r)$',
               #yscale='log',xscale='log')
               yscale='linear',xscale='linear')
        #a24.set_yscale('log')
        #x0a24.set_xscale('log')
        

    if 0:
        a24.legend(loc=0)

    fig2.savefig('p52_convolve_4.png')
    fig3.savefig('p52_convolve_5.png')
    print('saved')
