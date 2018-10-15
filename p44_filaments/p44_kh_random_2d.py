if 'ef' not in dir():
    execfile('go')
plt.close('all')

import fourier_tools.fourier_filter as Filter
def sym2(b2):
    Nx=b2.shape[0]
    b3 = np.zeros_like(b2)
    b3[1:Nx/2,:] = b2[1:Nx/2,:]
    b3[Nx:Nx/2:-1,:]=b2[1:Nx/2,:].conj()
    return b3
class rtool():
    def __init__(self,field, dk=None, bin_size=None):
        self.field=field
        self.radial_average(dk,bin_size)
    
    def radial_average(self, dk=None, bin_size=None):
        self._kk = []
        field=self.field
        ones = np.ones_like(field)
        l = len(field.shape)
        if dk is None:
            dk = np.ones_like(field.shape)
        for i,dim in enumerate(field.shape):
            sl = i*(1,)+(dim,)+(l-i-1)*(1,)
            #k = np.fft.fftfreq(dim)
            k=np.arange(0,dim*dk[i], dk[i])
            k.resize(sl)
            self._kk.append(k)
        self.kradius = na.sqrt(na.sum((k**2 for k in self._kk)))
        self.max_radius=self.kradius.max()
        if bin_size is None:
            bin_size=min(dk)
        self.bin_centers = np.arange(0, self.max_radius,bin_size)
        self.zones_per_bin = np.zeros_like(self.bin_centers)
        self.average = np.zeros_like(self.bin_centers)
        for ib,b in enumerate(self.bin_centers):
            L = b-0.5*bin_size
            R = b+0.5*bin_size
            select = ((self.kradius >= L) & (self.kradius <= R))
            self.zones_per_bin[ib] = ones[select].sum()
            self.average[ib] = field[select].sum()/self.zones_per_bin[ib]
        self.averagem = self.average[:min(field.shape)+1]
    



Delta = 1 #0.15 # the physical pixel size in units of length
Delta_xi = Delta; Delta_xj = Delta
N = 200 # the number of pixels


Nharm = N/2 + 1 # number of harmonics for real fft
Deltak = 1.0# 2*np.pi/N/Delta # size of pixels in harmonic space
Delta_ki = Deltak; Delta_kj = Deltak

Nxi = N; Nxj = N
Nki = Nxi; Nkj = Nxj/2+1
x = np.arange(N)*Delta # the x-values
k = np.arange(Nharm)*Deltak # the k values
x_j, x_i = np.mgrid[0:Nxi*Delta_xi:Delta_xi, 0:Nxj*Delta_xj:Delta_xj]
x = (x_j**2+x_i**2)**0.5
#ki_mg = Nki-1 #for the endpoints of the mgrid
#kj_mg = Nkj-1 #for the endpoints of the mgrid
#k_j, k_i = np.mgrid[0:ki_mg*Delta_ki:Delta_ki, 0:kj_mg*Delta_kj:Delta_kj]
k_j, k_i = np.mgrid[0:Nki, 0:Nkj]*Deltak
k2d = (k_j**2+k_i**2)**0.5
k2d = sym2(k2d)
        
#rt = rtool(k2d, dk = [Delta_ki, Delta_kj])
###
# These parameters specify the lognormal field rho: the mean and the power spectrum
###
if 1:
    rhobar = 0.25 # mean of lognormal field
    sigma_k=0.333
    Prhok = np.exp(-k**2/2./sigma_k**2) # power spectrum of lognormal field
    Prhok2d = np.exp(-k2d**2/2./sigma_k**2) # power spectrum of lognormal field
    Prhok2d=sym2(Prhok2d)
    #plt.clf()
    #dumb_plt(plt,rt.average, Prt.average,'k','Prho','p44_kh_Prhok2d.png',scale=('linear','log'),c='r')
    #dumb_plt(plt,k, Prhok,'k','Prho','p44_kh_Prhok2d.png',scale=('linear','log'),c='k')
    #plt.clf()
    #plt.plot(Prhok,c='r')
    #plt.plot(Prhok2d[k_j==k_i],c='g')
    #plt.savefig('tmp.png')



###
# Set up ancillary statistics, overdensity delta and Gaussianized field A
###


    Pdeltak = Prhok/rhobar**2 # power spectrum of overdensity field
    Pdeltak2d=Prhok2d/rhobar**2

    Xidelta = N*np.fft.irfft(Pdeltak)#*Deltak/2/np.pi  # iFT to get correlation function of overdensity, careful of numpy 1/N normalization on iFT (differs from eg FFTW).
    Xidelta_2d = N*np.fft.irfft2(Pdeltak2d)#*(Deltak/2/np.pi)  # iFT to get correlation function of overdensity, careful of numpy 1/N normalization on iFT (differs from eg FFTW).
    vx = np.var(Xidelta)
    vp = np.var(Pdeltak)
    print("var xi", vx)
    print("var Pk", vp)
    print("var xi/Pk", vp/vx)
    vx2 = np.var(Xidelta_2d)
    vp2 = np.var(Pdeltak2d)
    print("var xi2", vx2)
    print("var Pk2", vp2)
    print("var xi/Pk", vp2/vx2)
    


    XiA = np.log(1 + Xidelta) # correlation function of Gaussianized field
    XiA_2d = np.log(1 + Xidelta_2d) # correlation function of Gaussianized field

    sig2A = XiA[0] # variance of Gaussianized field
    sig2A_2d = XiA_2d[0,0] # variance of Gaussianized field

    Abar = np.log(rhobar) - 0.5*sig2A # mean of Gaussianized field
    Abar_2d = np.log(rhobar) - 0.5*sig2A_2d # mean of Gaussianized field

    PAk =  np.fft.rfft(XiA)#*Delta # power spectrum of gaussianized field
    PAk_2d =  np.fft.rfft2(XiA_2d)#*Delta # power spectrum of gaussianized field



###
# Make a realization (no changes to 2d until here)
###
#TODO: Nharm in 2d
#     
# create some unit variance complex gaussian deviates
# xi, xj are real space X and Y
    np.random.seed(123123)
    uharm = (np.random.randn(Nharm) + 1.0j*np.random.randn(Nharm) )/np.sqrt(2)
    uharm2d = (np.random.randn(Nki, Nkj) + 1.0j*np.random.randn(Nki, Nkj) )/np.sqrt(2) 
#uharm[0] = 0.0  # zero the DC component

    norm1 = 1# 2*np.pi/Deltak
    PAkAmp = np.sqrt(norm1* PAk)  #make this 2d, with the right symmetry; 2pi/Deltak needs (2pi)^2/(dki dkj)
    norm2 = 1 # 2*np.pi/Deltak
    PAkAmp2d = np.sqrt(norm2* PAk_2d)  #make this 2d, with the right symmetry; 2pi/Deltak needs (2pi)^2/(dki dkj)

    Aharm = uharm * PAkAmp # Compute the harmonics of the Gaussianized field
    A = Abar + Deltak/2/np.pi * N*np.fft.irfft(Aharm) # Compute real-space Gaussianized field

    Prhok2d = np.exp(-k2d**2/2./sigma_k**2) # power spectrum of lognormal field
    PAkAmp2d =sym2(Prhok2d)
    Aharm_2d = uharm2d * PAkAmp2d # Compute the harmonics of the Gaussianized field
    norm3 = Deltak/2/np.pi
    A_2d = Abar_2d + norm3 * N*np.fft.irfft(Aharm_2d) # Compute real-space Gaussianized field
    plt.clf()
    plt.hist(A_2d.flatten())
    plt.savefig('tmp.png')

    rho = np.exp(A) # make lognormal field
    rho_2d = np.exp(A_2d) # make lognormal field

    rhoharm = np.fft.rfft(rho)*Delta  # compute the harmonics of the lognormal field to check power spectrum
    rhoharm_2d = np.fft.rfft(rho_2d)*Delta  # compute the harmonics of the lognormal field to check power spectrum
    rho_rt = rtool( (rhoharm_2d*rhoharm_2d.conj()).real )
if 0:
    plt.clf()
    plt.imshow(rho_2d)
    plt.savefig('tmp_2d.png')
    plt.clf()
    plt.plot(rho_rt.average)
    plt.yscale('log')
    plt.savefig('tmp.png')


if 0:
# some plots
    #fig = figure()
    plt.clf()
    plt.plot(k,Prhok)
    plt.plot(k,abs(rhoharm)**2*Deltak/2/np.pi)
    plt.xscale('log');plt.yscale('log')
    plt.title('Prhok input vs output')
    plt.savefig('p44_kh_k_P.png')

    #figure()
    plt.clf()
    plt.plot(x,Xidelta,label='Xidelta')
    plt.plot(x,XiA,label='XiA')
    plt.legend()
    plt.savefig('p44_kh_Xi.png')

    #figure()
    plt.clf()
    plt.plot(k,PAk)
    plt.plot(k,abs(Aharm)**2*Deltak/2/np.pi)
    plt.title('PAk model vs output')
    plt.xscale('log');plt.yscale('log')
    plt.savefig('p44_kh_model.png')

    plt.clf()
    plt.plot(x,A)
    plt.title('A')
    plt.savefig('p44_kh_A.png')

    #figure()
    plt.clf()
    plt.plot(x,rho)
    plt.title('rho')
    plt.savefig('p44_kh_rho.png')

    #figure()
    plt.clf()
    plt.hist(rho,bins=100)
    plt.title('hist rho')
    plt.savefig('p44_kh_hist_rho.png')

    #figure()
    plt.clf()
    plt.hist(A,bins=50)
    plt.title('hist A')
    plt.savefig('p44_kh_hist_A.png')

