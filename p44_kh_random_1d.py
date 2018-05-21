if 'ef' not in dir():
    execfile('go')
plt.close('all')

Delta = 1.0 #0.15 # the physical pixel size in units of length
N = 8192 # the number of pixels
#N = 100 # the number of pixels


Nharm = N/2 + 1 # number of harmonics for real fft
Deltak = 2*np.pi/N/Delta # size of pixels in harmonic space

x = np.arange(N)*Delta # the x-values
k = np.arange(Nharm)*Deltak # the k values


###
# These parameters specify the lognormal field rho: the mean and the power spectrum
###

rhobar = 0.25 # mean of lognormal field
sigma_rho=0.5
Prhok = np.exp(-k**2/2./sigma_rho**2) # power spectrum of lognormal field


###
# Set up ancillary statistics, overdensity delta and Gaussianized field A
###


Pdeltak = Prhok/rhobar**2 # power spectrum of overdensity field

Xidelta = N*np.fft.irfft(Pdeltak)*Deltak/2/np.pi  # iFT to get correlation function of overdensity, careful of numpy 1/N normalization on iFT (differs from eg FFTW).

XiA = np.log(1 + Xidelta) # correlation function of Gaussianized field

sig2A = XiA[0] # variance of Gaussianized field

Abar = np.log(rhobar) - 0.5*sig2A # mean of Gaussianized field

PAk =  np.fft.rfft(XiA)*Delta # power spectrum of gaussianized field



###
# Make a realization
###

# create some unit variance complex gaussian deviates
np.random.seed(123123)
uharm = (np.random.randn(Nharm) + 1.0j*np.random.randn(Nharm) )/np.sqrt(2)
#uharm[0] = 0.0  # zero the DC component


Aharm = uharm * np.sqrt(2*np.pi/Deltak* PAk)  # Compute the harmonics of the Gaussianized field
A = Abar + Deltak/2/np.pi * N*np.fft.irfft(Aharm) # Compute real-space Gaussianized field

rho = np.exp(A) # make lognormal field

rhoharm = np.fft.rfft(rho)*Delta  # compute the harmonics of the lognormal field to check power spectrum
if 1:
# some plots
    #fig = figure()
    plt.clf()
    plt.plot(k,Prhok, c='r')
    plt.plot(k,abs(rhoharm)**2*Deltak/2/np.pi, c='g')
    plt.xscale('linear');plt.yscale('log')
    plt.title('Prhok input vs output')
    plt.savefig('p44_kh_k_P.png')

    #figure()
    plt.clf()
    #plt.plot(x,Xidelta,label='Xidelta')
    plt.plot(x,Xidelta,label='Xidelta')
    #plt.plot(x,XiA,label='XiA')
    plt.yscale('log')
    plt.legend()
    plt.savefig('p44_kh_Xi.png')

    #figure()
    plt.clf()
    plt.plot(k,PAk,c='r')
    plt.plot(k,abs(Aharm)**2*Deltak/2/np.pi,c='g')
    plt.title('PAk model vs output')
    plt.xscale('linear');plt.yscale('log')
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

