

import fourier_tools.fourier_filter as Filter
from numpy.fft import fftn, ifftn, fftfreq

def thing1(field):
    FTdivV = fftn(divV)
    N = FTdivV.shape[1]
    # assumes 3d periodic domain with L = 1
    kx = fftfreq(N)
    # is actually k/N
    kx, ky, kz = np.meshgrid(kx,kx,kx)
    # discrete fourier representation of -div grad based on consecutive 2nd order first derivatives 
    denom = -1/2. * N**2. * (np.cos(4.*np.pi*kx) + np.cos(4.*np.pi*ky) + np.cos(4.*np.pi*kz) - 3.)
    # these are 0 in the nominator anyway, so set this to 1 to avoid division by zero
    denom[denom == 0.] = 1.    

    FTdivV /= denom

def power1(field):
    fhat = fftn(field)/field.size
    power = fhat*np.conj(fhat)
    ffd = Filter.FourierFilter(power)
    kspace=ffd.get_shell_k()
    mk = 1.0 #kspace[ kspace != 0].min()
    power_1d = np.array([power[ffd.get_shell(bin)].sum() for bin in range(ffd.nx)])
    plt.clf()
    plt.plot((kspace/mk)[1:], power_1d[1:],c='r')
    #plt.xscale('log') 
    plt.xlabel(r'$k/k_{\rm{min}}$')
    #plt.yscale('log') 
    plt.ylabel(r'Power')
    return kspace, power

def power_all(data):
    #data = np.random.rand(301) - 0.5
    ft = np.fft.fft(data)
    ps = np.abs(ft)**2

    time_step = 1./data.size #1. / 30
    freqs = np.fft.fftfreq(data.size, time_step)
    idx = np.argsort(freqs)

    plt.plot(freqs[idx], ps[idx],c='k')
    return freqs[idx], ps[idx], ft

def power3(data):
    import pylab as pl
    rate = 30.0
    t = np.arange(0, 10, 1/rate)
    x = np.sin(2*np.pi*4*t) + np.sin(2*np.pi*7*t) + np.random.randn(len(t))*0.2
    p = 20*np.log10(np.abs(np.fft.rfft(x)))
    f = np.linspace(0, rate/2, len(p))
    plt.plot(f, p)
    return f,p

def power_half(data):
    #data = np.random.rand(301) - 0.5
    ft = np.fft.fft(data)/np.sqrt(data.shape)
    N=data.size
    ps = np.abs(ft)**2
    time_step = 1./data.size #1. / 30
    freqs = np.fft.fftfreq(data.size, time_step)
    idx = np.argsort(freqs)
    k=freqs[0:N/2]
    f=ps[0:N/2]
    return k,f,ft


plt.clf()
rate = 1./30
x = np.arange(0,1,rate) #np.arange(0,2*np.pi,2*np.pi/300)
signal1 = np.zeros_like(x)
modes = np.zeros(10)
modes[4] = 1.0
modes[7] = 2.0
k_used=np.arange(len(modes))
for k,ak in zip(k_used,modes):
    signal1 += ak*np.sin(2*np.pi*k*x)

out=power_half(signal1)
f=out[0]
p=out[1]
ft=out[2]
var_s = (signal1**2).sum()
var_f = (np.abs(ft)**2).sum()
print "squares signal %0.2e ft %0.2e relerr %0.2e"%( var_s, var_f, 1-var_f/var_s)

#plt.plot(k_used,modes,c='g')
#plt.plot(k,f,c='r')
plt.plot(f,p,c='b')
plt.savefig('test3.pdf')
