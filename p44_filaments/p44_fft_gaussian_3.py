#From https://stackoverflow.com/questions/5398304/fourier-transform-of-a-gaussian-is-not-a-gaussian-but-thats-wrong-python
import matplotlib.pyplot as plt
import numpy as np

def symmetric(v):
    #for real fft, Ak = A(-k)^*; the negative-phase is the conjugate
    #(that way when you sum, the imaginary parts cancel.)
    #This lets us take an arbitrary K-space signal and ensure it's inverse-fft is real.
    s2=np.zeros_like(v)
    Nx=v.size
    s2[1:Nx/2] = v[1:Nx/2]
    s2[Nx:Nx/2:-1] = v[1:Nx/2].conj()
    s2[0]=v[0]
    return s2
test=0
if 1:
    #this works.
    N = 512
    kstart=-5.
    kend=5.
    var=0.333
    k0=0
    fk_norm_target=np.sqrt(2*np.pi*var)
    test=0
if 0:
    Nx = 512
    N = 512 #Nx/2+1
    dk = 10./512
    kstart=0
    kend=N*dk
    var=0.333
    k0=0
    k = np.linspace(0,kend,N)
    fk_norm_target=np.sqrt(2*np.pi*var)/2 #half a gaussian.
    test=1


dk = (kend-kstart)/(N)
k = np.arange(kstart,kend,dk)
k = np.fft.fftfreq(N,dk)
fk = np.exp(-(k-k0)**2/(2*var))
#fk = symmetric(fk)
fx=np.fft.ifft(fk)
fx_shift = np.fft.fftshift(np.abs(fx))* np.sqrt(2 * N)

print "From Plancherel", np.sum(np.abs(fk)**2)/(2*N) /np.sum(np.abs(fx)**2)
fk_norm_actual= fk.sum()*dk
print "norms: fk",fk_norm_target, fk_norm_actual, (fk_norm_target-fk_norm_actual)/fk_norm_target 
print "norms: ||fx||", np.sum(np.abs(fx))

plt.clf()
plt.plot(k,fk,label='fk',c='r')
if test in [0]:
    plt.plot(k,fx_shift,label='fx_shift',c='g')
    plt.plot(k,fx.real,label='fx.real',c='b')
if test in [1]:
    #x2 = np.linspace(0,1/(2.*dk),N/2)
    plt.plot(fx/fx.sum())
plt.legend(loc=0)
plt.savefig('p44_gaussian_test.pdf')

if test in [0]:
#what I think I should get
    x = k
#xstart = 1./kstart
#xend = 1./kend
#dx = 1./(N*dk)
#x = np.arange(xstart,xend,dx)
#should_x = (2*np.pi*var)**0.5/N*np.exp(-(np.pi)**2*2*var/(N**2)*(x-x0)**2)
#norm = (2*np.pi*var)**0.5/N #what I think it should be (not right)
#norm = fx_shift.max() #kludge
#newvar = (np.pi)**2*2*var/(N**2)  #from math that I got wrong
##newvar = 1/0.247                  #from the variance of the other one
#sigma_fx_shift = np.sum(fx_shift*k**2)/np.sum(fx_shift)
#newvar = 1/(2.*sigma_fx_shift)
#print sigma_fx_shift
    norm=np.sqrt(2*N)
    newvar = 1/(2*N*var)
    should_x_unnorm = np.exp(-(x-x0)**2/(2*newvar))
    should_x = should_x_unnorm * np.sqrt(2*N)/should_x_unnorm.sum()
    var_fx = np.sum(fx_shift*x**2)/np.sum(fx_shift)
    var_sh = np.sum(should_x*x**2)/np.sum(should_x)
    print "var_fx", var_fx
    print "var_sh", var_sh
    plt.clf()
    plt.plot(k,fx.real,label='fx.real',c='b')
    plt.plot(k,fx_shift,label='fx_shift',c='g')
    plt.plot(x,should_x,label='should',c='k')
    plt.yscale('log')
    plt.ylim(1e-16,100)
    plt.savefig('p44_gaussian_test_2.pdf')
