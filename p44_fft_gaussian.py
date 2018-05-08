#From https://stackoverflow.com/questions/5398304/fourier-transform-of-a-gaussian-is-not-a-gaussian-but-thats-wrong-python
import matplotlib.pyplot as plt
import numpy as np

if 1:
    #works
    N = 128
    #k = np.arange(-5,5,10./(2*N))
    k = np.arange(-5,5,10./(2*N))+5
    var=0.25
    k0=0
    fk = np.exp(-(k-k0)**2/(2*var))
    fx=np.fft.ifft(fk)
    fx_shift = np.fft.fftshift(np.abs(fx))* np.sqrt(2 * N)

    print "From Plancherel", np.sum(np.abs(fk)**2)/(2*N) /np.sum(np.abs(fx)**2)

    plt.clf()
    plt.plot(k,fk,label='fk',c='r')
    plt.plot(k,fx_shift,label='fx_shift',c='g')
    plt.plot(k,fx.real,label='fx.real',c='b')
    plt.legend(loc=0)
    plt.savefig('p44_gaussian_test.pdf')

if 0:
    N = 256
    DeltaX = 10.
    dx = DeltaX/N
    x = np.arange(N)*dx
    Nk = N/2+1
    dk = 1/(2.*N*dx)
    k = np.arange(Nk)*dk
    var=0.25
    k0=0
    fk = np.exp(-(k-k0)**2/(2*var))
#fk = np.zeros(2*Nk)

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

    fx=np.fft.irfft(fk)
    plt.clf()
    plt.close('all')
    fig, axes = plt.subplots(1,2)
    ax_k = axes[0]
    ax_x = axes[1]
    ax_k.plot(k,fk, c='r',label='fk')

    fx_shifted = np.fft.fftshift(fx)*np.sqrt(N)

    var_fk = np.sum(np.abs(fk)**2)
    var_fx = np.sum(np.abs(fx)**2)
    print "From Plancherel (2)", np.sum(np.abs(fk)**2)/(2*N) /np.sum(np.abs(fx)**2)

    ax_x.plot(x,fx, c='g',label='fx')
    ax_x.plot(x,fx_shifted, c='b',label='fx_shifted')
    ax_k.legend(loc=0)
    ax_x.legend(loc=0)
    fig.savefig('p44_gaussian_test_2.pdf')


    if 0:
        fx_shift = np.fft.fftshift(np.abs(fx))* np.sqrt(2 * N)
        var_measure = np.sum(fk*(k-k0)**2)/np.sum(fk)
        print("var %f var_measure %f"%(var,var_measure))

        plt.clf()
        plt.plot(k,fk,label='fk',c='r')
        plt.plot(k,fx_shift,label='fx_shift',c='g')
        plt.plot(k,fx.real,label='fx.real',c='b')
        plt.legend(loc=0)
        plt.savefig('p44_gaussian_test_2.pdf')
