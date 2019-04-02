#https://flothesof.github.io/Fourier-series-rectangle.html
if '/home/dcollins/repos/p49c/p49_eigenmodes' not in sys.path:
    sys.path.append('/home/dcollins/repos/p49c/p49_eigenmodes')
if './p49c' not in sys.path:
    sys.path.append('./p49c')
from p49_print_tools import *

from go import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scatter_fit
from scipy.special import gamma
size =10
dx=1./size
n_partial_sine = size//2
x = np.arange(0,1,dx)
mask=x>0
tp=np.pi*2
k=np.fft.fftfreq(size)*size
sft=np.fft.fftshift
pi=np.pi
plt.close('all')
plt.clf()
output_dir = "p60_dft_square/plots/"
sin_n = lambda n: np.sin(2 * np.pi * n * x)
F_n = lambda n: (1/2. + np.sum([2./np.pi/(2*i+1) * np.sin(2 * pi * (2*i+1) * x) for i in range(n+1)], axis=0))

sample_freq = x.size / 1.
my_k=np.linspace(0, 1, size) * sample_freq

if 0:
    plt.plot(x, x<1/2.)
    for i in [0, 1, 2, 3, 5, 10]:
        plt.plot(x, F_n(i))
    plt.ylim(-1.1, 2.1)
    plt.savefig("%s/square_direct_1.png"%output_dir)

if 0:
    plt.clf()
    plt.plot(x, x<1/2.)
    plt.plot(x, F_n(n_partial_sine))
    plt.ylim(-1.1, 2.1)
    plt.title('square wave approximation with %i terms' % n_partial_sine)
    plt.savefig('%s/square_direct_2.png'%output_dir)

if 0:
    plt.clf()
    rect_fft = np.fft.fft(x<1/2.)
    plt.plot(my_k, np.abs(rect_fft) ,c='k')
    plt.plot(my_k, np.real(rect_fft),c='r')
    plt.plot(my_k, np.imag(rect_fft),c='g')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('amplitude')
    plt.xlim(0, sample_freq / 2)
    plt.savefig('%s/square_direct_3.png'%output_dir)

if 0: 
    #examine the DFT of the sine series
    plt.clf()
    fhat=np.fft.fft(F_n(x.size//2))
    my_k=np.linspace(0, 1, fhat.size) * sample_freq
    plt.plot(my_k,  np.abs(fhat) ,c=[0.5]*3)
    plt.plot(my_k, np.real(fhat),c='m')
    plt.plot(my_k, np.imag(fhat),c='y')

    plt.plot(my_k, np.abs(rect_fft) ,c='k')
    plt.plot(my_k, np.real(rect_fft),c='r')
    plt.plot(my_k, np.imag(rect_fft),c='g')
    plt.xlim(0, sample_freq / 2)

    plt.savefig('%s/square_direct_4.png'%output_dir)

if 0:
    #plot the deviations from the actual DFT and the sine series
    fig,ax=plt.subplots(1,2)
    ax0=ax[0]
    ax1=ax[1]
    ax0.plot(rect_fft.imag-fhat.imag,c='g')
    ax0.plot(rect_fft.real-fhat.real,c='r')
    ax1.plot(np.abs(rect_fft)-np.abs(fhat),c='k')
    di=rect_fft.imag-fhat.imag
    dr=rect_fft.real-fhat.real
    da=np.abs(rect_fft)-np.abs(fhat)
    fig.savefig('%s/square_direct_5.png'%output_dir)

if 0:
    #Look at the amplitudes of the sine series
    #My fourier amplitudes for fhat should be 2/(pi k)*n
    #Which they are.
    plt.clf()
    f2=F_n(1) #1/2.+2/np.pi*np.sin(2*pi*x)
    f2hat = np.fft.fft(f2)
    plt.plot(x,f2)
    plt.plot(x,F_n(0))
    plt.savefig('%s/square_direct_6.png'%output_dir)

if 0:
    #It's really a string of delta functions.
    #What do delta functions do?
    f2=np.zeros_like(x)
    delta_list = range(5)
    for d in delta_list:
        f2[d]=1
    f2hat = np.fft.fft(f2)

    fig,ax=plt.subplots(1,2)
    ax0=ax[0];ax1=ax[1]
    #ax0.plot(x,f2)
    ax0.set_title('fourier')
    ax0.plot(x,f2hat.real,c='r')
    ax0.plot(x,f2hat.imag,c='g')
    expect_real=sum([np.cos(-2*np.pi*d*x) for d in delta_list])
    expect_imag=sum([np.sin(-2*np.pi*d*x) for d in delta_list])
    dr=f2hat.real-expect_real
    di=f2hat.imag-expect_imag
    ax1.plot(x,dr,c='r')
    ax1.plot(x,di,c='g')
    ax1.set_title('difference with sums of sines.')
    fig.savefig('%s/square_direct_6.png'%output_dir)

if 1:
    #Why is the real part non-zero?  
    #It seems to be one on the odds, zero on evens, 
    #This comes from the lack of end point.
    # sum( cos(x), {x,0,pi}) = 0,
    # but if the last end point is excluded, = 1.
    
    delta_list = range(5)
    rm=rainbow_map(len(delta_list))

    x2=np.arange(0,1,1./100)
    x=np.arange(0,1,1./size)
    k=x+0  #k is also integer valued.
    kint = x*size
    expect_real_full=sum([np.cos(-2*np.pi*m*k) for m in delta_list])
    expect_imag_full=sum([np.sin(-2*np.pi*m*k) for m in delta_list])

    plt.clf()
    partials=[]
    test1=0j
    args=[]
    kmy=kint
    kmy = np.exp(-2*np.pi*1j*kint/size)
    #othertest = (1-kmy**(m+1))/(1-kmy)
    for m in delta_list:#range(1,50):#range(10,60,10):
        args.append(-2*np.pi*m/size)
        partials.append(np.cos(args[-1]*kint))
        test1 = kmy**m + test1
        expect_real_part=sum(partials)#range(nd)])
        #plt.plot(kint,expect_real_part,marker='*',c=rm(m))
        #plt.plot(kint,partials[-1],marker='*',c=rm(m),linestyle=':',label='m=%d'%m)
        #plt.savefig('%s/square_direct_7.png'%output_dir)
        #time.sleep(1)

    
    partials = nar(partials)
    othertest = (1-kmy**(m+1))/(1-kmy)
    print("other", othertest)
    print("test",test1)
    print("ratio",othertest/test1)
    args=nar(args)
    n=delta_list[-1]
    plt.plot(kint,expect_real_full,c='k')

    plt.legend(loc=0)

    plt.xlabel('k')
    plt.savefig('%s/square_direct_7.png'%output_dir)

if 0:
    #Why is the imaginary part not quite right?
    
    delta_list = range(5)
    rm=rainbow_map(len(delta_list))

    x2=np.arange(0,1,1./100)
    x=np.arange(0,1,1./size)
    k=x+0  #k is also integer valued.
    kint = x*size
    Ny = size//2
    expect_real_full=sum([np.cos(-2*np.pi*m*k) for m in delta_list])
    expect_imag_full=sum([np.sin(-2*np.pi*m*k) for m in delta_list])

    plt.clf()
    partials=[]
    args=[]
    for m in delta_list:#range(1,50):#range(10,60,10):
        args.append(-2*np.pi*m/size)
        delta = Ny-m
        fakery = -np.cos(-tp*Ny*kint/size)*np.sin(-tp*delta*kint/size)
        #partials.append(np.sin(args[-1]*kint))
        partials.append(fakery)
        expect_imag_part=sum(partials)#range(nd)])
        plt.plot(kint,expect_imag_part,marker='*',c=rm(m))
        plt.plot(kint,partials[-1],marker='*',c=rm(m),linestyle=':',label='m=%d'%m)
        #plt.savefig('%s/square_direct_7.png'%output_dir)
        #time.sleep(1)
    partials = nar(partials)
    args=nar(args)



if 0:
    #can i sort it out using nyquist?
    x2=np.arange(0,1,1./100)
    x=np.arange(0,1,1./size)
    k=x+0  #k is also integer valued.
    kint = x*size
    Ny = size//2
    def integral0(kint):
        dx=1./size
        kny = 2./dx
        dk = kny-kint
        a = -2*np.pi*(dk)/(kny**2-dk**2)
        return a
    def integral1(kint,xi):
        dx=1./size
        kny = 2./dx
        dk = kny-kint
        a = -2*np.pi*(dk)/(kny**2-dk**2)
        a *= -np.cos(-tp*Ny*xi/size)*np.cos(-tp*dk*xi/size)
        return a
    plt.clf()
    plt.plot(kint,f2hat.imag,c='k')
    plt.plot(kint[1:],take1,c='k',linestyle=':')
    take1=-size/np.pi/kint[1:]
    take1ish = take1[::2]
    take2 = integral1(kint, size/2)-integral1(kint,0)
    #plt.plot(kint[1:],take1*np.cos(2*np.pi*(kny-kint[1:])/size),c='g')
    plt.savefig('%s/square_direct_8.png'%output_dir)
    right=f2hat.imag[1::2]

if 0:
    #messing around.
    plt.plot(kint,f2hat.imag,c='k')
    kny = size//2
    take1=-size/np.pi/kint[1:]
    plt.plot(kint[1:],take1,c='k',linestyle=':')
    plt.plot(kint[1:],take1*np.cos(2*np.pi*(kny-kint[1:])/size),c='g')
    plt.legend(loc=0)

    plt.xlabel('k')
    plt.savefig('%s/square_direct_7.png'%output_dir)
    plt.clf()
    expect_full = 1j*take1+1
    plt.plot(kint,np.abs(f2hat),c='k')
    plt.plot(kint,np.abs(f2hat.real),'k--')
    plt.plot(kint,np.abs(f2hat.imag),'k:')
    plt.plot(kint[1:],np.abs(take1),c='g')
    plt.plot(kint[1:],np.abs(expect_full),c='c')
    #plt.plot(kint[1:],np.abs(f2hat[1:])/np.abs(expect_full),c='m')
    plt.plot(kint[1:],np.angle(f2hat[1:])/-np.pi,c='m')
    plt.savefig('%s/square_direct_9.png'%output_dir)

if 0:
    #Using Lagrange's trig identities
    def many_sin(N,theta):
        #sum_n=1^N sin(n T) = 1/2 cot(T/2)-cos( (N+1/2)T))/(2 sin(T/2))
        out = 0.5/np.tan(theta/2)
        out += np.cos( (N+1./2)*theta)/( 2*np.sin(theta/2.))
        return out
    del size
    size1 = 8
    x=np.arange(0.,1.,1./size1)
    kint = x*size1
    x2=np.arange(0.,1.,1./size2)
    args1=[]
    partials_i=[]
    partials_r=[]
    nk = 3 #int(f.sum())
    fig9, ax9 = plt.subplots(1,1)
    rm = rainbow_map(nk)
    fig,ax = plt.subplots(1,1)
    for m in range(nk):
        args1.append(-2*np.pi*m/size1)
        partials_i.append(np.sin(args1[-1]*kint))
        partials_r.append(np.cos(args1[-1]*kint))
        #ax.plot(kint,partials_i[-1],c=rm(m),marker='*')
    

    ok = slice(1,None)
    tots = sum(partials_i)[ok]
    dur = many_sin(nk,-2*np.pi*kint[ok]/size1)
    ax.plot(tots,dur)
    ax.set_xlabel('sum'); ax.set_ylabel('lagrange')
    #ax.plot(kint,dur,c='k')
    fig.savefig('%s/direct_lagrange.png'%output_dir)
    plt.close(fig)

if 0:
    #more messing around.
    #bigger sizes.
    del size
    size1 = 8
    size2 = 16#//2+1#2*size1
    x=np.arange(0.,1.,1./size1)
    x2=np.arange(0.,1.,1./size2)
    f=(x<0.5).astype('float')
    fhat = np.fft.fft(f)
    fhat2 = np.zeros(size2,dtype='complex')
    fhat2[:size1]=fhat
    fhat2[size2:size2//2:-1] = fhat2[1:size2//2].conj()

    f2 = np.fft.ifft(fhat2)
    k1 = np.fft.fftfreq(size1,d=1/size1)
    k2 = np.fft.fftfreq(size2,d=1/size2)
    k1=k1#[k1>=0]
    k2=k2#[k2>=0]

    fig, ax = plt.subplots(2,2)
    ax0=ax[0][0]
    ax1=ax[0][1]
    ax2=ax[1][0]
    ax3=ax[1][1]

    ax0.plot(f2)
    ax1.plot(np.fft.ifft(fhat))


    args1=[]
    args2=[]
    partials1=[]
    partials2=[]
    nk = int(f.sum())
    fig9, ax9 = plt.subplots(1,1)
    rm = rainbow_map(nk)
    for m in range(nk):
        args1.append(-2*np.pi*m/size1)
        args2.append(-2*np.pi*m/size2)
        #delta = Ny-m
        #fakery = -np.cos(-tp*Ny*kint/size)*np.sin(-tp*delta*kint/size)
        #partials.append(fakery)
        partials1.append(np.exp(1j*args1[-1]*k1))
        partials2.append(np.exp(1j*args2[-1]*k2))
        #ax2.plot(partials2[-1],c=rm(m))
        #ax3.plot(partials1[-1],c=rm(m))
        ax9.plot(x2,partials2[-1],c=rm(m),linestyle=':',marker='*')
        ax9.plot(x,partials1[-1],c=rm(m),marker='*')
#       expect_imag_part=sum(partials)#range(nd)])

    ax3.plot(sum(partials1))
    ax2.plot(sum(partials2))
    fig.savefig('%s/direct_10.png'%output_dir)
    fig9.savefig('%s/direct_11.png'%output_dir)
    plt.close('all')
