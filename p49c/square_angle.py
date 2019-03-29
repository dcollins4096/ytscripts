
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
size =30
dx=1./size
Ny = size//2
n_partial_sine = size//2
x = np.arange(0,1,dx)
kint = x*size
mask=x>0
tp=np.pi*2
tpi = np.pi*2j
kfreq=np.fft.fftfreq(size)*size
sft=np.fft.fftshift
pi=np.pi
plt.close('all')
plt.clf()


if 0:
    f2 = (x<0.5).astype('float')
    f2hat = np.fft.fft(f2)
    theta = np.angle(f2hat[1::2])
    print(Ny-theta/-np.pi*size)
    print(np.tan(theta)*f2hat.real[1::2]/f2hat.imag[1::2])
    print(f2hat.imag[1::2] - np.tan(theta))
    print(f2hat.imag[1::2] / np.tan(-(Ny-kint[1::2])/size*np.pi  ))
    test_real = np.zeros_like(x)
    test_real[1::2]=1
    test_real[0] = size/2
    test_imag = np.zeros_like(x)
    theta_good =  -np.pi*1/size*(Ny-kint[1::2])
    test_imag[1::2]= np.tan(theta_good)
    test = test_real+1j*test_imag

    plt.plot(kint[1::2], f2hat.imag[1::2],c='k')
    plt.plot(kint[1::2], test_imag[1::2],c='r')
    plt.savefig('p49c_plots/square_correct.png')
    plt.clf()
    err=f2hat.imag[1::2]-test_imag[1::2]
    plt.plot(kint[1::2], err,c='k')
    plt.savefig('p49c_plots/square_correct_error.png')

if 0:
    size=30
    x=np.arange(0,1,1./size)
    f3a = (x<0.5).astype('float')
    f3ahat = np.fft.fft(f3a)
    Ny = size//2
    k = np.fft.fftfreq(size,d=1./size)
    kodd = k[1::2]
    f3hat = np.zeros(size,dtype='complex')
    f3hat[0]=size/2
    theta3=  -np.pi*1./size*(Ny-np.abs(kodd))*np.sign(kodd)
    f3hat[1::2] = 1+1j*np.tan(theta3)
    f3 = np.fft.ifft(f3hat)
    plt.clf()
    plt.plot(f3.real,c='k')
    plt.plot(f3.imag,c='m')
    plt.savefig('p49c_plots/square_sorted.png')



if 1:
    f2 = np.zeros_like(x)
    #f2 = np.logical_and(x>0.25,x<0.75).astype('float')
    start = x.size//3
    f2[start:x.size//2+start] = 1
    f2hat = np.fft.fft(f2)
    theta = np.angle(f2hat[1::2])
    print(Ny-theta/-np.pi*size)
    print(np.tan(theta)*f2hat.real[1::2]/f2hat.imag[1::2])
    print(f2hat.imag[1::2] - np.tan(theta))
    print(f2hat.imag[1::2] / np.tan(-(Ny-kint[1::2])/size*np.pi  ))
    test_real = np.zeros_like(x)
    test_real[1::2]=1
    test_real[0] = size/2
    test_imag = np.zeros_like(x)
    theta_good =  -np.pi*1/size*(Ny-kint[1::2])
    test_imag[1::2]= np.tan(theta_good)
    test = test_real+1j*test_imag

    fig,ax = plt.subplots(2,2)
    ax0=ax[0][0]
    ax1=ax[1][0]
    ax2=ax[0][1]
    ax3=ax[1][1]
    ax0.plot(kint,np.abs(f2hat),c=[0.5]*4)
    ax0.plot(kint,np.real(f2hat),c='r')
    ax0.plot(kint,np.imag(f2hat),c='g')
    ax0.set_title('mag,real,imag,theta')

    ax0.plot(kint,np.angle(f2hat)/np.pi,c='c')
    #plt.plot(kint,np.angle(f2hat),c='c')


    theta_good=np.zeros_like(f2hat)
    theta_good[1:] =  -np.pi*1/size*(Ny-kint[1:])
    #ax1.set_title('tan,cos,sin theta')
    #ax1.plot(kint,np.tan( theta_good),'c:')
    #ax1.plot(kint,np.cos( theta_good),'r:')
    #ax1.plot(kint,np.sin( theta_good),'g:')

    mask = np.logical_or(np.abs(f2hat.imag)>1e-11 , np.abs(f2hat.real)>1e-11)
    mask = np.logical_and(mask,kint<Ny)
    ax2.plot(kint[mask],np.abs(f2hat)[mask],c=[0.5]*4)
    ax2.plot(kint[mask],np.real(f2hat)[mask],c='r')
    ax2.plot(kint[mask],np.imag(f2hat)[mask],c='g')
    ax2.set_title('mag,real,imag,theta')
    ax2.plot(kint[mask],np.angle(f2hat[mask])/np.pi*size,c='c')

    mask = np.logical_or(np.abs(f2hat.imag)>1e-11 , np.abs(f2hat.real)>1e-11)
    model = np.zeros_like(f2hat)
    model[mask] = 1 + 1j*np.tan(theta_good[mask])
    shift = np.zeros_like(f2hat)
    shift[mask] = np.cos(-2*np.pi*kint[mask]*start/size)+1j*np.sin(-2*np.pi*kint[mask]*start/size)
    model[mask] *= shift[mask]

    #ax3.plot(kint[mask],model[mask].real,c='r')
    #ax3.plot(kint[mask],model[mask].imag,c='g')
    f3b=np.fft.ifft(model)
    ax3.plot(x,f3b,c='k')
    print(shift[mask])



    fig.savefig('p49c_plots/square_almost.png')
