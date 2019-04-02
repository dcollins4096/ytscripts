
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

output_dir = "p60_dft_square/plots/"

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
    #dope this works.
    kmy = np.exp(-2*np.pi*1j*kint/size)
    test2 = (1-kmy**(f2.sum()))/(1-kmy) 
    test=test2
    test2_imag=test2.imag
    test2_real=test2.real

    test_imag=test2.imag
    test_real=test2.real

    print("=========================")
    print(test2_imag-test.imag)

    plt.plot(kint[1::2], f2hat.imag[1::2],c='k')
    plt.plot(kint[1::2], test_imag[1::2],c='r')
    plt.savefig('%s/square_correct.png'%output_dir)
    plt.clf()
    err=f2hat.imag[1::2]-test_imag[1::2]
    plt.plot(kint[1::2], err,c='k')
    plt.savefig('%s/square_correct_error.png'%output_dir)

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
    plt.savefig('%s/square_sorted.png'%output_dir)



if 1:
    f2 = np.zeros_like(x)
    #f2 = np.logical_and(x>0.25,x<0.75).astype('float')
    start = 5 #x.size//3
    width = 15#x.size//3
    end = width+start
    f2[start:end] = 1
    f2hat = np.fft.fft(f2)
    theta = np.angle(f2hat[1::2])
    angle=np.angle(f2hat)
    scaled_angle= angle/np.pi*size
    #print(Ny-theta/-np.pi*size)
    #print(np.tan(theta)*f2hat.real[1::2]/f2hat.imag[1::2])
    #print(f2hat.imag[1::2] - np.tan(theta))
    #print(f2hat.imag[1::2] / np.tan(-(Ny-kint[1::2])/size*np.pi  ))
    test_real = np.zeros_like(x)
    mask = slice(None)
    test_real[1::2]=1
    test_real[0] = size/2
    test_imag = np.zeros_like(x)
    theta_good =  -np.pi*1/size*(Ny-kint[1::2])
    test_imag[1::2]= np.tan(theta_good)
    test = test_real+1j*test_imag

    # 4/2/19
    #
    # THIS WORKS EVERYTHING ELSE IS A DISTRATION.
    # it's a geometric series.
    kmy = np.exp(-2*np.pi*1j*kint/size)
    test2 = (1-kmy**(end))/(1-kmy) 
    test2 -= (1-kmy**(start))/(1-kmy) 
    test2[0] = f2.sum()

    fig,ax = plt.subplots(2,2)
    ax0=ax[0][0]
    ax1=ax[1][0]
    ax2=ax[0][1]
    ax3=ax[1][1]
    #ax0.plot(kint,np.abs(f2hat),c=[0.5]*4)
    #ax0.plot(kint,np.real(f2hat),c='r')
    #ax0.plot(kint,np.imag(f2hat),c='g')
    ax0.plot(kint, np.abs(test2),linestyle=":",c=[0.5]*4)
    ax0.plot(kint,np.real(test2),linestyle=":",c='r')
    ax0.plot(kint,np.imag(test2),linestyle=":",c='g')
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
    theta_good[1:] =  -np.pi*1/size*(-(1-width)*kint[1:])
    wut=4/np.pi
    theta_good = -1/wut*np.arccos(np.cos(wut*theta_good))
    model[mask] = 1 + 1j*np.tan(theta_good[mask])
    shift = np.zeros_like(f2hat)
    shift[mask] = np.cos(-2*np.pi*kint[mask]*start/size)+1j*np.sin(-2*np.pi*kint[mask]*start/size)
    model[mask] *= shift[mask]
    ax1.plot(theta_good)
    ax1.plot(angle,c='k')

    #ax3.plot(kint[mask],model[mask].real,c='r')
    #ax3.plot(kint[mask],model[mask].imag,c='g')
    f3b=np.fft.ifft(test2)
    ax3.plot(x,f3b,c='k')
    #ax3.plot(x,f2,c=[0.5]*4)
    print(shift[mask])



    fig.savefig('%s/square_almost.png'%output_dir)
