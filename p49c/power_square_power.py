

#Why are my fourier transforms not working?
#https://math.stackexchange.com/questions/2173780/computing-fourier-transform-of-power-law
#This is pretty straight forward; rewrite 1/r^a as a Gamma, but make the Gamma look like
#a Gaussian, which we can fourier tranform  easily.
#
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from scipy.special import gamma
size =32
dx=1./size
x = np.arange(0,1,dx)
mask=x>0
k=np.fft.fftfreq(size)*size
sft=np.fft.fftshift

x1 = 1-2./size # 0.75#+1./size
x0 = 2./size #0.25#-1./size

mask = np.logical_and(x>=x0,x<=x1)
y = np.zeros_like(x)
#y[mask]=1 #x[mask]**-2.5
y[mask]=x[mask]
z = np.fft.fft(y)
k = np.fft.fftfreq(size,d=dx)

#expect = size/(2*np.pi*k)*(np.cos(-2*np.pi*x1*k/size)
#                          -np.cos(-2*np.pi*x0*k/size)
#                          -1j* np.sin(-2*np.pi*x1*k/size)
#                          +1j* np.sin(-2*np.pi*x0*k/size))
expect = -size/(2*np.pi*k*1j)*(np.exp(-2*np.pi*1j*k*x1)-np.exp(-2*np.pi*1j*k*x0))
def inter(xi):
    return  -size*(1/(2*np.pi*1j))**2*(1/k**2+2*np.pi*1j*xi/k)*np.exp(-2*np.pi*1j*k*xi)

expect = inter(x1)-inter(x0)

expect[0] = np.sum(y)
#expect[0]=np.sum(y)
plt.clf()
def dumb(Q,k,do_plot=False):
    x = np.arange(Q.size)
    N = Q.size
    dumbest =  Q*np.exp(-2*np.pi*1j*k*x/N)
    if do_plot:
        plt.plot(1.*x/N,dumbest)
    #ndumb = len(glob.glob('p49_dumb/*png'))
    #tots=np.sum(dumbest)
    return dumbest
more_dumb = np.zeros(x.size,dtype='complex')
moredumb_k=np.arange(x.size)
for ki in moredumb_k:
    more_dumb[ki]=np.sum(dumb(y,ki,do_plot=False))
plt.savefig('p49c_plots/square_dumb.png')



#for n,v in enumerate(z):
#    print("%5d : %9.2f + %9.2f"%(n,v.real,v.imag))
fig,ax = plt.subplots(1,1)
#ax0=ax[0]
ax1=ax#[1]
#ax0.plot(x,y)
ax1.plot(z.real,c='r',label='real')
ax1.plot(z.imag,c='g',label='imag')
#ax1.plot(moredumb_k,more_dumb.real,c='y')
#ax1.plot(moredumb_k,more_dumb.imag,c='c')
ax1.plot(expect.real,c='m',label='xr')
ax1.plot(expect.imag,c='c',label='xi')
#ax1.plot(k[k>0],size/k[k>0]/np.pi,c='k')
#ax1.plot(np.abs(z),c=[0.5]*4)
ax1.legend(loc=0)
fig.savefig('p49c_plots/power_square.png')
plt.close(fig)

#why are you stupid?
fig2,ax2 = plt.subplots(1,1)

z2 = z+0
z2[np.abs(z)==1]=0.+0.j
y2 = np.fft.ifft(z2)
y3 = np.fft.ifft(expect)
ax2.plot(y2,c='r',marker='*')
ax2.plot(y,c='g',marker='*')
ax2.plot(y3,c='c',marker='*')
ax2.plot(y2-1./size,c='k',marker='*')

fig2.savefig('p49c_plots/power_square_oneiszero.png')
