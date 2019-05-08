

#Why are my fourier transforms not working?
#https://math.stackexchange.com/questions/2173780/computing-fourier-transform-of-power-law
#This is pretty straight forward; rewrite 1/r^a as a Gamma, but make the Gamma look like
#a Gaussian, which we can fourier tranform  easily.
#
from go import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scatter_fit
from scipy.special import gamma, gammaincc
size =8
dx=1./size
x = np.arange(0,1,dx)
mask=x>0
k=np.fft.fftfreq(size)*size
sft=np.fft.fftshift

x1 = 0.5# 1-2./size # 0.75#+1./size
x0 = 0.0# 2./size #0.25#-1./size

mask = np.logical_and(x>=x0,x<=x1)
y = np.zeros_like(x)
#y[mask]=1 #x[mask]**-2.5
power=0.0
y[mask]=x[mask]**power
z = np.fft.fft(y)
k = np.fft.fftfreq(size,d=dx)

#expect = size/(2*np.pi*k)*(np.cos(-2*np.pi*x1*k/size)
#                          -np.cos(-2*np.pi*x0*k/size)
#                          -1j* np.sin(-2*np.pi*x1*k/size)
#                          +1j* np.sin(-2*np.pi*x0*k/size))
if power==1.:
    def inter(xi):
        return  -size*(1/(2*np.pi*1j))**2*(1/k**2+2*np.pi*1j*xi/k)*np.exp(-2*np.pi*1j*k*xi)
elif power==0.:
    def inter(xi):
        return  -size*(1/(2*np.pi*1j))*(1./k)*np.exp(-2*np.pi*1j*k*xi)
else:
    pass
#   def inter(xi):
#       a = 1/size**power
#       a*= 1/(-2*np.pi*1j*k/size)**power
#       return a*gammaincc(1+power,

shift=0.5/size
shift=0
expect = inter(x1+shift)-inter(x0-shift)

expect[0] = np.sum(y)

#expect_real=sum([np.cos(-2*np.pi*d*x) for d in [nd-2,nd-1]])#range(nd)])
expect=sum([np.exp(-2*np.pi*d*x*1j) for d in range(size//2+1)])#range(nd)])

plt.clf()
plt.plot(x,expect.real,c='r')
plt.plot(x,expect.imag,c='g')
x2 = np.arange(0,1,0.01*dx)
expect2=sum([np.exp(-2*np.pi*d*x2*1j) for d in range(size//2+1)])#range(nd)])
plt.plot(x,expect.real,c='r')
plt.plot(x,expect.imag,c='g')
plt.plot(x2,expect2.real,'r:')
plt.plot(x2,expect2.imag,'g:')
plt.savefig('p49c_plots/power_square_3.png')
###
tp=np.pi*2
x3 = np.arange(0,1,1./10)
size3=x3.size
x4 = np.arange(0,1,0.1/size3)
npoints=x3.size//2+1

prange=np.arange(npoints)
rm = rainbow_map(npoints)
partials=[]
fig1,ax1=plt.subplots(1,1)
fig2,ax2=plt.subplots(1,1)
Ny = x3.size//2
for d in prange:
    partials.append(-np.sin(tp*d*x3))
    other_thing = -np.sin(tp*(Ny-d)*x4)-np.sin(tp*(Ny+d)*x4)
    #ax1.plot(x4,other_thing,linestyle=":",c=rm(d))
    ax1.plot(x3,partials[-1],marker='*',c=rm(d),label='%d'%d)
    delta = Ny-d
    fakery = -np.cos(-tp*Ny*x3)*np.sin(-tp*delta*x3)
    ax1.plot(x3,fakery,c='c')
    ax2.plot(x4,-np.sin(tp*d*x4),marker='*',c=rm(d),label='%d'%d)
    #plt.plot(x2,np.cos(tp*d*x2),marker='*')
partials=nar(partials)
expect3=sum([np.exp(-2*np.pi*d*x3*1j) for d in range(npoints)])#range(nd)])
expect4=sum([np.exp(-2*np.pi*d*x4*1j) for d in range(npoints)])#range(nd)])
ax1.plot(x3,expect3.imag,c='k',label='sum')
ax2.plot(x4,expect4.imag,c='k',label='sum')
#ax1.legend(loc=0)
fig1.savefig('p49c_plots/square_tired.png')
fig2.savefig('p49c_plots/square_tired_h.png')
ok =  expect3.imag< -1e-11
xok = x3[ok]
yok = -expect3.imag[ok]
integrals= [(np.cos(np.pi*d)-1)/(tp*d/size3) for d in prange]
plt.clf()
#plt.plot(xok,yok,marker='*',c='r')
#plt.plot(xok,100/np.pi*((xok/xok[0]))**-1,marker='*')
plt.plot(xok,-yok+(100/np.pi*((xok/xok[0]))**-1),marker='*')
#plt.xscale('log');
#ppltlt.yscale('log')
plt.savefig('p49c_plots/nyquist_imag.png')

###


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
#more_dumb = np.zeros(x.size,dtype='complex')
#moredumb_k=np.arange(x.size)
#for ki in moredumb_k:
#    more_dumb[ki]=np.sum(dumb(y,ki,do_plot=False))
##plt.savefig('p49c_plots/square_dumb.png')



#for n,v in enumerate(z):
#    print("%5d : %9.2f + %9.2f"%(n,v.real,v.imag))
fig,ax = plt.subplots(1,2)
ax0=ax[1]
ax1=ax[0]
#ax0.plot(x,y)
ax1.plot(z.real,c='r',label='real')
ax1.plot(z.imag,c='g',label='imag')
ax1.plot(np.abs(z),c='k',label='abs')
##ax1.plot(moredumb_k,more_dumb.real,c='y')
##ax1.plot(moredumb_k,more_dumb.imag,c='c')
ax1.plot(expect.real,c='m',label='xr')
ax1.plot(expect.imag,c='c',label='xi')
ax1.plot(np.abs(expect),c=[0.5]*3,label='xa')

ok = slice(None)
ok=np.abs(z.real)>0
ax0.scatter((x*size)[ok],z.real[ok]-expect.real[ok])
ok = np.abs(z.imag) > 0
ax0.plot((x*size)[ok],z.imag[ok]-expect.imag[ok])
ok = np.abs(z.imag)>0
ax0.plot((x*size)[ok],np.abs(z)[ok]-np.abs(expect)[ok])
#ax0.plot(x*size,-0.1*np.sin(2*np.pi*x),c='k')
ax1.legend(loc=0)
fig.savefig('p49c_plots/power_square.png')
plt.close(fig)
"""
fig,ax = plt.subplots(1,2)
ax0=ax[1]
ax1=ax[0]
diff = z.imag-expect.imag
ok = np.logical_and(k>0,diff>0)
kok=k[ok]
dok=diff[ok]
fit=scatter_fit.scatter_fit(ax1 ,kok,dok)
ax0.plot(kok,dok)
fig.savefig('p49c_plots/square_error.png')
plt.close(fig)

fig,ax = plt.subplots(1,2)
ax0=ax[1]
ax1=ax[0]
ok=k>=0
ax0.plot(np.abs(expect)[ok],c=[0.5]*3,label='xa')
ax0.plot(np.abs(z)[ok],c='k',label='abs')
ok2 = np.logical_and(k>0, k<10)
fitk=k[ok2]
fitabs=np.abs(z)[ok2]
scatter_fit.scatter_fit(ax0,fitk,fitabs,plot_text=True)
powerline(ax0,fitk[0],fitk[-1],fitabs[0],power-1)
ax0.set_xscale('log'); ax0.set_yscale('log')
#ax1.plot(k[k>0],size/k[k>0]/np.pi,c='k')
#ax1.plot(np.abs(z),c=[0.5]*4)
ax1.legend(loc=0)
fig.savefig('p49c_plots/power_square_slopes.png')
plt.close(fig)

#why are you stupid?
fig2,ax2 = plt.subplots(1,1)

z2 = z+0
z2[np.abs(z)==1]=0.+0.j
y2 = np.fft.ifft(z2)
y3 = np.fft.ifft(expect)
#ax2.plot(y2,c='r',marker='*',)
ax2.plot(y,c='g',marker='*',label='actual')
ax2.plot(y3,c='c',marker='*',label='expect')
ax2.legend(loc=2)

#ax2.plot(y2-1./size,c='k',marker='*')

fig2.savefig('p49c_plots/power_square_oneiszero.png')
"""
