
plt.close('all')
a = nar([0,0,0,0,1,1,1,1,0,0,0,0])
Nx=32
a= np.zeros(Nx)
a[Nx//4:3*Nx//4] = 1
ahat = np.fft.fft(a)
cr = np.convolve(a,a)
cr_hat_1 = np.fft.fft(cr)

if 0:
    fig,axes=plt.subplots(1,2)
    ax, ax1=axes
    ax.plot(a,'k')
    ax.plot(cr,'r')

    #ax1.plot(cr_hat_1,'k')
    #ax1.plot(cr_hat_1*cr_hat_1.conj(),'r')

fig,axes=plt.subplots(1,2)
ax0, ax1=axes

X=[1.+1.j, 2.+3.j, 1.+2.j, 3.+2.j]
H=[2.+3.j, 1.+1.j, 3.+3.j, 4.+5.j]
X=a
H=a
#t=np.arange(0,2*np.pi,0.1)
#X = np.sin(2*t)
#H = np.sin(3*t)
t01=0
sigma1=2
X = 1./np.sqrt(np.pi*sigma1**2)*np.exp(-(t-t01)**2/(2*sigma1**2))
H = 1./np.sqrt(np.pi*sigma1**2)*np.exp(-(t-t01)**2/(2*sigma1**2))


import scipy
if 0:
    ax0.plot(X)
    ax0.plot(H)
    Y=np.array(X)*np.array(H)
    y=np.fft.ifft(Y)
    x=np.fft.ifft(X)
    h=np.fft.ifft(H)
    k = np.fft.fftfreq(Y.size)
    new_x=np.hstack((x[1:], x))
    ax0.plot(x,'g')
    new_y = scipy.convolve(new_x,h, 'valid')
    ax1.plot(new_y,'r')
    ax1.plot(y,'k')
    delta = new_y-y
    ax1.plot(delta)
    fig.savefig('p56_convolve_2.png')

if 1:
    ax0.clear()
    ax1.clear()
    ax0.plot(X,'k')
    ax0.plot(H,'b')
    new_X = np.hstack( (X[1:],X))
    AC = scipy.correlate(new_X, H)
    ax0.plot(AC,c='r')
    
    x = np.fft.ifft(X)
    h = np.fft.ifft(H)
    hdag = h.conj()
    XH = np.fft.fft(x*hdag)
    ax1.plot(x,'k')
    ax1.plot(h,'b')
    ax1.plot(x*h,'r')
    ax1.plot(XH,'k:')
    AC_from_ft = XH*Nx
    renorm=(XH*Nx).max()/AC.max()
    renorm=0.5
    AC_from_ft = AC_from_ft/renorm
    ax0.plot(AC_from_ft,'k:')
    

    L = np.sum(AC)/AC[0]
    import matplotlib.patches as patches
    rect=patches.Rectangle((0,0),L,AC[0],facecolor=[0.5]*3)
    ax0.add_patch(rect)


    fig.savefig('p56_convolve_3.png')

