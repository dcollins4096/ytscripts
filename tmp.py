if 'ef' not in dir():
    execfile('go')

#import p44_random_tools
#reload(p44_random_tools)
def symmetric2d(v):
    #for real fft, Ak = A(-k)^*; the negative-phase is the conjugate
    #(that way when you sum, the imaginary parts cancel.)
    #This lets us take an arbitrary K-space signal and ensure it's inverse-fft is real.
    s2=np.zeros_like(v)
    Nx=v.size
    s2[1:Nx/2] = v[1:Nx/2]
    s2[Nx:Nx/2:-1] = v[1:Nx/2].conj()
    s2[0]=v[0]
    return s2
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
N1 = 100
b = np.zeros(N1)
b[3]=1
b12=symmetric(b)
b12[0]=0
bhat = np.fft.irfft(b12)
print("mean bhat",np.mean(bhat))
print("var bhat",np.var(bhat))
print("var b12",np.var(b12))
print("var bhat/b12",2*N1*np.var(bhat)/np.var(b12))

Nx = 2000
b2 = np.zeros([Nx,Nx/2+1])
b2[3,3]=1
b3 = np.zeros_like(b2)
b3[1:Nx/2,:] = b2[1:Nx/2,:]
b3[Nx:Nx/2:-1,:]=b2[1:Nx/2,:].conj()
bhat2=np.fft.irfft2(b2)
bhat3=np.fft.irfft2(b3)
print("...")
print("var b3 2d", np.var(b3))
print("var bhat 2d", np.var(bhat3))
R=np.var(bhat3)/np.var(b3)
print("var hat/3",Nx*Nx*R)
plt.clf()
plt.imshow(bhat3)
plt.savefig('tmp.png')



if 0:
    plt.clf()
    plt.plot(k,PAk,c='r')
    plt.plot(k,abs(Aharm)**2*Deltak/2/np.pi,c='g')
    plt.title('PAk model vs output')
    plt.yscale('log');plt.xscale('linear')
    plt.savefig('p44_kh_model.png')
if 0:
    Deltai=0.15; Deltaj=0.15
    Delta_x = nar([Deltai,Deltaj])
    Nxi = 10; Nxj = 10
    Total_xi = Nxi*Deltai; Total_xj = Nxj*Deltaj
    np.random.seed(123123)
    X = np.random.randn(Nxi, Nxj)# + 1.0j*np.random.randn(Nxi, Nxj) 

    Delta_ki = 2*np.pi/N/Delta_x # size of pixels in harmonic space
    Nharm_i = N; Nharm_j = N
    #kj, ki = np.mgrid[0:10:
