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

size = [64,64]
k = np.zeros(size)*1j
y,x = np.mgrid[0:1:1./64, 0:1:1./64]
k_unit = nar([8.,14])
kp_unit = k_unit*2*np.pi
#k_unit /= (k_unit**2).sum()**0.5
#rho = np.sin(2*np.pi*k_unit[0]*x)+np.sin(2*np.pi*k_unit[1]*y)
#rho=( np.exp( 1j*(k_unit[0]*x+k_unit[1]*y))).real
#rho=( np.exp( 1j*(k_unit[0]*x+k_unit[1]*y))).real
#rho = np.cos(2*np.pi*k_unit[0]*x)
#rho = (np.exp(1j*2*np.pi*k_unit[0]*x)).real
rho = (np.exp(1j*(kp_unit[0]*x+kp_unit[1]*y))).imag
plt.clf()
plt.imshow(rho)
plt.savefig('p49b_test.png')

rhohat=np.fft.fftn(rho)
plt.clf()
plt.imshow(rhohat.real)
plt.savefig('p49b_real.png')
plt.clf()
plt.imshow(rhohat.imag)
plt.savefig('p49b_imag.png')

def su(x):
    y=np.abs(x).sum()
    if y < 1e-8:
        return "--"
    return "%0.1e"%y
cn = rainbow_map(64)
plt.clf()
for n in range(64):
    plt.plot(rhohat.real[:,n],c=cn(n))
    print( "n %3d r %8s i %8s"%(n,su(rhohat.real[:,n]), su(rhohat.imag[:,n])))
plt.savefig('p49b_rhohat_real.png')
plt.clf()
for n in range(64):
    plt.plot(rhohat.imag[:,n],c=cn(n))
plt.savefig('p49b_rhohat_imag.png')

#k[1,2] = 1+0j
#khat = np.fft.ifft(k)
#x
