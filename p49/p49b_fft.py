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

if 0:
    size = [64,64]
    k = np.zeros(size)*1j
    y,x = np.mgrid[0:1:1./64, 0:1:1./64]
if 0:
    size = 16
    #k = np.zeros(size)*1j
    y,x = np.mgrid[0:1:1./size, 0:1:1./size]
if 1:
    size = 16
    #k = np.zeros(size)*1j
    if 1:
        dx = 1./size
        end = 1+dx
    if 0:
        dx = 1./(size-1)
        end = 1
    y,x = np.mgrid[0:end:dx,0:end:dx]
kp_unit = nar([0,1])*2*np.pi
rho = (np.exp(1j*(kp_unit[0]*x+kp_unit[1]*y))).imag
if 0:
    kp_unit = nar([2,1])*2*np.pi
    rho += (np.exp(1j*(kp_unit[0]*x+kp_unit[1]*y))).imag


rhohat=np.fft.rfftn(rho)
plt.imshow(rhohat.real)
plt.imshow(rhohat.imag)
import turb_quan
reload(turb_quan)
turb_quan.plotter2([rho,rho,rho],
                   labs=['rho','rhohat','rhohat i','rhohat real'],
                   norm='ind',
                  fname = 'p49c_plots/2dgames_xspace.png', share=False)
turb_quan.plotter2([np.abs(rhohat)],npx=1,zmin=1e-11,
                   labs=['rhohat','rhohat i','rhohat real'],
                   norm='positive',
                  fname = 'p49c_plots/2dgames_kspace.png', share=True)
if 0:
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
