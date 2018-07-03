ef=execfile
import p49_eigen
import numpy as np
nar = np.array
reload(p49_eigen)

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

twod = True
if twod:
    size = nar([64,64])
    dx=1./size
    kall= p49_eigen.make_k_freqs_2d(size[0])*size[0]

    k = np.zeros(size)*1j
    x,y = np.mgrid[0:1:dx[0], 0:1:dx[1]]
    z=0.
    k_unit = nar([8.,-14, 0])
    kp_unit = k_unit*2*np.pi
else:
    size = nar([64,64,64])
    dx=1./size
    kall= p49_eigen.make_k_freqs(size[0])*size[0]

    k = np.zeros(size)*1j
    x,y,z = np.mgrid[0:1:dx[0], 0:1:dx[1],0:1:dx[2]]
    k_unit = nar([8.,14,-5])
    kp_unit = k_unit*2*np.pi
#k_unit /= (k_unit**2).sum()**0.5
#rho = np.sin(2*np.pi*k_unit[0]*x)+np.sin(2*np.pi*k_unit[1]*y)
#rho=( np.exp( 1j*(k_unit[0]*x+k_unit[1]*y))).real
#rho=( np.exp( 1j*(k_unit[0]*x+k_unit[1]*y))).real
#rho = np.cos(2*np.pi*k_unit[0]*x)
#rho = (np.exp(1j*2*np.pi*k_unit[0]*x)).real
rho = (np.exp(1j*(kp_unit[0]*x+kp_unit[1]*y+ kp_unit[2]*z))).real
cutname = '_y0'
def cut(arr):
    return arr[:,0,:].reshape(64,64)
def nz(arr):
    return  np.where(np.abs(arr) > 1e-9)
def nonzero(arr):
    return arr[ nz(arr)]
def pnz(arr):
    print arr[ nz(arr)]
    nza = nz(arr)
    print nza
    if twod:
        print zip(nza[0],nza[1])
    else:
        print zip(nza[0],nza[1],nza[2])
#plt.clf()
#plt.imshow(cut(rho))
#plt.savefig('p49b_fft3_rho%s.png'%cutname)

def do_more_stuff(msg,arr):
    print("===== %s ======"%msg)
    print(arr.shape)
    print("real sum %0.2e"%np.abs(arr.real).sum())
    print("imag sum %0.2e"%np.abs(arr.imag).sum())
    pnz(arr)
    mask = nz(arr)

    if twod:
        k_nonzero = zip(kall[0,...][mask],kall[1,...][mask])
    else:
        k_nonzero = zip(kall[0,...][mask],kall[1,...][mask],kall[2,...][mask])
    print k_nonzero
rhohat_full=np.fft.fftn(rho)
rhohat_real=np.fft.rfftn(rho)
rhohat_len=np.fft.rfftn(rho,rho.shape)
#i still don't understand what this does.
back_nolen = np.fft.irfftn(rhohat_real)
back_len = np.fft.irfftn(rhohat_real,rho.shape)

do_more_stuff("full fft",rhohat_full)
do_more_stuff("real fft",rhohat_real)
do_more_stuff("real fft, len",rhohat_len)
