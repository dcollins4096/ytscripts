from go_lite_pyt3 import *
import p49_eigen
import numpy as np
nar = np.array
reload(p49_eigen)

cf = '({0.real:5.2e} + {0.imag:5.2e}i)'
def ampl_str(arr):
    out = ""
    for v in arr:
        out += cf.format(v)
    return out
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
    #k_unit = nar([8.,14, 0])
    k_unit = nar([-2.,0, 0])
    kp_unit = k_unit*2*np.pi
else:
    size = nar([64,64,64])
    dx=1./size
    kall= p49_eigen.make_k_freqs(size[0],real=False)*size[0]

    k = np.zeros(size)*1j
    x,y,z = np.mgrid[0:1:dx[0], 0:1:dx[1],0:1:dx[2]]
    k_unit = nar([8.,14,-5])
    #k_unit = nar([1.,1.,1.])
    k_unit = nar([1.,-1.,0.])
    kp_unit = k_unit*2*np.pi
#k_unit /= (k_unit**2).sum()**0.5
#rho = np.sin(2*np.pi*k_unit[0]*x)+np.sin(2*np.pi*k_unit[1]*y)
#rho=( np.exp( 1j*(k_unit[0]*x+k_unit[1]*y))).real
#rho=( np.exp( 1j*(k_unit[0]*x+k_unit[1]*y))).real
#rho = np.cos(2*np.pi*k_unit[0]*x)
#rho = (np.exp(1j*2*np.pi*k_unit[0]*x)).real
rho = (np.exp(1j*(kp_unit[0]*x+kp_unit[1]*y+ kp_unit[2]*z))).real
cutname = '_y0'
strfmt="%20s "
def cut(arr):
    return arr[:,0,:].reshape(64,64)
def nz(arr):
    return  np.where(np.abs(arr) > 1e-9)
def nonzero(arr):
    return arr[ nz(arr)]
def pnz(arr):
    print(strfmt%"nonzero values" +ampl_str( arr[ nz(arr)]))
    nza = nz(arr)
    print(strfmt%"nonzero ind" +str(nza))
    if twod:
        print(strfmt%"nonzero k"+str(list(zip(nza[0],nza[1]))))
    else:
        print(strfmt%"nonzero k"+str(list(zip(nza[0],nza[1],nza[2]))))
if True and twod:
    plt.clf()
    plt.imshow(rho)
    plt.savefig('p49b_fft3_rho_2d_mk.png')
if 0:
    #do some plots
    plt.clf()
    plt.imshow(cut(rho))
    plt.savefig('p49b_fft3_rho%s.png'%cutname)

if 1:
    #ask about non-zero axes
    def do_more_stuff(msg,arr):
        print("===== %s ======"%msg)
        print(strfmt%"shape" + str(arr.shape))
        print(strfmt%"real sum" + "%0.2e"%np.abs(arr.real).sum())
        print(strfmt%"imag sum" + "%0.2e"%np.abs(arr.imag).sum())
        pnz(arr)
        mask = nz(arr)

        if twod:
            k_nonzero = list(zip(kall[0,...][mask],kall[1,...][mask]))
        else:
            k_nonzero = list(zip(kall[0,...][mask],kall[1,...][mask],kall[2,...][mask]))
        print(strfmt%"k vec"+str(k_nonzero))
    rhohat_full=np.fft.fftn(rho)/rho.size
    rhohat_real=np.fft.rfftn(rho)/rho.size*2
    rhohat_len=np.fft.rfftn(rho,rho.shape)/rho.size
#i still don't understand what this does.
    back_nolen = np.fft.irfftn(rhohat_real)
    back_len = np.fft.irfftn(rhohat_real,rho.shape)

    do_more_stuff("full fft",rhohat_full)
    do_more_stuff("real fft",rhohat_real)
    do_more_stuff("real fft, len",rhohat_len)

