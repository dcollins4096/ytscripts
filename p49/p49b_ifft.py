
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
#k[14,8] = -2048j #does what the other thing does.
#k[-14,-8]=2048j  #
#k[-14,8] = 2048j # pi rotation
#k[14,-8] = -2048j# pi rotation
khathat = np.fft.ifftn(k).real
plt.clf()
plt.imshow(khathat)
plt.savefig('p49b_ifft.png')

d=1.
vx=vy=vz=0
Bx, By, Bz = 1.0,1.41421,0.5
P = 0.6
Gamma=1.6666666667
Egas=P/((Gamma-1)*d)
E = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(Bx*Bx+By*By+Bz*Bz)/d + Egas
right = p49_eigen.eigen(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)
speeds = p49_eigen.speeds(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)

B0 = nar([1.,0.,0.])
size = [64,64,64]
kx, ky, kz = np.mgrid[0.:shape[0],0.:shape[1],0.:shape[2]]


"""
a_unit = khat
sinTheta b_unit = Bvec - cosTheta a_unit
b_unit = b_uni/||b_unit||
c_unit = a_unit cross b_unit
"""
