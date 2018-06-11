import enzo_write
reload(enzo_write)
import p49_eigen
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
B_unit = B0/(B0**2).sum()**0.5
size = nar([4,4,4]) #[64,64,64]
dx = 1./size
k_all = np.mgrid[0.:size[0],0.:size[1],0.:size[2]]
knorm = (k_all[0,...]**2+k_all[1,...]**2+k_all[2,...]**2)**0.5
ok=knorm>0
a_unit = np.zeros_like(k_all)
for dim in range(len(size)):
    a_unit[dim,...][ok] = k_all[dim,...][ok]/knorm[ok]

#b_unit
b_unit = np.zeros_like(k_all)
B0_dot_a = B_unit[0]*a_unit[0,...]+B_unit[1]*a_unit[1,...]+B_unit[2]*a_unit[2,...]
for dim in range(len(size)):
    b_unit[dim,...] = -B0_dot_a*a_unit[dim,...]+B0[dim]
b_norm = (b_unit[0,...]**2+b_unit[1,...]**2+b_unit[2,...])**0.5
ok=b_norm>0
for dim in range(len(size)):
    b_unit[dim,...][ok] /= b_norm[ok]

along = np.logical_and( np.abs(a_unit[0,...]-b_unit[0]) < dx[0],\
                        np.abs(a_unit[1,...]-b_unit[1]) < dx[1])
along = np.logical_and( np.abs(a_unit[2,...]-b_unit[2]) < dx[2],\
                        along)
print("should be 64", along.sum())


c_unit = np.zeros_like(k_all)
#0 = 1 2 - 2 1
#1 = 2 0 - 0 2
#2 = 1 0 - 0 1
#c_unit[0,...] = a_unit[1,...]*b_unit[2,...]-a_unit[2,...]*b_unit[1,...]
#c_unit[1,...] = a_unit[2,...]*b_unit[0,...]-a_unit[0,...]*b_unit[2,...]
#c_unit[2,...] = a_unit[0,...]*b_unit[1,...]-a_unit[1,...]*b_unit[0,...]

plt.clf()
plt.quiver(k_all[0,:,:,0],k_all[1,:,:,0],a_unit[0,:,:,0],a_unit[1,:,:,0])
plt.quiver(k_all[0,:,:,0],k_all[1,:,:,0],b_unit[0,:,:,0],b_unit[1,:,:,0])
plt.savefig('p49_units.png')
"""
a_unit = khat
sinTheta b_unit = Bvec - cosTheta a_unit
b_unit = b_uni/||b_unit||
c_unit = a_unit cross b_unit
"""
