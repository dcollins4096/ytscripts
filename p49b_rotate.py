if 'ef' not in dir():
    execfile('go')
    for i in range(3):
        print("====================")
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



this_system = p49_eigen.waves(bx=1.0,by=0.0,bz=0.0,P=0.6,this_wave='a+')
#right_fast.perturb(base_size=nar([16]*3),pert= 1e-6,wave= 6,directory= '/Users/dcollins/scratch/Paper49b_play/r17c_check')
#k1 = [1,1,0]
#k2 = [0,0,1]
k_test = nar([[2.,1.],[1.,1.],[0.,1]])
ampl = nar([1e-6,1e-6])
#k_test = nar([[1.],[1.],[0.]])
kint = k_test.astype(int)
#this_system.rotate(k_test)
this_system.rot_write(k_test,base_size=nar([16]*3),pert=ampl,directory='/Users/dcollins/scratch/Paper49b_play/r101_two')


"""
size = [16,16,16]
all_hats = {}
field_list = ['d','vx','vy','vz','bx','by','bz','e']
for f in field_list:
    all_hats[f]=np.zeros(size)*1j
    all_hats[f][kint[0,...],kint[1,...],kint[2,...]] = this_system.rot[f]*ampl
    all_hats[f][-kint[0,...],-kint[1,...],-kint[2,...]] = (all_hats[f][kint[0,...],kint[1,...],kint[2,...]] ).conj()

af={}
cutname = '_z0'
def cut(arr):
    return arr[:,:,0].reshape(16,16)
def nz(arr):
    return  np.where(np.abs(arr) > 1e-9)
def nonzero(arr):
    return arr[ nz(arr)]
for f in ['bz']: # field_list:
    tmp=np.fft.ifftn(all_hats[f])
    real_mean = np.mean(np.abs(tmp.real))
    imag_mean = np.mean(np.abs(tmp.imag))
    if imag_mean/real_mean > 1e-9:
        print("Warning: large imaginary component")
    af[f]=tmp.real

    plt.clf()
    plt.imshow(cut(af[f]))
    outname = 'p49b_rot_%s%s.png'%(f,cutname)
    plt.savefig(outname)
    print(outname)

"""


#I don't know why you would need this, but this is what you'd do.
#k_all = np.mgrid[0.:size[0],0.:size[1],0.:size[2]]
#for dim in range(3):
#    k_all[dim,kint[0,...],kint[1,...],kint[2,...]] = 1


#k[14,8] = -2048j #does what the other thing does.
#k[-14,-8]=2048j  #
#k[-14,8] = 2048j # pi rotation
#k[14,-8] = -2048j# pi rotation
#khathat = np.fft.ifftn(k).real
#plt.clf()
#plt.imshow(khathat)
#plt.savefig('p49b_ifft.png')

"""
B0 = this_system['b0']
B0_mag = (B0**2).sum()**0.5
B_unit = B0/B0_mag
size = nar([4,4,4]) #[64,64,64]
dx = 1./size
#k_all = np.mgrid[0.:size[0],0.:size[1],0.:size[2]]
knorm = (k_all[0,...]**2+k_all[1,...]**2+k_all[2,...]**2)**0.5
ok=knorm>0
a_unit = np.zeros_like(k_all)
for dim in range(len(size)):
    a_unit[dim,...][ok] = k_all[dim,...][ok]/knorm[ok]

#b_unit
b_unit = np.zeros_like(k_all)
B0_dot_a = B_unit[0]*a_unit[0,...]+B_unit[1]*a_unit[1,...]+B_unit[2]*a_unit[2,...]
for dim in range(len(size)):
    b_unit[dim,...] = -B0_dot_a*a_unit[dim,...]+B_unit[dim]


b_norm = (b_unit[0,...]**2+b_unit[1,...]**2+b_unit[2,...]**2)**0.5
ok=b_norm>0
b_unit_tmp = np.zeros_like(b_unit)
for dim in range(len(size)):
    b_unit[dim,...][ok] = b_unit[dim,...][ok]/ b_norm[ok]

along = np.logical_and( np.abs(a_unit[0,...]-b_unit[0]) < dx[0],\
                        np.abs(a_unit[1,...]-b_unit[1]) < dx[1])
along = np.logical_and( np.abs(a_unit[2,...]-b_unit[2]) < dx[2],\
                        along)
#print("should be 64", along.sum())

#c_unit = a X b
c_unit = np.zeros_like(k_all)
c_unit[0,...] = a_unit[1,...]*b_unit[2,...]-a_unit[2,...]*b_unit[1,...]
c_unit[1,...] = a_unit[2,...]*b_unit[0,...]-a_unit[0,...]*b_unit[2,...]
c_unit[2,...] = a_unit[0,...]*b_unit[1,...]-a_unit[1,...]*b_unit[0,...]
#
B0hat = np.zeros_like(a_unit)
B0hat[0,...] = B_unit[0]*a_unit[0,...]+B_unit[1]*a_unit[1,...]+B_unit[2]*a_unit[2,...]
B0hat[1,...] = B_unit[0]*b_unit[0,...]+B_unit[1]*b_unit[1,...]+B_unit[2]*b_unit[2,...]
B0hat[2,...] = B_unit[0]*c_unit[0,...]+B_unit[1]*c_unit[1,...]+B_unit[2]*c_unit[2,...]

#this_system = p49_eigen.waves(bx=1.0,by=0.0,bz=0.5,P=0.6,this_wave='a+')
scalars={}
for field in ['d','P','e']:
    scalars[field]=this_system.quan[field]

hat_system = p49_eigen.waves(bx=B0hat[0,...], by=B0hat[1,...], bz=B0hat[2,...],**scalars)
wave='a+'
dvx  = a_unit[0,...]*this_system.right[wave]['vx']
dvx += b_unit[0,...]*this_system.right[wave]['vy']
dvx += c_unit[0,...]*this_system.right[wave]['vz']

dvy  = a_unit[1,...]*this_system.right[wave]['vx']
dvy += b_unit[1,...]*this_system.right[wave]['vy']
dvy += c_unit[1,...]*this_system.right[wave]['vz']

dvz  = a_unit[2,...]*this_system.right[wave]['vx']
dvz += b_unit[2,...]*this_system.right[wave]['vy']
dvz += c_unit[2,...]*this_system.right[wave]['vz']

dbx  = a_unit[0,...]*this_system.right[wave]['bx']
dbx += b_unit[0,...]*this_system.right[wave]['by']
dbx += c_unit[0,...]*this_system.right[wave]['bz']

dby  = a_unit[1,...]*this_system.right[wave]['bx']
dby += b_unit[1,...]*this_system.right[wave]['by']
dby += c_unit[1,...]*this_system.right[wave]['bz']

dbz  = a_unit[2,...]*this_system.right[wave]['bx']
dbz += b_unit[2,...]*this_system.right[wave]['by']
dbz += c_unit[2,...]*this_system.right[wave]['bz']


#print("a_unit",a_unit.flatten())
#print("b_unit",b_unit.flatten())
#print("c_unit",c_unit.flatten())

"""
"""
plt.clf()
plt.quiver(k_all[0,:,:,0],k_all[1,:,:,0],a_unit[0,:,:,0],a_unit[1,:,:,0])
plt.quiver(k_all[0,:,:,0],k_all[1,:,:,0],b_unit[0,:,:,0],b_unit[1,:,:,0])
plt.savefig('p49_units.png')
"""
#a_unit = khat
#sinTheta b_unit = Bvec - cosTheta a_unit
#b_unit = b_uni/||b_unit||
#c_unit = a_unit cross b_unit
