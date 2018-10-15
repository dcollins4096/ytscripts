
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

size = nar([64,64,64])
N = 32
size = nar([N,N,N//2+1])
k = np.zeros(size)*1j
#k[mask]=1 #this comes from the fft prior
k[5,2,4] = 1+0.2j

#How do yu convert compressible to solenoidal?
##k[[ 5, 59], [14, 50],[ 8, 56]] = 1

#k[14,8] = -2048j #does what the other thing does.
#k[-14,-8]=2048j  #
#k[-14,8] = 2048j # pi rotation
#k[14,-8] = -2048j# pi rotation
#khathat = np.fft.ifftn(k).real
khathat = np.fft.irfftn(k).real
dumbhat = np.fft.fftn(khathat)
cutname = '_y0'
def cut(arr):
    return arr[:,0,:].reshape(N,N)
def nz(arr):
    return  np.where(np.abs(arr) > 1e-9)
def nonzero(arr):
    return arr[ nz(arr)]
plt.clf()
plt.imshow(cut(khathat))
outname = 'p49b_ifft3_%s.png'%cutname
plt.savefig(outname)
print(outname)



"""
a_unit = khat
sinTheta b_unit = Bvec - cosTheta a_unit
b_unit = b_uni/||b_unit||
c_unit = a_unit cross b_unit
"""
