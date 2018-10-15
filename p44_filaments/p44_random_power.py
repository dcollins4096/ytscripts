if 'ef' not in dir():
    execfile('go_lite')
def gaussian(x,p):
    #p = {mean, width, norm, offset}
    return p[2]/np.sqrt(p[1]*2*np.pi)*np.exp(-( (x-p[0])/p[1])**2 ) #-p[3]

def fwhm(arr, x=None):
    if x is None:
        x = np.arange(arr.size)
    half_max = 0.5*(arr.max())
    over_half = arr > half_max
    max_x = x[over_half].max()
    min_x = x[over_half].min()
    max_pos = x[arr==arr.max()]
    if len(max_pos) > 0:
        max_pos = max_pos[0]
    return half_max, max_pos, 1.0*(max_x-min_x)
def L1(data):
    return np.abs(data).sum()
def L2(data):
    return (np.abs(data)**2).sum()
class spectral_slope():
    def __init__(self,field,k, name=""):
        self.field = field
        self.absfield=np.abs(field)
        #pdb.set_trace()
        self.name=name
        self.k=k
        mask = np.logical_and( k > 0, field>0)
        self.mask = mask
        self.mask[mask] = np.logical_and( self.mask[mask], np.log(k[mask]) != 0.0)
        self.mask = np.logical_and( self.mask, np.log10(self.absfield/self.absfield.max()) > -9)
        self.slope = np.zeros_like(field)
        self.slope[self.mask] = np.log(self.field[self.mask])/np.log(self.k[self.mask])
        tmpthing = self.slope[self.mask]
        self.mean_slope = np.mean(tmpthing)
        self.var  = np.var(self.slope[self.mask])
        print self
    def __str__(self):
        return " mean %0.2e var %0.2e %s"%( self.mean_slope, self.var,self.name)

def ratios(c, string):
    print "%0.2e %s"%(np.abs(c.imag).sum()/np.abs(c.real).sum(), string)
def iratios(c, string):
    print "%0.2e %s"%(np.abs(c.real).sum()/np.abs(c.imag).sum(), string)
def symmetric(v):
    s2=np.zeros_like(v)
    Nx=v.size
    s2[1:Nx/2] = v[1:Nx/2]
    #s2[Nx/2+1:] = v[1:Nx/2] #.conj()
    s2[Nx:Nx/2:-1] = v[1:Nx/2].conj()
    s2[0]=v[0]
    return s2
def make_random_power(Nx,k_min,k_max,slope):
    k=np.arange(Nx*1.0)
    k_power=np.arange(Nx*1.0)
    n = np.zeros((Nx,), dtype=complex)
    plt.clf()
#  0 for 1d, 0.5 for 2d, 1 for 3d.  (target_slope - (D-1))/2
    dk=1
    k_filter =  np.logical_and( k >= k_min, k<k_max)
    k_power[k_filter] =  (np.arange(k_min,k_max))**(slope*0.5)
    if 'all_phasesc' not in dir():
        all_phasesc = np.exp(1j*np.random.uniform(0, 2*np.pi, (Nx),))
        all_phases3 = np.exp(1j*np.random.uniform(0, 2*np.pi, (3),))
#all_phases=np.ones(Nx)
    all_phases=all_phasesc
    n[k_filter] = k_power[k_filter]*all_phases[k_filter]
    ns = symmetric(n) #needs to be self-adjoint in the right way.
    #for i in range(Nx/2,Nx):
    #    print ns[i] == ns[-i].conj()
    s = np.fft.ifft(ns)
    return s

class test1():
    def __init__(self,Nx=30,target_slope=-1.0, k_min=2,k_max=None):
        self.Nx=Nx
        self.target_slope=target_slope
        self.k_min = k_min
        self.k_max= k_max
        if self.k_max is None:
            self.k_max=self.Nx
        self.s = make_random_power(self.Nx,self.k_min,self.k_max,self.target_slope)
        #this is the right fft and spectra to take.
        #Since s is now real, only half the values matter.
        self.k=np.arange(self.Nx*1.0)
        self.shat = np.fft.fft(self.s)  
        self.abs_slope=spectral_slope(np.abs(self.shat[:self.Nx/2]),self.k[:self.Nx/2],'slope s') 
        self.power=spectral_slope(np.abs(self.shat[:self.Nx/2])**2,self.k[:self.Nx/2],'slope s') 
        #plt.clf()
        #plt.plot(self.s)
        #plt.savefig('math_tmp2.pdf')

"""
ratios(n,"N should be complex now")
ratios(s,"S should be real now")
#s2=np.zeros_like(s)
#s2[1:Nx/2] = s[1:Nx/2]
#s2[Nx:Nx/2:-1] = s[1:Nx/2].conj()
#s2[0]=s[0]

shat = np.fft.fft(s)
sabs = np.abs(shat)
s_v4a=spectral_slope(sabs,k,'right slope, v complex') 
ratios(sabs,'right slope, v complex')
if 1:
    v1 = np.fft.fft(s)
    v1a = np.abs(v1)
    s_v1a=spectral_slope(v1a,k,'right slope, v complex') 
    ratios(v1,'right slope, v complex')
    plt.plot(k, s_v1a.slope,c='r',label='s complex')
    plt.plot(k, s_v4a.slope,c='g', label= 'symmetric')
    plt.legend(loc=0)
    plt.savefig('math_dumb1.pdf')
#plt.clf()
#plt.plot(k, s_v4a.slope)
#plt.savefig('math_dumb1.pdf')
if 0:
#correct slope, but s is complex.

    v1 = np.fft.fft(s)
    v1a = np.abs(v1)
    spectral_slope(v1a,k,'right slope, v complex') 
    ratios(v1,'right slope, v complex')
#correct slope, but s is complex.
    v2 = np.fft.fft(np.abs(s))
    v2a = np.abs(v2)
    spectral_slope(v2,k,'not-quite-right slope, v real') 
    ratios(v2,'not-quite-right slope, v real')
if 0:
    sabs = spectral_slope(np.abs(s),k,'sabs')
    v = np.fft.fft(s)
    vabs = spectral_slope(np.abs(v),k,'vabs')
    vpow = spectral_slope(np.abs(v)**2,k,'vpow')
    vr = v.real
    vrabs = spectral_slope(np.abs(vr),k,'vrabs')
    vrpow = spectral_slope(np.abs(vr)**2,k,'vrpow')

    sr = np.abs(s)
    v2 = np.fft.fft(sr)
    v2s = spectral_slope(v2,k,'sabs')

#print nslope
#print "k ", np.log(k_power[ok])/np.log(k[ok])
#print "v ", np.log(np.abs(v)[ok])/np.log(k[ok])
#print "v2", np.log((np.abs(v)**2)[ok])/np.log(k[ok])
#
#print "vh ", np.log(np.abs(v)[ok])/np.log(k[ok])
#print "v2", np.log((np.abs(v)**2)[ok])/np.log(k[ok])
#
#oot=power_half(v)
#power=oot[1]
#print np.log(power)/np.log(np.arange(len(power)))
"""
