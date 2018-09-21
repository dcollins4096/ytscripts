from scipy import signal
from scipy.optimize import curve_fit

def collect_extrema(a,b=None):
    """b collects the extrema of a"""
    if b is None:
        b=np.array([min(a),max(a)])
    b[0] = min([min(a),min(b)])
    b[1] = max([max(a),max(b)])
    return b
def collect_extrema_n(a,b=None):
    """b collects the extrema of a"""
    if b is None:
        b=np.array([a.min(),a.max()])
    b[0] = min([a.min(),b.min()])
    b[1] = max([a.max(),b.max()])
    return b
def gauss_me(x, off,a, b):
    return off * np.exp(-(((x-a)**2)/(2*(b**2))))
def gauss_fit(centers,this_spectra, fit_fwhm_only=True):
    total = this_spectra.sum()
    vbar_est= (this_spectra*centers).sum()/total
    sigma2 = (this_spectra*(centers-vbar_est)**2).sum()/total
    sigma_est = (sigma2)**0.5
    norm_est = this_spectra.max()
    popt, pcov = curve_fit(gauss_me, centers, this_spectra, p0=[norm_est,vbar_est,sigma_est] ) 
    return {'vbar_est':vbar_est,'total':total,'sigma_est':sigma_est,'fit_norm':popt[0],'fit_center':popt[1],'fit_width':popt[2]}

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
def make_random_power(Nx,k_min,k_max,slope, var=-1,mean=None, phases = None):
    """Makes a random 1d real-space signal of *Nx* bins, with power spectra of *slope* 
    between *k_min* and *k_max*, and zero power elsewhere.
    *phases* can be passed in for repeated trials with the same random realization."""
    k=np.arange(Nx*1.0)
    k_power=np.arange(Nx*1.0)
    n = np.zeros((Nx,), dtype=complex)
    dk=1
    k_filter =  np.logical_and( k >= k_min, k<k_max)
    k_power[k_filter] =  (np.arange(k_min,k_max))**(slope*0.5)
    if phases is None:
        phases = np.exp(1j*np.random.uniform(0, 2*np.pi, (Nx),))
    n[k_filter] = k_power[k_filter]*phases[k_filter]
    ns = symmetric(n) #forces Ak=conjugate(A(-k)) so that the ifft is real.
    #analytic norm is of this form:  Do the numerical norm, though.
    #if np.abs(slope+1) < 1e-10:
    #    A = 0.5*np.sqrt(var/(np.log(k_max)-np.log(k_min)))
    #else:
    #    A = 0.5*np.sqrt(var/(k_max**(slope+1.)-k_min**(slope+1.))*(slope+1.))
    s = np.fft.ifft(ns).real
    if var > 0:
        s *= np.sqrt(var*s.size/(np.abs(s)**2).sum())
    if mean is not None:
        s += mean - np.mean(s)
    return phases, s
def power_spectrum(actual_powerlaw):
    """P = A k^g.  
    P/P(kmin) = k^g/kmin^g
    log(P/P(kmin)) = g log(k) - g log(kmin)
    g = log(P/P(kmin)/(1-log(kmin)
    """
    shat = np.fft.fft(actual_powerlaw)  
    shat = shat[:shat.size/2]     #only these K modes are present
    power = np.abs(shat)**2; 
    return power
def check_powerlaw(actual_powerlaw):
    """P = A k^g.  
    P/P(kmin) = k^g/kmin^g
    log(P/P(kmin)) = g log(k) - g log(kmin)
    g = log(P/P(kmin)/(1-log(kmin)
    """
    shat = np.fft.fft(actual_powerlaw)  
    shat = shat[:shat.size/2]     #only these K modes are present
    power = np.abs(shat)**2; 
    kmin = np.where( np.abs(power[2:]) > 1e-13*np.abs(power).sum() )[0][0]+2
    k = np.arange(shat.size) #skip the 0 and 1 modes because logs
    p2 = power/power[kmin]
    g = np.log(p2)/(np.log(k) - np.log(k[kmin]))
    return g[kmin+1:], kmin

class powerlaw_parameters():
    def __init__(self,var=1.0,mean=0.1,k_min=2,k_max=30,slope=-1):
        self.var=var
        self.mean=mean
        self.k_min=k_min
        self.k_max=k_max
        self.slope=slope
    def __str__(self):
        args=(self.var,self.slope,self.mean,self.k_min,self.k_max)
        string = r'$r=%0.1f k^{%0.1f}+%0.1f,\ k\in[%0.0f,%0.0f]$'%args
        return string

class amp_phase():
    def __init__(self,k,alpha):
        self.k=k
        self.alpha=alpha

    def __call__(self,a="A", p="R", **kwargs):
        if a == "A":
            self.A(**kwargs)
        if a == "B":
            self.B(**kwargs)
        if a == "C":
            self.C(**kwargs)
        if a == "D":
            self.D(**kwargs)
        if a == "E":
            self.E(**kwargs)
        if a == "F":
            self.F(**kwargs)
        if p == "R":
            self.P1()
        if p == "U":
            self.P2()
    def P1(self):
        Nx = self.k.size
        self.phases = np.exp(1j*np.random.uniform(0, 2*np.pi, (Nx),))
    def P2(self):
        Nx = self.k.size
        self.phases = np.exp(1j*np.zeros( (Nx),))
    def D(self, sigma=0.913):
        """random distribution with with *sigma*, fixed with k"""
        k = self.k
        self.sigma = sigma
        u1 = np.random.uniform(size=k.size)
        u2 = np.random.uniform(size=k.size)
        self.u1=u1
        self.u2=u2
        self.p1=np.sqrt(-2*np.log(u1))
        self.p2 = np.cos(2*np.pi*u2)
        mean = 0
        ph = sigma * self.p1*self.p2
        self.amp = mean + ph
        return 

    def A(self):
        alpha = self.alpha
        k = self.k
        sigma=np.zeros_like(k)
        sigma[1:] = k[1:]**(-0.5*alpha)
        sigma[0]=0
        u1 = np.random.uniform(size=k.size)
        u2 = np.random.uniform(size=k.size)
        self.u1=u1
        self.u2=u2
        self.p1=np.sqrt(-2*np.log(u1))
        self.p2 = np.cos(2*np.pi*u2)
        mean = 0
        ph = sigma * self.p1*self.p2
        self.amp = sigma
        return 
    def E(self,sigma=0.712):
        k=self.k
        u1 = sigma*np.random.normal(size=k.size, scale=sigma)
        u2 = sigma*np.random.normal(size=k.size, scale=sigma)
        self.amp =u1+1j*u2
    def F(self,k_target=0.712):
        k=self.k
        dk = k[1:]-k[:-1]
        omega=np.zeros(k.size)
        #omega[ np.abs(k-k_target) < 3*dk[0] ] = 1.
        omega[ 1 ] = 0.5
        #u1 = sigma*np.random.normal(size=k.size, scale=sigma)
        #u2 = sigma*np.random.normal(size=k.size, scale=sigma)
        self.amp = omega
    def B(self):
        sigma = self.k**(-0.5*self.alpha)
        u1 = nar([np.random.normal( scale=s) + 1j*np.random.normal( scale=s) for s in sigma])
        #ph = np.sqrt(-1*np.log(u1))*np.cos(2*np.pi*u2)
        self.amp = u1
    def C(self,k, alpha):
        sigma = k**(-0.5*alpha)
        k_min = 2; k_max = 5
        k_filter = np.logical_and(k>=k_min, k<=k_max)
        ok = np.ones(k.size)[k_filter]
        u1 = ok*np.random.normal(size=k.size, scale=2.0)
        u2 = ok*np.random.normal(size=k.size, scale=2.0)
        mean = 0
        ph = np.sqrt(-1*np.log(u1))*np.cos(2*np.pi*u2)
        self.amp=u1+1j*u2
        #ph=1.
        #return mean + sigma*ph
