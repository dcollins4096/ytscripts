if 'ef' not in dir():
    execfile('go')
from scipy.interpolate import interp1d
from scipy import interpolate
import fourier_tools.fourier_filter as Filter
from numpy.fft import fftn, ifftn, fftfreq


#set up function to be evaluated and points
misc_from_matlabtricks = nar([0, 1, 10, 2, 0, 0, 4, 5, 3, 1, 0])
x1 = np.arange(len(misc_from_matlabtricks))
x1_fine = np.linspace(0,x1.max(),100)


k = np.arange(0.0,1,0.1)
k_fine = np.arange(0.0,1,0.01)
PS = (k+1e-2)**(-2)
PS[0]=0 #needs to start at zero, or the cdf gets messed up.
#sample_and_plot(k,PS,k_fine,'math_test_PS.pdf')
if 1:
    x=k
    x_fine=k_fine
    prob=PS
    outname = 'math_misc1.pdf'
if 0:
    x=x1
    x_fine=x1_fine
    prob = misc_from_matlabtricks
    outname = 'math_misc1.pdf'
#set up a general pdf.  The important bit is the 
#function "sample_me"
class make_sampler():
    def __init__(self,x,prob,x_fine=None):
        #makes a tool that samples PDF prob, numerically defined at points x.
        #x_fine is a linspace that oversamples prob.
        if x_fine is None:
            x_fine = x
        self.x_fine=x_fine
        self.prob=prob
        self.x=x
        self.spline_tool = interpolate.splrep(x,prob)
        self.pdf = interpolate.splev(x_fine, spline_tool)
        self.pdf=self.pdf/self.pdf.sum()
        neg = pdf<0
        self.pdf[neg]=0
        #set up the cumulative function and inverse
        self.cdf = np.cumsum(self.pdf)
        self.cdf/=self.cdf.max() #needs to go to one.
        self.val,self.ind = np.unique(self.cdf,return_index=True)
        self.sample_me=interp1d(self.val,self.x_fine[ind])
    def __call__(self,rans):
        return self.sample_me(rans)

def hist_plotter(plt,data,bins=100):
    oot=np.histogram(data,bins=100)
    bincenters=0.5*(oot[1][1:]+oot[1][:-1])
    bins = oot[0]*1.0/oot[0].sum()
    plt.plot(bincenters,bins,c='r')
    return oot

def sample_and_plot(x,prob,x_fine,outname='math_test.pdf'):
    plt.clf()
    #sm = sample_me #(x,prob,x_fine)
    sm = make_sampler(x,prob,x_fine)
    #oot=np.histogram(sm( np.random.uniform(size=1200)) ,bins=100)
    data = sm(np.random.uniform(size=1200))
    oot = hist_plotter(plt,data,bins=100)
    dumb_plt(plt,x_fine,pdf,'x','pdf',outname)

#sample_and_plot(x,misc_from_matlabtricks,x_fine)
#sample_and_plot(k,PS,k_fine,'math_test_ps_simple.pdf')

k = np.arange(0.0,1,0.1)
k_fine = np.arange(0.0,1,0.01)
PS = np.zeros_like(k)
PS[1:] = (k[1:])**(-2) #PS needs to start with zero, or the cdf gets mesed up.

sampler = make_sampler(k,PS,k_fine)
these_rands=np.random.uniform(size=120)
rho_power = sampler(these_rands)
if 1:
    rho_hat = np.sqrt(rho_power)
    N=rho_power.size
    rho_hat_hat = np.fft.ifft(rho_hat)
    rho_power_hat = np.fft.ifft(rho_power)
    rho = rho_hat_hat[:N/2]
    dumb_plt(plt,None, rho,'k', 'rho hat', 'math_rhohat.pdf')

    plt.clf()
    dumb_plt(plt,np.arange(len(rho)),rho,'x (arbitrary)','rho','math_rho.pdf')
def power_half(data):
    #data = np.random.rand(301) - 0.5
    ft = np.fft.fft(data)/np.sqrt(data.shape)
    N=data.size
    ps = np.abs(ft)**2
    time_step = 1./data.size #1. / 30
    freqs = np.fft.fftfreq(data.size, time_step)
    idx = np.argsort(freqs)
    k=freqs[0:N/2]
    f=ps[0:N/2]
    return k,f,ft

out=power_half(rho_power_hat)
f=out[0]
p=out[1]
plt.clf()
dumb_plt(plt,f,p,'k','rho(k) back','math_rho_k_return.pdf')



