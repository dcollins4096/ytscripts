
Nx=300
k=np.arange(Nx)
slope=1.
t = np.arange(Nx)
n = np.zeros((Nx,), dtype=complex)
k_power = np.zeros(Nx)
k_target = np.zeros(Nx)
k_min = 2
k_max= Nx
dk=1
k_filter =  np.logical_and( k >= k_min, k<k_max)
ReducedDim=-slope*0.5-1 #3d
ReducedDim=-slope*0.5   #1d
k_power[k_filter] = (np.arange(k_min,k_max))**(ReducedDim)
#k_power[k_filter]/=k_power[k_filter].max()
k_target[k_filter] = (np.arange(k_min,k_max))**(-slope)
#n[k_max] = np.exp(1j*(np.pi/180*43)) #np.random.uniform(0, 2*np.pi)) #this works.
if 'all_phases' not in dir():
    all_phases_c = np.exp(1j*np.random.uniform(0, 2*np.pi, (Nx),))
    all_phases_o =np.ones(Nx)
all_phases=all_phases_o
n[k_filter] = k_power[k_filter]*all_phases[k_filter]
#n=k_power
#n[k_filter] = k_power[k_filter]
s = np.fft.ifft(n)
v = np.abs(s)
#plt.clf()
#plt.plot(t, s.real, 'b-', t, s.imag, 'r--',t,v,'k:')
#plt.legend(('real', 'imaginary', 'norm'))
#plt.savefig('math_power2_vis.pdf')

def power_half(data):
    #data = np.random.rand(301) - 0.5
    ft = np.fft.fft(data) #/np.sqrt(data.shape)
    N=data.size
    ps = np.abs(ft)**2
    time_step = 1./data.size #1. / 30
    freqs = np.fft.fftfreq(data.size, time_step)
    idx = np.argsort(freqs)
    #k=freqs[0:N/2]
    #p=ps[0:N/2]
    k=freqs
    p=ps
    return k,p,ft

ef('davetools.py')
#fig, (ax1, ax2) = plt.subplots( 1,2) #, sharey=True)
fig, (ax1, ax2,ax3) = plt.subplots( 3) #, sharey=True)
#fig,  ax2 = plt.subplots( 1, sharex=True)
out=power_half(v)
f=out[0]
print "P max", out[1].max()
p=out[1]# * n.max()/out[1].max()#/out[1].max()
ft=out[2]

#plt.plot(k[1:],k_power[1:])
#ax2.plot(f,p)
#ax2.figure.savefig('math_power2_good.pdf')
#fig.savefig('math_power2_good.pdf')
ax2_scale = ('linear','linear')
#ok = np.logical_and(np.abs(p) < 1e12)
ok = np.abs(p) > 1e-15

this_x1=f[ok]
this_y1=np.log(np.abs(p[ok]))/np.log(f[ok])
dumb_plt(ax1,this_x1,this_y1,'k','log(p)/log(k)','math_power2_good.pdf',scale=('linear','linear'),label='P(v)',c='g',marker="*")
ok = np.logical_and(np.log(k) > 0, np.abs(n) > 1e-3)

this_x3,this_y3=k[ok],np.log(np.abs(k_target[ok]))/np.log(k[ok])
dumb_plt(ax1,this_x3,this_y3,'k','log(p)/log(k)','math_power2_good.pdf',scale=('linear','linear'),label='P(v)',c='b',marker="*")

this_x2,this_y2=k[ok],np.log(np.abs(n[ok]))/np.log(k[ok])
dumb_plt(ax1,this_x2,this_y2,'k','log(p)/log(k)','math_power2_good.pdf',scale=('linear','linear'),label='P(v)',c='r',marker="*")

dumb_plt(ax2,f,np.abs(p)*np.abs(n).max()/np.abs(p).max(),'k','rho(k) back','math_power2_good.pdf',scale=('linear','log'),label='P(v)',c='g',marker="*")
dumb_plt(ax2,None,n.real,'k','rho(k) back','math_power2_good.pdf',scale=('linear','log'),label='n',c='r',linestyle='--')
dumb_plt(ax2,None,n.imag,'k','rho(k) back','math_power2_good.pdf',scale=('linear','log'),label='n',c='r',linestyle=':')
dumb_plt(ax2,None,np.abs(n),'k','rho(k) back','math_power2_good.pdf',scale=('linear','log'),label='|n|',c='r')
#ax3.scatter([Nx/k_max],[0.0] )
dumb_plt(ax3,None,v,'x','v(x)','math_power2_good.pdf')
#ax2.set_ylim(1e-3,0.2)
#ax2.set_xlim(0,50)
ax1.set_xlim(0,50)
ax1.set_ylim(-20,0)
dumb_plt(ax2,None,k_target*p.max()/k_target.max(),'k','rho(k) back','math_power2_good.pdf',scale=('linear','log'),label='target',c='k')
#dumb_plt(ax1,None,v,'x','v','math_power2_good.pdf',scale=('linear','linear'))
plt.close(fig)
