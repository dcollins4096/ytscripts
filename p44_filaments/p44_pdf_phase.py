
import p44_random_tools as RT
reload(RT)
Nx = 100
k  = np.arange(Nx,dtype='complex')
x = np.arange(0,1,1./Nx)

Ak = np.zeros(Nx,dtype='complex')
odds = k.real%2 == 1
# -Nx keeps the normalization right, so F(x) = 1,0
# 1/np.pi is for the actual series.
# 1j makes it a sign series.
# Ak[0]=50 keeps the zero-point right (otherwise it's +- 1/2)

Ak[odds] = -Nx/np.pi/k[odds]*1j #1/2 + 2./np.pi/k[odds]
Ak[0] += 50
Ak = RT.symmetric(Ak)
ax = np.fft.ifft(Ak)
#tp = np.abs(ax)
plt.clf()
dumb_plt(plt,None,ax.real,'x','square','p44_square_test.pdf',c='r')
dumb_plt(plt,None,ax.imag,'x','square','p44_square_test.pdf',c='g') #should be zero.

Ak_phase = Ak*np.exp((np.pi/2+np.random.random(Nx))*1j)
Ak_phase = RT.symmetric(Ak_phase)

ax_phase = np.fft.ifft(Ak_phase)
plt.clf()
dumb_plt(plt,None,ax_phase.real,'x','square','p44_square_test_phase.pdf',c='r')
dumb_plt(plt,None,ax_phase.imag,'x','square','p44_square_test_phase.pdf',c='g') #should be zero.

plt.clf()
plt.hist(ax.real,histtype='step',color='r')
plt.hist(ax_phase.real,histtype='step',color='g')
outname='p44_square_test_hist.pdf'
plt.savefig(outname)
print(outname)

plt.clf()
dumb_plt(plt,None,np.abs(Ak_phase[odds]),'k','Ampl(r)angle(g)','p44_square_test_phase_ampl.pdf',c='r')
dumb_plt(plt,None,np.angle(Ak_phase[odds]),'k','Ampl(r)angle(g)','p44_square_test_phase_ampl.pdf',c='g')
