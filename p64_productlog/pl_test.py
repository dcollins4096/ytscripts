
from go import *
from scipy.special import lambertw  
plt.clf()
x = np.arange(-1,1,0.01)
#plt.plot(x, x/lambertw(x,0),c='k')
#plt.plot(x, x/lambertw(x,-1).real,c='r')
#plt.plot(x, x/lambertw(x,-1).imag,c='g')
#plt.plot(x, x*np.log(x))
W0y=x/lambertw(x,k=0)
W1y=x/lambertw(x,k=-1)
plt.plot(x,W0y.real,c='y')
plt.plot(x,W1y.real,c='m')
plt.plot(x,W1y.imag,c='c')
#plt.plot(x,W0y.imag,c='b')
plt.grid(True)
plt.savefig('../PigPen/p64_test.png')


