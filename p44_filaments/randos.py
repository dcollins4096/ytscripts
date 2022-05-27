from go import *


def gaa(x,mu=0,sigma=None):
    original_curve= np.exp(-(x-mu)**2/(2*sigma**2))
    original_curve/=original_curve.sum()
    return original_curve

x=np.linspace(-100,100,1000)

r=np.random.random(x.size)

G1 = gaa(x,sigma=0.3)

G2 = x*0
G2[ 500-300:500+300]=1

plt.clf()
plt.plot(x)
plt.savefig(os.environ['HOME']+'/PigPen/p44_randos.png')
