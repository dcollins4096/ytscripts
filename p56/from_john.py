from go import *
from PIL import Image
im = Image.open('p56/autoself/circ.tif')
ar = np.array(im)
f = np.fft.fftn(ar)
f = f*np.conj(f)
f = np.fft.ifftn(f)
f = np.fft.fftshift(f)
f = f.real
f = f/f.max()
plt.clf()
plt.plot(f[33,33:])
plt.savefig('p56_auto.png')
