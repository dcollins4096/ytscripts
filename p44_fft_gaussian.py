#From https://stackoverflow.com/questions/5398304/fourier-transform-of-a-gaussian-is-not-a-gaussian-but-thats-wrong-python
import matplotlib.pyplot as plt
import numpy as np

N = 128
k = np.arange(-5,5,10./(2*N))
var=0.25
k0=0
fk = np.exp(-(k-k0)**2/(2*var))
fx=np.fft.ifft(fk)
fx_shift = np.fft.fftshift(np.abs(fx))* np.sqrt(2 * N)

plt.clf()
plt.plot(k,fk,label='fk',c='r')
plt.plot(k,fx_shift,label='fx_shift',c='g')
plt.plot(k,fx.real,label='fx.real',c='b')
plt.legend(loc=0)
plt.savefig('p44_gaussian_test.pdf')

