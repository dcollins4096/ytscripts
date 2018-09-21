# Number of samplepoints
N = 200
# sample spacing
T = 1.0 / 800.0
x = np.linspace(0.0, N*T, N)
y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
y = np.exp( -(x)**2/(2*0.025))
yf = np.fft.fft(y)
yf_shift = np.fft.fftshift(np.abs(yf))
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

plt.close('all')
fig, ax = plt.subplots()
#ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
ax.plot(yf_shift)

fig.savefig('p44_fft_test_2.pdf')
