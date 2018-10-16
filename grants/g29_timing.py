from go import *
b = np.fromfile('p10_r193_cycle_time', dtype='float', sep=' ')
b.shape =  (1186,2)
T = b[:,1]
N = b[:,0]
dt = T[1:] - T[:-1]
tavg = 0.5*(T[1:] + T[:-1])
Ncy = N[1:] - N[:-1]
dtdn = dt/Ncy
import davetools
reload(davetools)
plt.clf()
print(dt)
davetools.dumb_plt(plt,tavg,dtdn, 'T','<dt>/<n>', 'g29_timesteps.png' , scale=('linear','log'))
