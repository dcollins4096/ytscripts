

def g_disk(z):
    #kpc/Myr^2 net
    a1 = 1.42e-3
    a2 = 5.49e-4
    a3 = 5.0e-5
    z0 = 0.18 #kpc
    g = -a1*z*(z**2+z0**2)**(-0.5) - a2*z + a3*z*np.abs(z)
    return g


z = np.arange(0,3,0.01)
plt.clf()
plt.plot([0.18,0.18], [-0.0015, -0.0005])
dumb_plt(plt,z,g_disk(z), 'z [kpc]', 'g [kpc/Myr^2]', 'p06_vert_test.pdf', clobber=False)
