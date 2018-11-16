exec(open('p52_WD2/potentials.py').read(),globals(),locals())
pot = potential(a01)
plt.clf()
plt.close('all')
x=pot.read(1)
dR = np.mean(x.Rb[1:]-x.Rb[:-1])
rho = x.Ms/(x.Rb**2*4*np.pi*dR)
print(rho)
