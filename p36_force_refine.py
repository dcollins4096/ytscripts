

P = 0.6
rho = 10
gamma = 5/3
G = 100
dx = 1./16/2
S = 4

Lj  = np.sqrt(4*3.14*3.14*(gamma*P/rho)/( G*rho))


print S*dx > Lj
print S*dx, Lj
