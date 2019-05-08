from go import *
if 'car' not in dir():
    car = taxi.load('p19_sphere')
    car.frames=[1]
    car.fields=['GravPotential']
    car.load()
#car.plot()
#car.profile(['radius','GravPotential'], scales=['linear','linear'])
#car.profile(['radius','density'], scales=['linear','linear'])
x = car.ds.index.grids[0]['x'].v
y = car.ds.index.grids[0]['y'].v
z = car.ds.index.grids[0]['z'].v
r = np.sqrt( (x-0.5)**2+(y-0.5)**2+(z-0.5)**2)

plt.clf()
d = car.ds.index.grids[0]['density']
phi = car.ds.index.grids[0]['GravPotential']
plt.scatter(r, d)
plt.savefig('p19_sphere_r_d.png')
plt.clf()
def indef(a,r):
    #phi(a) = int rho(r)/(r-a) r^2 dr
    return -0.5* r* (2* a + r) - a**2*np.log(np.abs(-a + r))
def ddd(a):
    rmax=0.5
    dif = indef(a,rmax)-indef(a,0)
    return dif
plt.scatter(r,phi)
this_r = r.flatten()
this_r = this_r[ this_r>0]
#plt.plot(this_r, np.pi*this_r**2)
plt.plot(this_r,-1/(6*np.pi*np.sin(this_r)))
plt.plot(this_r,ddd(this_r))
plt.savefig('p19_sphere_r_phi.png')
