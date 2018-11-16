
if 'car' not in dir():
    car=taxi.taxi('p49d_a02')
#car.plot()
if 'frame' not in dir():
    frame=3
ds = car.load(frame)
ad = ds.all_data()
dv = ad['cell_volume']
vtotal = dv.sum()
mean_density=(ad['density']*dv).sum()/vtotal
var_density=( (ad['density']-mean_density)**2*dv).sum()/vtotal
sqr_density=( (ad['density'])**2*dv).sum()/vtotal
std_density = np.sqrt(var_density)
print("Mean density %0.7e"%(mean_density))
print("Msq  density %0.7e"%(mean_density**2))
print("var  density %0.7e"%(var_density))
print("std  density %0.7e"%(std_density))
print("sqr  density %0.7e"%(sqr_density))
print("std  density %0.7e"%( np.sqrt(sqr_density-mean_density**2)))
