from go import *

#@yt.particle_filter(requires=["dynamical_time"], filtered_type='all')
#def big_dust(pfilter, data):
#    filter = data[(pfilter.filtered_type, "dynamical_time")]  > 0.0005
#    return filter

if 'car' not in dir():
    #car = taxi.taxi('fd07')
    car = taxi.taxi('fgd01')
    car.particle_filter_names=['big_dust', 'small_dust']

def bd(pfilter, data):
    filter = data[(pfilter.filtered_type, "dynamical_time")] > 0.0005
    return filter
yt.add_particle_filter("big_dust", function=bd, filtered_type='all', requires=["dynamical_time"])
def sd(pfilter, data):
    filter = data[(pfilter.filtered_type, "dynamical_time")] < 0.0005
    return filter

yt.add_particle_filter("small_dust", function=sd, filtered_type='all', requires=["dynamical_time"])
for n in [20]:
    car.frames=[n]
    car.plot()
    car.the_plot.ds.periodicity = (True,True,True)
    car.the_plot.set_cmap('density','gray')
    #car.the_plot.set_zlim('density',0.9,1.0)
    car.the_plot.annotate_particles(1,ptype='small_dust',col='r')
    car.the_plot.annotate_particles(1,ptype='big_dust',col='g')
    print(car.the_plot.save('%s_dsb_%04d'%(car.outname,n)))

for n in [0,10,20]:
    ds = car.load(n)
    ad=ds.all_data()
    stat(ad['big_dust','particle_velocity_magnitude'],'big vy')
    stat(ad['small_dust','particle_velocity_magnitude'],'small vy')
    stat(ad['all','particle_velocity_magnitude'],'small vy')

#ds = car.load()
#ds.add_particle_filter('big_dust')
#ad=ds.all_data()
#print(ad['big_dust','particle_index'])

#creation_time dynamical_time metallicity_fraction typeia_fraction
