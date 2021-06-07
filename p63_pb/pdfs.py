from go import *
reload(taxi)

if 'car' not in dir():
    #car = taxi.load('b02')
    car = taxi.load('ca02')

def costhetaz(field,data):
    return data['magnetic_field_z']/data['magnetic_field_strength']
def thetaz(field,data):
    output = np.zeros_like( data['costhetaz'])
    output[:] = np.arccos(data['costhetaz'][:].v)
    return output


def add_cos_theta(obj):
    obj.add_field('costhetaz',costhetaz,units='dimensionless')
    obj.add_field('thetaz',thetaz,units='dimensionless')

#car.derived_fields['add_cos_theta']=add_cos_theta
##car.load(100)
##car.ds.all_data()['density']
#g=car.ds.index.grids[0]

if 1:
    #car.fields=['costhetaz']
    car.frames=[100]
    #car.plot()
    ybins=np.arange(-2,2,0.1)
    xbins=np.arange(1e-3,100)
    bins={'costhetaz':ybins,'magnetic_field_strength':xbins}
    #car.phase(['magnetic_field_strength','costhetaz','cell_volume'],phase_args={'override_bins':bins})
    #car.profile(['costhetaz','cell_volume'],scales=['linear','linear'])
    car.profile(['thetaz','cell_volume'],scales=['linear','linear'])

if 0:
    plt.clf()
    prof = car.last_prof
    bin_edges = prof.x_bins
    bins =0.5*(bin_edges[1:]+bin_edges[:-1])
    plt.plot(bins, prof['cell_volume'])
    plt.savefig('../PigPen/test.png')


