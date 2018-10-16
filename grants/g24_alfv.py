execfile('go')
def alf(field,data):
    b=data['magnetic_field_strength'].in_units('gauss')
    rho=(4*np.pi*data['density'].in_units('g/cm**3'))**0.5
    v = data['velocity_magnitude'].in_units('cm/s')
    return np.abs(b.v/rho.v) + np.abs(v.v)

yt.add_field('alfv',function=alf, units='')


fname = '/scratch1/dcollins/Paper08/B02/512/RS0080/restart0080'
#pos = [0.570190429688,0.720092773438, 0.987670898438] g1306 for B02/256/rs0080
ds = yt.load(fname)
max_alf = 0.0
for g in ds.index.grids:
    if max_alf < g['alfv'].max():
        max_alf = g['alfv'].max()
        ng = g
    print g, "%0.2e"%max_alf


print ds.index.grids[-1]['alfv']

