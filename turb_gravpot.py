

if 'eq44_grav' not in dir():
    eq44_grav=taxi.taxi('eq44_grav')

level=0
eq44_grav.fill(3)
dims=[512]*3
field_list = ['density','PotentialField']
cg = eq44_grav.ds.covering_grid(level,left_edge=[0.0,0.0,0.0],dims=dims,fields=field_list)
