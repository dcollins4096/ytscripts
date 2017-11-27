if 'ef' not in dir():
    execfile('go')
reload(taxi)
car=taxi.taxi('g18e_tight')
car.plot()
if 0:
    car.fields=['density'] #,'magnetic_field_strength']
    car.axis=[2]
    car.plot()
    L= car.left
    R= car.right
    C= car.center
    print L.in_units('kpc')
    print R.in_units('kpc')
    print C.in_units('kpc')
    car.right[2] = car.center[2] + car.ds.arr(3,'kpc')
    car.left[2] = car.center[2] - car.ds.arr(3,'kpc')
    car.callbacks += ['pol_frac']
    car.outname = 'g30_'+car.name+'_pfrac'
    print car.outname
    car.frames=[90]
    car.plot()
