from go import *
import p49_fields
reload(taxi)
if 'flt' not in dir():
    #pg00 = taxi.load('pg00')
    #this_car=pg00
    #ds=this_car.load()
    flt = taxi.fleet(['pg%d'%s for s in [3,4,5,6,7,8,9,10]])

    flt['frames']=[100]#range(10,100,10)
#this_car.fields=['density']
    flt['fields']=['magnetic_field_strength','density','velocity_magnitude']
    flt.plot()
if 1:
    import turb_quan
    reload(turb_quan)
    for car in flt:
        car.derived_fields['QU'] = p49_fields.add_QU
        qb = turb_quan.quan_box(car=car)
        qb.EBall()
        #qb()


#acceleration_x-512-C.hdf5  acceleration_z-512-C.hdf5     ce'pg1','pg2'll_centered_B_y-512-C.hdf5  density-512-C.hdf5   velocity_x-512-C.hdf5  velocity_z-512-C.hdf5
#acceleration_y-512-C.hdf5  cell_centered_B_x-512-C.hdf5  cell_centered_B_z-512-C.hdf5  pressure-512-C.hdf5  velocity_y-512-C.hdf5
