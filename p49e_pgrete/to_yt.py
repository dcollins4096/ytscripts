import p49_fields
reload(taxi)
if 'pg00' not in dir():
    pg00 = taxi.load('pg00')
    this_car=pg00
    ds=this_car.load()

if 1:
    import turb_quan
    reload(turb_quan)
    for car in [this_car]:
        car.derived_fields['QU'] = p49_fields.add_QU
        qb = turb_quan.quan_box(car=car)
        qb.EBall()
        qb()


#acceleration_x-512-C.hdf5  acceleration_z-512-C.hdf5     cell_centered_B_y-512-C.hdf5  density-512-C.hdf5   velocity_x-512-C.hdf5  velocity_z-512-C.hdf5
#acceleration_y-512-C.hdf5  cell_centered_B_x-512-C.hdf5  cell_centered_B_z-512-C.hdf5  pressure-512-C.hdf5  velocity_y-512-C.hdf5
