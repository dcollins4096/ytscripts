if 'quan_box' not in dir():
    ef('turb_quan.py')
    qb4 = fPickle.load('quan_box_eq44.pickle')
    qb4.stuff.keys()
    for n in qb4.stuff.keys():
        qb4.stuff[n] = nar(qb4.stuff[n])
#print qb4.stuff['ex']+qb4.stuff['ey']+qb4.stuff['ez'] - qb4.stuff['ke_tot']

if 'aw17' not in dir():
    aw17=taxi.taxi('aw17')
    aw17.frames=[100]
    aw17.region_type='all'
    reg = aw17.get_region()

if 0:
    qb = quan_box(aw17)
    y=qb(aw17)

if 0:
    for n in qb.stuff.keys():
        qb.stuff[n] = nar(qb.stuff[n])
    print qb.stuff['ex']+qb.stuff['ey']+qb.stuff['ez'] - qb.stuff['ke_tot']

if 1:
    dv = reg['cell_volume'].in_units('code_length**3')
    vx = reg['velocity_x'].in_units('code_velocity')
    vy = reg['velocity_y'].in_units('code_velocity')
    vz = reg['velocity_z'].in_units('code_velocity')
    rho= reg['density'].in_units('code_density')
    ke1 = 0.5*(rho*(vx*vx+vy*vy+vz*vz)*dv).sum()
    ke2 = (reg['kinetic_energy']*dv).in_units('code_mass*code_velocity**2').sum()
    print ke1
    print ke2


if 0:
    stuff={}
    all_fields = ['vx','vy','vz','mach','px','py','pz','ex','ey','ez','t','bx','by','bz','bx2','by2','bz2']
    all_fields +=['Bx','By','Bz','Bfield_strength','AlfMach','beta','AlfvenSpeed','frames']
    all_fields +=['ke_tot','ke_rel','grav_pot','grav_pot_2','gas_work']
    for k in all_fields:
        stuff[k]=[]
    ds=aw17.ds
    stuff['vx'].append(reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume').v)
    stuff['vy'].append(reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume').v)
    stuff['vz'].append(reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume').v)
    stuff['px'].append(reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume').v)
    stuff['py'].append(reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume').v)
    stuff['pz'].append(reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume').v)
    stuff['ex'].append(reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume').v)
    stuff['ey'].append(reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume').v)
    stuff['ez'].append(reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume').v)
    stuff['ke_tot'].append(reg.quantities['WeightedAverageQuantity']('kinetic_energy','cell_volume').v)
    reg.set_field_parameter('bulk_velocity',ds.arr([stuff['vx'][-1],stuff['vy'][-1], stuff['vz'][-1]],'code_velocity'))
    stuff['ke_rel'].append(reg.quantities['WeightedAverageQuantity']('rel_kinetic_energy','cell_volume').v)
    for n in stuff.keys():
        stuff[n] = nar(stuff[n])
    print stuff['ex']+stuff['ey']+stuff['ez'] - stuff['ke_tot']

if 0:
    dv = reg['cell_volume'].in_units('code_length**3')
    vx = reg['velocity_x'].in_units('code_velocity')
    vy = reg['velocity_y'].in_units('code_velocity')
    vz = reg['velocity_z'].in_units('code_velocity')
    rho= reg['density'].in_units('code_density')
    ke1 = 0.5*(rho*(vx*vx+vy*vy+vz*vz)*dv).sum()
    ke2 = (reg['kinetic_energy']*dv).in_units('code_mass*code_velocity**2').sum()
    print ke1
    print ke2

