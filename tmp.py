<<<<<<< /home1/00369/tg456484/yt3_scripts/tmp.py
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
=======
execfile('go')
reload(taxi)
execfile('clump_particles.py')
if 'fltc' not in dir():
    #fltc=taxi.fleet(['aj25', 'aj26']) #,'aj22'])
    #fltc=taxi.fleet(['aj29','aj30','aj31','aj32']) #,'aj22'])
    fltc=taxi.fleet(['aj37','aj38','aj39'])
    fltc['frames'] = [500]# range(0,500,15) + [500]
    fltc['callbacks']=['particles','nparticles']
    fltc['fields'] = ['density','temperature']
    fltc['operation'] = 'RegionProjection'

if 0:
    def _dbg(field,data):
        return data['DebugField']
    yt.add_field('SFdbg',function=_dbg,take_log=False)

    ef('p33_sfhunt.py')

    fltc['frames']=[18]
    fltc['fields'] = ['mjeans','bmass']
    fltc('car.outname = car.name+"_t2"')
    fltc['callbacks']=[]

#fltc['region_type'] = 'grid-1'
#a = fltc[0].get_region()
#fltc.find_extrema(fields=['bmass','mjeans'])
#fltc('print car.extrema["bmass"]')
#fltc('print car.extrema["mjeans"]')
#fltc.find_extrema(['bmass','mjeans'])
#fltc.phase(['bmass','mjeans','cell_volume'])

if 0:
    fltc['frames'] = [10,100,500]
    fltc.profile(['metal_accounting','cell_volume'],scales=['linear','linear'])

if 0:
    metalfields = [
    "Cooling_Time",
    "Density",
    "Electron_Density",
    "HII_Density",
    "HI_Density",
    "HeIII_Density",
    "HeII_Density",
    "HeI_Density",
    "Metal_Density",
    "Temperature",
    "TotalEnergy",
    'kinetic_energy']
    fields=['density','jeans_density']
    fltc['frames'] = [500]
    fltc.find_extrema(metalfields, manual_positive=True)
    fltc['frames'] = range(15,500,15) 
#fltc.phase(['density','jeans_density','cell_volume'])
if 1:
    for f2 in ['metal_accounting_2']: #metalfields[1:]:
        for frame in [10,100,500]: #range(30,500,15):
            fltc['frames']=[frame]
            fltc.phase(['density',f2,'cell_volume'])
    #fltc.plot()
#fltc.plot()

#level=1
#density=2
#divergence=3
#cooling time = 4
#jeans mass = 5
#fltc['axis']=[2]
#fltc['weight_field']=None
#fltc['frames']=range(10,110,10)
#fltc['fields']=['density','SFdbg']
#fltc['callbacks'] = ['nparticles','particles']# ['new_particles']
#fltc('car.callback_args["new_particles"]={"args":[1],"kwargs":{"col":"r"}}')
#fltc['operation']='CenterSlice'
#fltc('car.outname = car.name+"_disk"')
#fltc['normal'] = [0,0,1]
#fltc['region_type']='disk'
#fltc['height'] = (1./32,'code_length')
#party=[]
#for n in 0,1,2:
#    fltc[0].fill(n)
#    party.append(fltc[0].get_new_indices())
>>>>>>> /tmp/tmp.py~other.bQP_yc

