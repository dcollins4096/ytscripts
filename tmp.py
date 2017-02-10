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

