execfile('go')
reload(taxi)
#execfile('clump_particles.py')
fltc=taxi.fleet(['aj23', 'aj24']) #,'aj22'])

all_fields=["Cooling_Time", "DebugField", "Density", "Electron_Density", "GasEnergy", "HII_Density", "HI_Density", "HeIII_Density", "HeII_Density", "HeI_Density", "Metal_Density", "Temperature", "TotalEnergy", "x-velocity", "y-velocity", "z-velocity"]
metal_fields=["Cooling_Time", "Electron_Density", "GasEnergy", "HII_Density", "HI_Density", "HeIII_Density", "HeII_Density", "HeI_Density", "Metal_Density", "Temperature", "TotalEnergy"]
def _dbg(field,data):
    return data['DebugField']
yt.add_field('SFdbg',function=_dbg,take_log=False)

fltc['fields'] = ['Electron_Density']
fltc['frames'] =  range(0,20,5)
fltc['Colorbar']='Fixed'
fltc['restrict']=False
#fltc['slice_zlim']={'Electron_Density':[0.5e-20,1.5e-20]}
fltc['slice_zlim']={'Electron_Density':[0.01,1]}
fltc('car.outname = car.name + "_hax_4"')
fltc.plot()
#level=1
#density=2
#divergence=3
#cooling time = 4
#jeans mass = 5
if 0:
    fltc['axis']=[2]
#fltc['weight_field']=None
    fltc['frames']=[6] #range(200)
    fltc['fields']='SFdbg'
    fltc['restrict'] = True
#fltc['callbacks'] = ['new_particles']
    fltc('car.callback_args["new_particles"]={"args":[1],"kwargs":{"col":"r"}}')
    fltc('car.outname = car.name+"_disk"')
    fltc['normal'] = [0,0,1]
    fltc['region_type']='disk'
    fltc['height'] = (1./32,'code_length')
    fltc.plot()
#party=[]
#for n in 0,1,2:
#    fltc[0].fill(n)
#    party.append(fltc[0].get_new_indices())

