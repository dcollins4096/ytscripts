execfile('go')
reload(taxi)
execfile('clump_particles.py')
fltc=taxi.fleet(['aj23', 'aj24']) #,'aj22'])

def _dbg(field,data):
    return data['DebugField']
yt.add_field('SFdbg',function=_dbg,take_log=False)

#level=1
#density=2
#divergence=3
#cooling time = 4
#jeans mass = 5
fltc['axis']=[2]
#fltc['weight_field']=None
fltc['frames']=range(200)
fltc['fields']='SFdbg'
fltc['callbacks'] = ['new_particles']
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

