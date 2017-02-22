#!/usr/bin/env python
#import yt


execfile('go')
reload(taxi)
execfile('clump_particles.py')
if 'fltc' not in dir():
    #fltc=taxi.fleet(['aj25', 'aj26']) #,'aj22'])
    #fltc=taxi.fleet(['aj29','aj30','aj31','aj32']) #,'aj22'])
    fltc=taxi.fleet(['aj37','aj38','aj39'])
    fltc=taxi.fleet(['aj31','aj43','aj41','aj45'])
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

if 1:
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
    for f2 in metalfields[1:]: #['metal_accounting_2']: #
        for frame in [10,20,30]: #range(30,500,15):
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
# old version 
if 0:
    ef('p33_sims.py')
    sim_list = ['i02','i04']
    field2 = 'Cooling_Time'
    for frame in [44]: #[4]: #range(4):
        T_min = 1e9
        T_max = 0
        den_min =1e9
        den_max = 0
        for sim in sim_list:
            name  = '%s/%s/DD%04d/DD%04d'%(bd,sim_base_dir[sim],frame,frame)
            ds = yt.load(name)
            stat(ds.all_data()['density'], 'rho %s'%sim)
            stat(ds.all_data()[field2], 'T   %s'%sim)
            ad = ds.all_data()
            den_min = min(ad['density'].min(), den_min)
            den_max = max(ad['density'].max(), den_max)
            T_min =  min(ad[field2].min(), T_min)
            #T_min = max(1e-15,T_min)
            T_max = max(ad[field2].max(), T_max)

        for sim in sim_list:
            name  = '%s/%s/DD%04d/DD%04d'%(bd,sim_base_dir[sim],frame,frame)
            ds = yt.load(name)
            stat(ds.all_data()['density'], 'rho %s'%sim)
            stat(ds.all_data()[field2], 'T   %s'%sim)
            ad = ds.all_data()
            den = ad['density']
            T = ad[field2]
            phase = yt.create_profile(ds.all_data(),bin_fields=['density',field2], fields=['cell_mass'],weight_field=None,
                                 extrema={'density':[den_min, den_max], field2:[T_min,T_max]}, n_bins=[32,32])
            pp = yt.PhasePlot.from_profile(phase)
            pp.set_xlabel('density')
            pp.set_ylabel(field2)
            print pp.save("%s_%04d"%(sim,frame))
#   phase = yt.PhasePlot(ds.all_data(),'density','magnetic_energy','cell_mass',weight_field=None)
#   phase.save("%s_%04d"%(dirname,frame))
#end
