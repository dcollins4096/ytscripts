#!/usr/bin/env python
#import yt
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
