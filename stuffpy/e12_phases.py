#!/usr/bin/env python
import yt
execfile('/Users/dcollins/yt3_scripts/davetools.py')
thelist = ['c21_ct_live',    'c22_ded_live','c20_ppm_live', 'c23_ded_trunk']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 'c33_ct_live_field_nodef','c34_ded_live_field_nodef',
          'c24_ded_live_def_field_riemann0']
for frame in range(0,40,10):
    T_min = 1e9
    T_max = 0
    den_min =1e9
    den_max = 0
    for dirname in thelist:
        ds = yt.load("%s/DD%04d/DD%04d"%(dirname,frame,frame))
        stat(ds.all_data()['density'], 'rho %s'%dirname)
        stat(ds.all_data()['Temperature'], 'T   %s'%dirname)
        ad = ds.all_data()
        den_min = min(ad['density'].min(), den_min)
        den_max = max(ad['density'].max(), den_max)
        T_min = min(ad['Temperature'].min(), T_min)
        T_max = max(ad['Temperature'].max(), T_max)

    for dirname in ['c20_ppm_live', 'c30_ppm_nodef']: #thelist:

        ds = yt.load("%s/DD%04d/DD%04d"%(dirname,frame,frame))
        stat(ds.all_data()['density'], 'rho %s'%dirname)
        stat(ds.all_data()['Temperature'], 'T   %s'%dirname)
        ad = ds.all_data()
        den = ad['density']
        T = ad['Temperature']
        phase = yt.create_profile(ds.all_data(),bin_fields=['density','Temperature'], fields=['cell_mass'],weight_field=None,
                             extrema={'density':[den_min, den_max], 'Temperature':[T_min,T_max]}, n_bins=[32,32])
        pp = yt.PhasePlot.from_profile(phase)
        pp.set_xlabel('density')
        pp.set_ylabel('Temperature')
        print pp.save("%s_%04d"%(dirname[0:6],frame))
#   phase = yt.PhasePlot(ds.all_data(),'density','magnetic_energy','cell_mass',weight_field=None)
#   phase.save("%s_%04d"%(dirname,frame))
#end
