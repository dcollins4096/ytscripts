#!/usr/bin/env python
import yt
from e22_def_fields import *
execfile('/Users/dcollins/yt3_scripts/go_lite')
execfile('/Users/dcollins/yt3_scripts/davetools.py')
thelist = ['c21_ct_live',    'c22_ded_live','c20_ppm_live', 'c23_ded_trunk']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 
            'c33_ct_live_field_nodef','c34_ded_live_field_nodef', 'c24_ded_live_def_field_riemann0']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 'c33_ct_live_field_nodef','c34_ded_live_field_nodef',
          'c24_ded_live_def_field_riemann0']
#thelist = ['c34_ded_live_field_nodef', 'c33_ct_live_field_nodef']

#thelist = ['c01_ppm'] #,'c60_ppm_allthings','c61_ct_allthings', 'c02_ct']
thelist = ['c69_hydro3_rad12','c66_ded_bsmall_Rad12','c67_ct2_bsmall_rad12','c68_ppm_rad12'] #RadiationFieldType = 12
thelist = ['c71_ct_new_bsmall_rad12','c67_ct2_bsmall_rad12','c66_ded_bsmall_Rad12','c72_ct_new_bsmall_rad12_PRtest']

#thelist += ['c61_ct_allthings','c62_ded_allthings','c69_hydro3_rad12']

#thelist = ['c01b_ppm_hll','c02b_ct_hll']
#y_field = 'magnetic_energy'
x_field = 'density'
y_field = 'temperature'

#y_field = 'HI_Density'
basename = '/Users/dcollins/scratch/EnzoProjects/E12/E12b_Cosmology/'
prefix='PR_prior'
for frame in [9]:# + range(10,40,10):
    y_min = 1e9
    y_max = 0
    x_min =1e9
    x_max = 0
    for dirname in thelist:
        setname = "%s/%s/RD%04d/RD%04d"%(basename,dirname,frame,frame)
        #setname = "%s/%s/DD%04d/DD%04d"%(basename,dirname,frame,frame)
        print setname
        ds = yt.load(setname)
        ad = ds.all_data()
        stat(ad[x_field], '%s %s'%(x_field,dirname))
        stat(ad[y_field], '%s %s'%(y_field,dirname))
        x_min = min(ad[x_field].min(), x_min)
        x_max = max(ad[x_field].max(), x_max)
        y_min = min(ad[y_field].min(), y_min)
        y_max = max(ad[y_field].max(), y_max)
        #y_min = 14
        #y_max = 30

    for dirname in thelist:

        #setname = "%s/%s/DD%04d/DD%04d"%(basename,dirname,frame,frame)
        setname = "%s/%s/RD%04d/RD%04d"%(basename,dirname,frame,frame)
        ds = yt.load(setname)
        #stat(ds.all_data()[x_field], 'rho %s'%dirname)
        #stat(ds.all_data()[y_field], 'T   %s'%dirname)
        ad = ds.all_data()
        if 1:
            if 1:
                print "=========", y_min, y_max
                phase = yt.create_profile(ad,bin_fields=[x_field,y_field], 
                                          fields=['cell_mass'],weight_field=None,
                                          extrema={x_field:[x_min, x_max], y_field:[y_min,y_max]}, n_bins=[128,128], 
                                          logs={x_field:True,y_field:True}) #, n_bins=[32,32])
                pp = yt.PhasePlot.from_profile(phase)
            else:
                pp = yt.PhasePlot(ad, x_field, y_field, 'cell_mass',weight_field=None)
            pp.set_xlabel(x_field)
            pp.set_ylabel(y_field)
            print pp.save("%s_%s_DD%04d"%(prefix,dirname[0:6],frame))
        if 0:
            prof = yt.ProfilePlot(ad,y_field,'cell_mass',weight_field=None)
            print prof.save("%s_%s_RD%04d"%(prefix,dirname[0:6],frame))

#        den = ad['density']
#        T = ad['Temperature']
        
#   phase = yt.PhasePlot(ds.all_data(),'density','magnetic_energy','cell_mass',weight_field=None)
#   phase.save("%s_%04d"%(dirname,frame))
#end
