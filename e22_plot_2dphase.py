#!/usr/bin/env python
import yt
from def_fields import *
execfile('/Users/dcollins/yt3_scripts/go_lite')
execfile('/Users/dcollins/yt3_scripts/davetools.py')
thelist = ['c21_ct_live',    'c22_ded_live','c20_ppm_live', 'c23_ded_trunk']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 
            'c33_ct_live_field_nodef','c34_ded_live_field_nodef', 'c24_ded_live_def_field_riemann0']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 'c33_ct_live_field_nodef','c34_ded_live_field_nodef',
          'c24_ded_live_def_field_riemann0']
#thelist = ['c34_ded_live_field_nodef', 'c33_ct_live_field_nodef']
thelist = ['x01_def_unigrid','x02_nodef_unigrid']
#y_field = 'magnetic_energy'
x_field = 'density'
y_field = 's_dave'
prefix='entropy_2d'
for frame in [1,2,3,4,5]:# + range(10,40,10):
    y_min = 1e9
    y_max = 0
    x_min =1e9
    x_max = 0
    for dirname in thelist:
        setname = "%s/DD%04d/DD%04d"%(dirname,frame,frame)
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

        ds = yt.load("%s/DD%04d/DD%04d"%(dirname,frame,frame))
        #stat(ds.all_data()[x_field], 'rho %s'%dirname)
        #stat(ds.all_data()[y_field], 'T   %s'%dirname)
        ad = ds.all_data()
        if 1:
            if 1:
                print "=========", y_min, y_max
                phase = yt.create_profile(ad,bin_fields=[x_field,y_field], 
                                          fields=['cell_mass'],weight_field=None,
                                          extrema={x_field:[x_min, x_max], y_field:[y_min,y_max]}, n_bins=[128,128], 
                                          logs={x_field:True,y_field:False}) #, n_bins=[32,32])
                pp = yt.PhasePlot.from_profile(phase)
            else:
                pp = yt.PhasePlot(ad, x_field, y_field, 'cell_mass',weight_field=None)
            pp.set_xlabel(x_field)
            pp.set_ylabel(y_field)
            print pp.save("%s_%s_%04d"%(prefix,dirname[0:6],frame))
        if 0:
            prof = yt.ProfilePlot(ad,y_field,'cell_mass',weight_field=None)
            print prof.save("%s_%s_%04d"%(prefix,dirname[0:6],frame))

#        den = ad['density']
#        T = ad['Temperature']
        
#   phase = yt.PhasePlot(ds.all_data(),'density','magnetic_energy','cell_mass',weight_field=None)
#   phase.save("%s_%04d"%(dirname,frame))
#end
