#!/usr/bin/env python
#import yt
from e22_def_fields import *
import dsdiff_helpers

base_dir = '/scratch1/dcollins/EnzoProjects/E22_DEF_NonCosmo/'
thelist = ['c21_ct_live',    'c22_ded_live','c20_ppm_live', 'c23_ded_trunk']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 
            'c33_ct_live_field_nodef','c34_ded_live_field_nodef', 'c24_ded_live_def_field_riemann0']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 'c33_ct_live_field_nodef','c34_ded_live_field_nodef',
          'c24_ded_live_def_field_riemann0']
#thelist = ['c34_ded_live_field_nodef', 'c33_ct_live_field_nodef']
#thelist = ['f01_ppm_nodef', 'f02_ct_nodef', 'f04_ct_new', 'f03_ded_nodef']
thelist = ['f01_ppm_nodef',  'f03_ded_nodef', 'f05_ct_newversion_newparam','f01b_ppm_problem0','f01c_ppm_probM1']
thelist = ['f01_ppm_nodef',  'f03_ded_nodef', 'f05_ct_newversion_newparam']
thelist += ['f08_ppm_prgio']
#thelist += [ 'f06_ppm_hll']
#y_field = 'magnetic_energy'
x_field = 'density'
y_field = 'temperature'
x_log = False
y_log = False
prefix='RhoT_j'
for frame in [50]: #[1]+ range(10,60,10): #[2,3,100,200,268]:# + range(10,40,10):
    for dirname in thelist:
        if 1:
            setname = dsdiff_helpers.get_ds_name("%s/%s"%(base_dir,dirname),frame)
            ds =  yt.load(setname)
            if 0:
                y_min = 1e9
                y_max = 0
                x_min =1e9
                x_max = 0
            elif 0:
                y_min = (0.994,'K')
                y_max = (1.010,'K')
                x_min = 0.985
                x_max = 1.010
            else:
                y_min = (0.90,'K')
                y_max = (1.1,'K')
                x_min = 0.90
                x_max = 1.1
            print setname
            ad = ds.all_data()
            if 0:
                stat(ad[x_field], '%s %s'%(x_field,dirname))
                stat(ad[y_field], '%s %s'%(y_field,dirname))
                x_min = min(ad[x_field].min(), x_min)
                x_max = max(ad[x_field].max(), x_max)
                y_min = min(ad[y_field].min(), y_min)
                y_max = max(ad[y_field].max(), y_max)

    for dirname in thelist:
        setname = dsdiff_helpers.get_ds_name("%s/%s"%(base_dir,dirname),frame)
        print setname
        ds = yt.load(setname)
        #stat(ds.all_data()[x_field], 'rho %s'%dirname)
        #stat(ds.all_data()[y_field], 'T   %s'%dirname)
        ad = ds.all_data()
        if 1:
            if 1:
                print "=========", y_min, y_max
                phase = yt.create_profile(ad,bin_fields=[x_field,y_field], 
                                          fields=['cell_mass'],weight_field=None,
                                          extrema={x_field:[x_min, x_max], y_field:[y_min,y_max]}, 
                                          n_bins=[128,128], 
                                          logs={x_field:x_log,y_field:y_log}) #, n_bins=[32,32])
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
