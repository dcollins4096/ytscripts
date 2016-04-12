#!/usr/bin/env python

#LightweigtPlotter.py frame name magnetic vectors

import yt
import sys
import dsdiff_helpers


basedir = '/scratch1/dcollins/EnzoProjects/E22_DEF_NonCosmo/'
dirname = 'f02_ct_nodef'

frame = 266
name = 'f02'
#if len(sys.argv) < 3:
#    print "plotter.py frame name"
#    sys.exit(0)
#frame = int( sys.argv[1])
#name = sys.argv[2]

field_list  = ['density','temperature']#,'metallicity']
do_the_vectors=False
if len(sys.argv) >= 4:
    print "her"
    field_list += ['magnetic_energy']
if len(sys.argv) >= 5:
    do_the_vectors=True
ds =  yt.load(dsdiff_helpers.get_ds_name("%s/%s"%(basedir,dirname),frame))
for field in field_list:
    print field
    try:
        proj = ds.proj(field,2)
        pw = proj.to_pw() #width=(0.1,'code_length'))
        if do_the_vectors:
            pw.annotate_streamlines('Bx','By')
        print pw.save('%s_zoom_%04d'%(name, frame))
    except:
        print "Problem with Field", field
    

#end

#end
