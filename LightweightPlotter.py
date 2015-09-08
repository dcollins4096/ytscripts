#!/usr/bin/env python

import yt
import sys
if len(sys.argv) < 3:
    print "plotter.py frame name"
    sys.exit(0)
frame = int( sys.argv[1])
name = sys.argv[2]

field_list  = ['density','temperature']#,'metallicity']
do_the_vectors=False
if len(sys.argv) == 4:
    print "her"
    field_list += ['magnetic_energy']
    do_the_vectors=True
for field in field_list:
    print field
    ds =  yt.load('DD%04d/DD%04d'%(frame,frame))
    proj = ds.proj(field,2)
    pw = proj.to_pw() #width=(0.1,'code_length'))
    if do_the_vectors:
        pw.annotate_streamlines('Bx','By')
    pw.save('%s_zoom_%04d'%(name, frame))
    

#end

#end
