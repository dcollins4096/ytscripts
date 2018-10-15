#!/usr/bin/env python

execfile('go')
#yt.enable_parallelism()
if 0:
    if len(sys.argv) < 3:
        print "plotter.py frame name"
        sys.exit(0)
    frame = int( sys.argv[1])
    name = sys.argv[2]

field_list  = ['density','temperature']#,'metallicity']
eq42 = taxi.taxi(directory='/scratch/00369/tg456484/Paper42_NewAK/eq42_m9_grav_512_p59_ppm_L4_J32',name='eq42')

car=eq42
car.fill(0)
car.frames=car.frame_dict.keys()
car.fields=['density']
car.Colorbar='monotonic'
car.callbacks=['grids']
car.plot()

if 0:
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
