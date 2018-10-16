if 'ef' not in dir():
    execfile('go')
import tube
reload(tube)


sim1 = '/Users/dcollins/Enzo/dave_vs_enzo_COSMOLOGY_UP/run/Cosmology/MHDZeldovichPancake'
sim2 = '/Users/dcollins/Enzo/dave_vs_enzo_COSMOLOGY_UP/run/Cosmology/MHDCTZeldovichPancake'
coord=(0.505,0.505)
frame_list = [0]
#ds_list = [yt.load("%s/DD%04d/data%04d"%(sim,frame,frame)) for frame in frame_list]
#rays = [ds.ortho_ray(0,coord) for ds in ds_list]
for this_frame in [0,20]:
    frame_list = [this_frame]
    ds_list = [yt.load("%s/DD%04d/data%04d"%(sim1,frame,frame)) for frame in frame_list]
    ds_list += [yt.load("%s/DD%04d/mhdct%04d"%(sim2,frame,frame)) for frame in frame_list]
    outname = 'with_a_in_pressure_%d.pdf'%frame
    fields = ['density','pressure','By','x-velocity']
    y=tube.tube(ds_list,   fields=fields, filename = outname,legend=True) #, renorm=renorm,filename = "sine"+".png")
if 0:
    print outname
#print rays[0]['magnetic_energy']/rays[2]['magnetic_energy']
    print rays[0]['magnetic_energy']/rays[2]['magnetic_energy']
    print "print rays[1]['magnetic_energy']/rays[3]['magnetic_energy']"
    x= rays[0]['density']/rays[2]['density']
    print "densty 0", x.min(), x.max()
