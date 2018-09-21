basedir = '/scratch1/dcollins/EnzoProjects/E24_JonathanTanProjects/E24a_Ben/'
aa01 = ''
aa02 = ''
aa03 = 'aa03_l02'
aa04 = 'aa04_fixed_proj'
sim=aa04
frame = 0
setname='test_4L64_init_interp4_'
dsname = "%s/%s/DD%04d/%s%04d"%(basedir,sim,frame,setname,frame)
ds = yt.load(dsname)
for n in range(64):
    #field = 'grid_indices'
    field = 'x-velocity'
    sl=yt.SlicePlot(ds,2,field,center=[0.5,0.5,(n+0.5)/64.+np.pi/2048])
    #sl.set_zlim('grid_indices',160,170)
    sl.annotate_grids(draw_ids=True)
    print sl.save('%s_sl%04d_n%04d'%(sim+"Offset",n,frame))
