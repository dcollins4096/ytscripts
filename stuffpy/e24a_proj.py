
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

for ng,grid in enumerate(ds.index.grids):
    if ng != 44:
        continue
    
    if grid.Level == 2:
        stat(np.abs(grid['x-velocity']),grid)
        g = grid['x-velocity']
        for k in range( g.shape[2]):
            plt.clf()
            plt.imshow( g[:,:,k].v, interpolation  = 'nearest', origin='lower')
            stat(g[:,:,k].v,k)
            outname = 'gslice2_g%04d_k%04d'%(ng,k)
            plt.savefig(outname)


