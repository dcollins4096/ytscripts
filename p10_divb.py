

output_template = '/data/astro12_3/student3/runs/%d/DD%04d/data%04d'
def stat(array,strin='', format='%0.16e'):
    template = '['+format+','+format+'] %s %s'
    print template%(array.min(),array.max(),array.shape,strin)

all_xbins = []
all_profiles = []
fields = ['Qplus','cell_volume']
for sim in [192]:
    for frame in [1000]:
        fname = output_template%(sim,frame,frame)
        ds = yt.load(fname)
        for g in ds.index.grids[0:1]:
            h5ptr = h5py.File(g.filename,'r')
            grid = h5ptr['Grid%08d'%g.id]
            BxF = grid['BxF'][:]
            ByF = grid['ByF'][:]
            BzF = grid['BzF'][:]
            h5ptr.close()
            #hdf5 rotates the indices, so if you take the data off disk, x and z axes are swapped.
            DeltaB = ((BxF[:,:,:-1]-BxF[:,:,1:])+   (ByF[:,:-1,:]-ByF[:,1:,:])+   (BzF[:-1,:,:]-BzF[1:,:,:]))
            B =  0.5*((BxF[:,:,:-1]+BxF[:,:,1:])**2+(ByF[:,:-1,:]+ByF[:,1:,:]**2)+(BzF[:-1,:,:]+BzF[1:,:,:])**2)**0.5
            db_b = np.abs(DeltaB/B)
            stat(db_b)
            

