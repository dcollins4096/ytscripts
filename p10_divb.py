

output_template = '/data/astro12_3/student3/runs/%d/DD%04d/data%04d'
def stat(array,strin='', format='%0.16e'):
    template = '['+format+','+format+'] %s %s'
    print template%(array.min(),array.max(),array.shape,strin)

all_xbins = []
all_profiles = []
if 'car' not in dir():
    car = taxi.taxi('aw21')
    car.fill(496)

def full_div_b(ds,grid_slice=slice(None), do=['extrema']):
    output = {}
    for n, g in enumerate(ds.index.grids[grid_slice]):
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
        if 'stat_all' in do:
            stat(db_b)
        if 'extrema' in do:
            if n == 0:
                min_divb = db_b.min()
                max_divb = db_b.max()
            else:
                min_divb = min([min_divb,db_b.min()])
                max_divb = max([max_divb,db_b.max()])
            output['extrema'] = [min_divb,max_divb]
    return output

ex= full_div_b(car.ds, do=['extrema','stat_all'])
print ex


                

