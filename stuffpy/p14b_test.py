if 'ef' not in dir():
    execfile('go')
import pyfits
ef('p14_clump_mask.py')


dataset = '/scratch1/dcollins/Paper08/B02/256/RS0050/restart0050'

if 'ds' not in dir():
    ds = yt.load(dataset)
    
if 'sphere' not in dir() and 0:
    ds = yt.load(dataset)
    max_density, position = ds.h.find_max('Density')
    sphere = ds.h.sphere(position, 0.05)


if 0:
    proj = yt.ProjectionPlot(ds,0,'Density')
    frb=proj._frb
    plt.clf()
    plt.imshow(np.log10(frb['Density']))
    plt.colorbar()
    plt.savefig('p14_atest.png')

if 0:
    field = 'Density'
    fits_name = 'p14b_atest.fits'
    hdu = pyfits.PrimaryHDU(frb[field])
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(fits_name)

clump_fname = 'p14_atest1_Mask.fits'
if 0:
    import clfind2d
    import pyfits
    if os.path.exists(clump_fname):
        os.remove(clump_fname)
    clfind_args = {'levels':[3], 'log':True,'nPixMin':3}
    clfind2d.clfind2d(fits_name,'p14_atest1',**clfind_args)

if 0:
    #twod_clump = np.transpose(pyfits.open(clump_fname)[0].data)
    twod_clump = pyfits.open(clump_fname)[0].data
    flat_clump = twod_clump.flatten()
    clump_axis = 0
    plt.clf()
    plt.imshow(twod_clump)
    plt.savefig('p14_atest_clump.png')

if 0:
    clump_axis = 0
    twod_clump = np.transpose(pyfits.open(clump_fname)[0].data)
    flat_clump = twod_clump.flatten()
    clump_stuff_tuple = (flat_clump,'xyz'[clump_axis],twod_clump)
    clump_stuff_tuple = (flat_clump,clump_axis,twod_clump)
    clump_stuff_dict={'clump_mask_stuff':clump_stuff_tuple}
    clump_ind = 5
    ef('p14_clump_mask.py')
    ef('p14_get_cut_region.py')

if 0:
    proj = ProjectionPlot(ds,0,'Density')
    proj.annotate_clumps([cut_region])
    print proj.save('p14b_atest_cut')

if 0 and False:
    #dave clumps doesn't really work.
    ef('p14_clumps.py')

if 0:
    master_clump = Clump(cut_region,"density")
    master_clump.add_validator("min_cells", 27)
    c_min=10
    c_max=cut_region['density'].max()
    step = 1.1
    find_clumps(master_clump, c_min, c_max, step)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.

if 1:
    import clump_properties
    reload(clump_properties)
    clump_stuff_list = []
    for clump in leaf_clumps:
        clump_stuff_list.append(clump_properties.clump_stuff(clump,ds,1))
    

if 1:
    my_subset = clump_subset.subset(clump_stuff_list)
    


#frb
#clumps from frb
#cut region from clumps
#3dclumps
#alpha on 3d clumps.  I don't need much else.
