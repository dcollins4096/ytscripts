import pyfits
import mw_stuff
import uuid
base = '/scratch1/dcollins/Paper08/B2/512'; outname = 'Mid_512'
#base = '/scratch1/dcollins/Paper08/B02/512'; outname = 'High_512'
if 'frame' not in dir():
    frame = 50
dsname = '%s/RS%04d/restart%04d'%(base,frame,frame)
def clump_list_sort(clump_list):
    """Returns a copy of clump_list, sorted by ascending minimum density.  This
    eliminates overlap when passing to
    yt.visualization.plot_modification.ClumpContourCallback"""
    minDensity = [c['density'].min() for c in clump_list]
    args = np.argsort(minDensity)
    list = nar(clump_list)[args]
    reverse = range(list.size-1,-1,-1)
    return list[reverse]

if 0:
    #make the covering grid and pseudo dataset
    size1 = 256
    nx,ny,nz = (256,256,256) # domain dimensions
    dims= [nx,ny,nz] #this could be cleaned up.
    ax = 0
    ds1 = yt.load(dsname)
    g_code = ds1.quan(ds1['GravitationalConstant']/(4.*np.pi), 'code_length*code_velocity**2/code_mass')
    cg = ds1.covering_grid(0,[0.0]*3, [size1]*3, fields=['velocity_x','density']) #['x-velocity','y-velocity','z-velocity','Density'])
    mindx = cg['dx'].min().in_units('code_length').v
    mindx = 1./256
    data = {}
    data["density"] =     (cg['density'],"code_mass/code_length**3")
    data["velocity_x"] =  (cg['velocity_x'], "code_velocity")
    data["velocity_y"] =  (cg['velocity_y'], "code_velocity")
    data["velocity_z"] =  (cg['velocity_z'], "code_velocity") # zero velocity in the z-direction
    bbox = np.array([[0.0,1.0],[0.0,1.0],[0.0,1.0]]) # bbox of width 1 on a side with center (0,0,0)
    ds = yt.load_uniform_grid(data, (nx,ny,nz), length_unit=(1,"code_length"), nprocs=1, bbox=bbox)
    ad = ds.all_data()

if 0:
    #Make clumps
    master_clump = Clump(ad,'density')
    master_clump.add_validator('min_cells',8)
    c_min = 10
    c_max = cg['density'].max()
    find_clumps(master_clump, c_min, c_max, 100)
    leaf_clumps = nar(get_lowest_clumps(master_clump)) #if both min_cells and grav_bound are used, this is empty.


if 0:
    g_code = ds1.quan(ds1['GravitationalConstant']/(4.*np.pi), 'code_length*code_velocity**2/code_mass')
    #g_code = ds1['GravitationalConstant']/(4.*np.pi)
    virial = np.zeros(len(leaf_clumps),dtype='float')
    for i, c in enumerate(leaf_clumps):
        ke = (c['kinetic_energy']*c['cell_volume']).sum().in_units('code_mass*code_velocity**2')
        R = (c['cell_volume'].sum())**(1./3)
        m = c['cell_mass'].sum().in_units('code_mass')
        #pretty sure this will work.
        #mw_temp_name = "./MW_Grav_Temp_%s"%(uuid.uuid1())
        #ge = -g_code*mw_stuff.run_mw_grav(mw_temp_name,data)
        binding = 3./5.*g_code*m**2/R
        virial[i] = ke/binding


axis = 0
axis1 = axis
axis2 = 1
axis3 = 2
x1 = 'x'
x2 = 'y'
x3 = 'z'
left = [0,0,0]
if 0:
    for ax in [0,1,2]:
        plot=yt.ProjectionPlot(ds,ax,'density')
        plot.set_cmap('density','Greys')
        #plot.annotate_clumps( clump_list_sort( leaf_clumps[ virial < 2] ) )
        plot.annotate_clumps( clump_list_sort( leaf_clumps ) )
        plot.save('p29_all')


if 0:
    #column density test image
    plt.clf()
    plt.imshow( np.log10(np.sum( cg['density'].v, axis=axis)), cmap='Greys' )
    plt.savefig('p29_column.png')




if 1:
    #full column density 
    plt.clf()
    colden = np.log10( np.sum(cg['density'].v, axis=axis)).flatten()
    vals,bins,obj=plt.hist( colden, histtype='step', bins=100, color='k')
    plt.yscale('log')
    plt.savefig('p29_hist.pdf')

if 1:
    spare_den = copy.copy( cg['density'].in_units('code_density').v) #.flatten()
    full_ind = np.arange(spare_den.size,dtype='int')
    for cl in leaf_clumps:
        #Index the 2d image plane for each pixel in the clump.
        i = list(((cl[x1].in_units('code_length').v- left[axis1])/mindx-0.5).astype('int'))
        j = list(((cl[x2].in_units('code_length').v- left[axis2])/mindx-0.5).astype('int'))
        k = list(((cl[x3].in_units('code_length').v- left[axis3])/mindx-0.5).astype('int'))
        #index = k + nz*(j+ ny*i)
        spare_den[i,j,k] = 0

    plt.clf()
    plt.imshow( np.log10( np.sum(spare_den,axis=0)), cmap='Greys')
    outname = 'p29_image_nobound.png'
    plt.savefig(outname)
    print outname




if 0:
    #
    # Make the column for each clump.
    # Does't really work.
    # 
    axis1 = axis
    axis2 = 1
    axis3 = 2
    x1 = 'x'
    x2 = 'y'
    x3 = 'z'
    left = [0,0,0]
    bin_dx = bins[1]-bins[0]
    collector = np.zeros(dims[axis2]*dims[axis3])
    for cl in leaf_clumps:
        #Index the 2d image plane for each pixel in the clump.
        j = (cl[x2].in_units('code_length').v- left[axis2])/mindx-0.5
        k = (cl[x3].in_units('code_length').v- left[axis3])/mindx-0.5
        index = j + dims[axis2]*k
        iunique = np.unique(index).astype('int')
        flat_col = np.zeros(iunique.shape)
        for nu, iu in enumerate(iunique):
            this_bool = index == iu
            collector[iu]  += (cl['density'][this_bool]*cl['d%s'%x1][this_bool]).sum().v
            #print collector[iu]
    plt.hist(np.log10(collector),histtype='step', color='b')

    plt.savefig('p29_hist2.pdf')
