
if hasattr(copy,'copy'):
    the_copy = copy.copy
else:
    the_copy = copy

if 1:
    dir1 = '/scratch1/dcollins/Paper05/OK4/'
    frame=700
    size = 128*2**4
    LOS = 'z'
    Blos = 'Bz'
    ds1 = yt.load(dir1+"/RS%04d/restart%04d"%(frame,frame))
    proj_flat = ds1.proj('density',LOS)
    proj_weight=ds1.proj(Blos, LOS, weight_field='cell_mass')

    frb_flat = proj_flat.to_frb(1,[size]*2)
    frb_weight = proj_weight.to_frb(1,[size]*2)

    density_flat = frb_flat['density']
    cell_mass_flat = frb_flat['cell_mass']
    Blos_weight = frb_weight[Blos]


    data_set = {'density':the_copy(density_flat.v.reshape(size,size,1)),
                Blos+"abs":np.abs(the_copy(Blos_weight.v.reshape(size,size,1))),
                'cell_mass':the_copy(cell_mass_flat.v.reshape(size,size,1))}
    bbox = np.array([[0.,1.]]*3)
    ds2 = yt.load_uniform_grid(data_set, [size,size,1], length_unit="Mpc", bbox=bbox)
    phase1=yt.PhasePlot(ds2.all_data(),'density',Blos+"abs",'cell_mass',weight_field=None)
    outname = 'p40_badone_abs.pdf'
    phase1.save(outname)
    print outname

if 1:
    ef('zeeman_measurements.py')
    #phase1.set_cmap('cell_mass','gray')
    ploot = phase1.plots[('gas','cell_mass')].axes
    ploot.scatter(ColumnDensity,np.abs(Bfield),c='r')
    ploot.scatter(TC[1],TC[0],marker='^', c='k',label='tC2008')
    ploot.scatter(FT[1],FT[0],marker='D',c='c', label = 'FT2008')
    outname = 'p40_2d3d.pdf'
    phase1.save(outname)
    print outname

