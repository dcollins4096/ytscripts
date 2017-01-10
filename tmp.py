if 'ef' not in dir():
    execfile('go')

if 'aj19' not in dir():
    aj19=taxi.taxi('aj19_sphere')
    aj20=taxi.taxi('aj20_sphere')
    flt=taxi.fleet([aj19,aj20])
    #flt['frames']=[0,3]
    #flt.plot()

if 0:
    vmin=-5e-2
    vmax=5e-2
    norm=mpl.colors.Normalize(vmin=minmin,vmax=maxmax)
    plt.clf()
    plot=plt.imshow(two_diff, interpolation='nearest',origin='lower')
    colorbar=plt.colorbar(plot, norm=norm)
    colorbar.set_clim(minmin,maxmax)
    plt.title('%s cube'%field)
    outname = 'aj18_diff_%s_h5.png'%(field)
    plt.savefig(outname)
    print outname

if 1:
    import dsdiff
    reload(dsdiff)
    ud = dsdiff.udiff(aj19.directory,aj20.directory,frames=[1],grids=[4], fields=['Density'])
    ud(output='z',output_range=[3], zlim=[0,10], grid_direct=False, unit_function=lambda x: x.in_units('code_density').v)
