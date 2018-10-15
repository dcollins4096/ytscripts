def cuber(fname):
    fname_list = ['/scratch1/dcollins/Paper33_galaxies/aj18_refinement_test/%s.grid0001'%fname]
    fname_list += ['/scratch1/dcollins/Paper33_galaxies/aj18_refinement_test/%s.grid0002'%fname]

    out={}
    for n, f in enumerate(fname_list):
        print fname
        fptr = h5py.File(f,'r')
        print fptr.keys()
        try:
            for key in fptr.keys():
                out[key] = fptr[key][:][3:-3,3:-3,3:-3]
        except:
            raise
        finally:
            fptr.close()


    return out

try:
    del oot
    del diff
    del two_diff
except:
    pass
#norm=mpl.colors.Normalize(vmin=-0.5, vmax=1.5)

oot=cuber('data1110000')
for field in ['density','Metal_Density']:
    rho_0 =g[field].in_units('code_density').v 
    rho_1 = oot[field]
    diff = (rho_0-rho_1)/rho_0

    plt.clf()
    two_diff = np.sum(diff,axis=0)
    stat(diff, "3d")
    stat(two_diff,"2d")
    plot=plt.imshow(two_diff, interpolation='nearest',origin='lower')
    plt.colorbar(plot)
    plt.title('%s difference'%field)
    outname = 'aj18_diff_%s.png'%(field)
    plt.savefig(outname)
    print outname

    for_cbar = sorted(rho_0.flatten())
    min_index = int(0.25*len(for_cbar))
    max_index = int( 0.75*len(for_cbar))
    minmin = {'density':-6, 'Metal_Density':-21}[field] #for_cbar[min_index] #-6. #min([rho_1.min(), rho_0.min()])
    maxmax = {'density':0, 'Metal_Density':-16}[field] #for_cbar[max_index] #0.0 #min([rho_1.max(), rho_0.max()]) #yes, i mean min
    print "CBAR: ",field, minmin, maxmax
    norm=mpl.colors.Normalize(vmin=minmin,vmax=maxmax)
    plt.clf()
    two_diff = np.log(np.sum(rho_1,axis=0))
    plot=plt.imshow(two_diff, interpolation='nearest',origin='lower')
    colorbar=plt.colorbar(plot, norm=norm)
    colorbar.set_clim(minmin,maxmax)
    plt.title('%s cube'%field)
    outname = 'aj18_diff_%s_h5.png'%(field)
    plt.savefig(outname)
    print outname

    plt.clf()
    two_diff = np.log(np.sum(rho_0,axis=0))
    plot=plt.imshow(two_diff, interpolation='nearest',origin='lower')
    plt.title('%s enzo'%field)
    colorbar=plt.colorbar(plot, norm=norm)
    colorbar.set_clim(minmin,maxmax)
    outname = 'aj18_diff_%s_ds.png'%(field)
    plt.savefig(outname)
    print outname
