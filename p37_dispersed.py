
import pyfits
import scatter_fit
reload(scatter_fit)


if 'dirname' not in dir():
    dirname = '/scratch1/share/fils'
    fname1 ='%s/%s'%(dirname,'b02_512_0020_smoothed_0256_density.fits_c0.1.up.NDskl.fits')
    fname2 ='%s/%s'%(dirname,'b02_512_0020_smoothed_0256_density.fits_c0.25.up.NDskl.fits')
    fname3 ='%s/%s'%(dirname,'b02_512_0020_smoothed_0256_density.fits_c0.5.up.NDskl.fits')
    fname4 ='%s/%s'%(dirname,'b02_512_0020_smoothed_0256_density.fits_c1.up.NDskl.fits')
    fname5 ='%s/%s'%(dirname,'b02_512_0020_smoothed_0256_density.fits.up.NDskl.fits')
    fname6 ='%s/%s'%(dirname,'b02_512_0070_smoothed_0256_density.fits_c0.25.up.NDskl.fits')
    fname7 ='%s/%s'%(dirname,'b02_512_0070_smoothed_0256_density.fits_c0.5.up.NDskl.fits')
    fname8 ='%s/%s'%(dirname,'b02_512_0070_smoothed_0256_density.fits_c1.up.NDskl.fits')
    fname9 ='%s/%s'%(dirname,'b02_512_0070_smoothed_0256_density.fits.up.NDskl.fits')

    set1_py = pyfits.open(fname4)
    set1_data = set1_py[0].data

    fname_full = '/scratch1/dcollins/Paper37_Filaments/b02_512_0020_smoothed_0256_density.fits'
    fullset = pyfits.open(fname_full)[0].data
    plt.imshow( np.log10(fullset))
    plt.savefig('sa.png')
if 1:
    """plot"""
    plt.clf()
    plt.imshow(set1_data)
    outname = 's4.png'
    plt.savefig(outname)
    print outname
    plt.clf()
    fname_full = '/scratch1/dcollins/Paper37_Filaments/b02_512_0020_smoothed_0256_density.fits'
    fullset = pyfits.open(fname_full)[0].data
    plt.imshow( np.log10(fullset),cmap='gray')
    dx = 1 #1./256
    left = 256.
    #yup.  cell centered, needs to start and end at half steps, use the exact number.
    #x, y = np.ogrid[0.5*dx:1-0.5*dx:256j, 0.5*dx:1-0.5*dx:256j]
    #x, y = np.mgrid[0.5*dx:1-0.5*dx:256j, 0.5*dx:1-0.5*dx:256j]
    y, x = np.mgrid[0.5*dx:left-0.5*dx:256j, 0.5*dx:left-0.5*dx:256j]
    distance = np.sqrt(x ** 2 + y ** 2)
    for nfil in np.unique(set1_data):
        if nfil not in [ 15]:
            continue
        this_fil = set1_data == nfil
        not_fil = set1_data != nfil
        x_fil = x[ this_fil]
        y_fil = y[ this_fil]
        plt.plot(x_fil,y_fil,c='r')
        cen_x = x_fil.sum()/x_fil.size
        cen_y = y_fil.sum()/y_fil.size
        if 1:
            """for imaging just this set"""
            #clone = copy.copy(set1_data)
            #clone[not_fil] = 0
            #print nfil, clone.sum()
            #plt.clf()
            #plt.imshow(clone) #, origin='lower')
            outname = 'this_fil_%02d.png'%nfil
            #plt.plot([93.5, 129.5], [20.701350979975842, 62.481034439683796],c='k')

            fit = scatter_fit.scatter_fit(None,x_fil,y_fil, plot_points=False)#x_fil,y_fil)

            print "fu"
            slope = fit['fit'][0]
            offset = fit['fit'][1]
            """pick a point, others possible"""
            ind = x_fil.size/2 
            thex = x_fil[ind]; they = y_fil[ind]
            """the transverse line"""
            width = 0.6/4.6*256
            perp_slope = -1./slope
            Theta = np.arctan(perp_slope)
            Deltax = width * np.cos(Theta) #there is a better way for this.
            Deltay = width * np.sin(Theta) #there is a better way for this.
            x0 = thex - 0.5*Deltax; x1 = thex + 0.5*Deltax
            y0 = they - 0.5*Deltay; y1 = they + 0.5*Deltay
            #plt.plot([x0,x1],[y0,y1],c='g') #fit ranges
            """now I need to get a range of dx"""
            x_i_logic = np.logical_and( x[0,:] >= x0, x[0,:] <= x1 )
            x_a = x[0,:][x_i_logic] #these are cell centered
            y_a = perp_slope*(x_i - thex) + they
            x_i = x_a.astype('int') #these are now indices
            left_edge = 0; dx=1
            y_i = ((y_a-left_edge)/dx).astype('int')
            plt.scatter(x_i,y_i, marker='o',c='y')
            coordinates = np.sqrt(x_a**2+y_a**2) #needs to be centered somehow.
            coordinates *= 4.6/256
            density = nar([fullset[ix,iy] for ix,iy in zip(x_i,y_i)])
        
            del fit

        if 1:
            plt.savefig(outname)
            print "fu", outname
        if 1:
            plt.clf()
            plt.plot(coordinates, density, marker='o',c='r')
            fname_profile='profile_center_%02d.png'%nfil
            plt.savefig(fname_profile)
            print 'fub', fname_profile

            



        
    

    
