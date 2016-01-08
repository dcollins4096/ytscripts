if 'ef' not in dir():
    execfile('go')

"""
Ok, for Friday (tomorrow is the proposal.)
1.) Do a bunch.  Hand cull the reasonable parts.
2.) Histogram for the widths.  Variance for now.
3.) Do it for the late time version.
"""
import pyfits
import scatter_fit
reload(scatter_fit)

use_high_res=False
setname = 'b02'
frame = 20
resolution = 256
skip_list={'b02':{20:[]}}
skip_list['b02'][20] = [0,5,30]
fig_list=[]
def newfig():
    fig_list.append(plt.figure())
    return fig_list[-1]
if 'dirname' not in dir():
    #dirname = '/scratch1/share/fils'
    dirname = '/Users/dcollins/RESEARCH2/Paper37_Philaments/2015-06-12-disperse/DATA'
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

    fname_full = '%s/%s_512_%04d_smoothed_%04d_density.fits'%(dirname,setname,frame,resolution)
    fullset = pyfits.open(fname_full)[0].data
    fname_full_high = '%s/%s_512_%04d_%04d_projection_density.fits'%(dirname,setname,frame,8192)
    fullset_high = pyfits.open(fname_full_high)[0].data

actual_resolution = fullset.shape[0]
if 1:
    #"""plot""" major image.  
    #filament_image = newfig()
    #filament_ax = filament_image.add_subplot(111)
    #d1=filament_ax.imshow(set1_data)
    #filament_image.colorbar(d1)
    #outname = 'filament_image_%s_n%04d_r%04d.png'%(setname, frame, resolution)
    #filament_image.savefig(outname)
    #print 'filaments', outname

    full_image = newfig()
    full_ax = full_image.add_subplot(111)
    image=full_ax.imshow( np.log10(fullset),cmap='gray',interpolation='nearest')
    full_image.colorbar(image)
    outname = 'full_image_%s_n%04d_r%04d.png'%(setname, frame, resolution)
    full_image.savefig(outname)
    print "full image", outname

    if use_high_res:
        full_high_image = newfig()
        full_high_ax = full_high_image.add_subplot(111)
        full_high_ax.imshow( np.log10(fullset_high),cmap='gray',interpolation='nearest')
        outname = 'full_high_image_%s_n%04d_r%04d.png'%(setname, frame, resolution)
        full_high_image.savefig(outname)
        print "full_high_image", outname
        full_image.savefig(outname)

    print outname
    dx = 1 #1./256
    left = actual_resolution
    Npoints = actual_resolution*1j
    #yup.  cell centered, needs to start and end at half steps, use the exact number.
    #x, y = np.ogrid[0.5*dx:1-0.5*dx:256j, 0.5*dx:1-0.5*dx:256j]
    #x, y = np.mgrid[0.5*dx:1-0.5*dx:256j, 0.5*dx:1-0.5*dx:256j]
    y, x = np.mgrid[0.5*dx:left-0.5*dx:Npoints, 0.5*dx:left-0.5*dx:Npoints]
    distance = np.sqrt(x ** 2 + y ** 2)

    if use_high_res:
        left=8192; Npoints=8192j
        yL, xL = np.mgrid[0.5*dx:left-0.5*dx:Npoints, 0.5*dx:left-0.5*dx:Npoints]
    for nfil in np.unique(set1_data):
        if nfil not in [18]: #in skip_list[setname][frame]:
            continue

        """for the profiles"""
        profile_image = newfig()
        profile_ax = profile_image.add_subplot(111)

        """for filament overlays"""
        fil_image = newfig()
        fil_ax = fil_image.add_subplot(111)
        fil_im = fil_ax.imshow( np.log10(fullset),cmap='gray',interpolation='nearest')
        fil_image.colorbar(fil_im)

        """Masks and positions"""
        this_fil = set1_data == nfil
        not_fil = set1_data != nfil
        x_fil = x[ this_fil]
        y_fil = y[ this_fil]
        fil_ax.scatter(x_fil.astype('int'),y_fil.astype('int'), marker='o',c='g', linewidths=0, s=1)

        """ lets see"""
        if use_high_res:
            full_high_ax.scatter(x_fil*32, y_fil*32)
            #full_high_image.savefig('test8192.png')
            fil_ax.scatter(x_fil,y_fil)

        """For actual profiles along the filament"""
        n_points = x_fil.size
        rmap = rainbow_map(n_points+1)

        """Fit the line"""
        fit = scatter_fit.scatter_fit(None,x_fil,y_fil, plot_points=False)#x_fil,y_fil)
        slope = fit['fit'][0]
        offset = fit['fit'][1]
        
        for ind in range(n_points): #point along the filament, at which the transverse measurement is done.

            #ind = x_fil.size/2 
            spine_point_x = x_fil[ind]; spine_point_y = y_fil[ind]
            """the transverse line"""
            width = 0.6/4.6*256
            perp_slope = -1./slope
            Theta = np.arctan(perp_slope)
            #Spatial extent in line, targeting 0.6 pc for the filament
            Deltax = width * np.cos(Theta) #there is a better way for this.
            Deltay = width * np.sin(Theta) #there is a better way for this.
            x0 = spine_point_x - 0.5*Deltax; x1 = spine_point_x + 0.5*Deltax
            y0 = spine_point_y - 0.5*Deltay; y1 = spine_point_y + 0.5*Deltay
            #plt.plot([x0,x1],[y0,y1],c='g') #fit ranges

            """now I need to get a range of dx"""
            #select all the x indices within the extents of "width",  then 
            #pick y indices for the x points on the line.
            x_i_logic = np.logical_and( x[0,:] >= x0, x[0,:] <= x1 )
            x_a = x[0,:][x_i_logic] #these are cell centered
            x_i = x_a.astype('int') #these are now indices
            y_a = perp_slope*(x_a - spine_point_x) + spine_point_y
            left_edge = 0; dx=1 #using code units.
            y_i = ((y_a-left_edge)/dx).astype('int')
            on_the_plot = np.logical_and(x_i < actual_resolution, y_i < actual_resolution)
            x_i = x_i[on_the_plot]; y_i = y_i[on_the_plot];
            x_a = x_a[on_the_plot]; y_a = y_a[on_the_plot];
            #coordinates = np.sqrt(x_a**2+y_a**2) #needs to be centered somehow.
            coordinates = np.sign(x_a-spine_point_x)*np.sqrt((x_a-spine_point_x)**2+(y_a-spine_point_y)**2) #This centers the profile on the Filament.
            coordinates *= 4.6/256
            density = nar([fullset[ix,iy] for ix,iy in zip(x_i,y_i)])
            profile_ax.plot(coordinates, density, marker='o',c=rmap(ind))
            profile_ax.scatter((spine_point_x-x_fil[0])*4.6/256,0.3,c=rmap(ind))
            fil_ax.scatter(x_i,y_i, marker='o',c=rmap(ind), linewidths=0, s=0.1)
            #fil_ax.scatter(spine_point_x,they,c=rmap(ind), s=1, linewidth=0)

            if use_high_res:
                """now repeat for high res"""
                thexL = 32*thex; theyL = 32*they
                x_i_logicL = np.logical_and( xL[0,:] >= 32*x0, xL[0,:] <= 32*x1 )
                x_aL = xL[0,:][x_i_logicL] #these are cell centered
                x_iL = x_aL.astype('int') #these are now indices
                y_aL = perp_slope*(x_aL - thexL) + theyL
                left_edge = 0; dx=1
                y_iL = ((y_aL-left_edge)/dx).astype('int')
                on_the_plot = np.logical_and(x_iL < 8192, y_iL < 8192)
                x_iL = x_iL[on_the_plot]; y_iL = y_iL[on_the_plot];
                x_aL = x_aL[on_the_plot]; y_aL = y_aL[on_the_plot];
                full_high_ax.scatter(x_iL,y_iL, marker='o',c='y',s=1)
          
                #coordinates = np.sqrt(x_a**2+y_a**2) #needs to be centered somehow.
                coordinates = np.sign(x_aL-thexL)*np.sqrt((x_aL-thexL)**2+(y_aL-theyL)**2) #needs to be centered somehow.
                coordinates *= 4.6/8192
                density = nar([fullset_high[ixL,iyL] for ixL,iyL in zip(x_iL,y_iL)])
                profile_ax.plot(coordinates, density, marker='o',c=rmap(ind))
                #profile_ax.scatter((thex-x_fil[0])*4.6/256,0.3,c=rmap(ind))



        fname_profile='profile_center_%s_n%04d_r%04d_f%02d.pdf'%(setname, frame, resolution, nfil)
        fil_name = 'fil_image_%s_n%04d_r%04d_f%02d.pdf'%(setname, frame, resolution, nfil)
        #fil_ax.plot([0,256],[0,256],c='b')
        fil_ax.set_xlim(0,actual_resolution)
        fil_ax.set_ylim(actual_resolution,0)
        fil_image.savefig(fil_name)
        print "skipping", fil_name
        #fname_profile='profile_center_unshift_%02d.png'%nfil
        profile_ax.set_ylabel('density')
        profile_ax.set_xlabel('r[pc]')
        profile_image.savefig(fname_profile)
        print "profile", fname_profile
        full_outname = 'full_fil_%s_n%04d_r%04d.png'%(setname,frame,resolution)
        full_image.savefig(full_outname)
        print full_outname
        outname = 'full_high_image_%s_n%04d_r%04d_f%02d.png'%(setname, frame, resolution,nfil)
        if use_high_res:
            full_high_image.savefig(outname)
        print "full_high_image, fil", outname

for fig in fig_list:
    plt.close(fig)
    del fig
