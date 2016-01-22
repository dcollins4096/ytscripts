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

def log_fun(field):
    return np.log10(field)

if 'use_high_res' not in dir():
    use_high_res=True
error_list = []
if 0:
    setname = 'b02'
    frame = 20
    resolution = 256
    filename_prefix = 'b02_512_0020_smoothed_0256_density'
if 1:
    setname = 'b02'
    frame = 70
    resolution = 256
    filename_prefix = 'b02_512_0070_smoothed_0256_density'
skip_list={'b02':{20:[]}}
skip_list['b02'][20] = [0,5,30]

filament_set_res = 256.
alternate_set_res = 8192.
box_size_pc = 4.6
extraction_width_pc = 0.6
conversion_to_alternate = int(alternate_set_res/filament_set_res)
n_points_to_fit = 4

fig_list=[]
def newfig():
    fig_list.append(plt.figure())
    return fig_list[-1]
def point_selector(x,y,x0,x1,y0,y1,dx,dy):
    #F(x,y) = (y1-y0)*x + (x0-x1)*y +(x1*y0-x0*y1)
    #F==0, (x,y) is on the line.  <0 above, >0 below
    #All above or below, then there's no intersection. 
    #Sign change means that two points are on opposite sides of the line.
    c1 = y1-y0
    c2 = x0-x1
    c3 = x1*y0-x0*y1
    a   = c1*(x      )+c2*(y      )+c3
    amm = c1*(x-dx/2.)+c2*(y-dy/2.)+c3
    amp = c1*(x-dx/2.)+c2*(y+dy/2.)+c3
    app = c1*(x+dx/2.)+c2*(y+dy/2.)+c3
    apm = c1*(x+dx/2.)+c2*(y-dy/2.)+c3
    keep = amm*amp <= 0
    keep = np.logical_or( keep, amp*app <= 0 )
    keep = np.logical_or( keep, app*apm <= 0 )
    keep = np.logical_or( keep, apm*amm <= 0 ) 
    #pdb.set_trace()
    return keep

if 'dirname' not in dir():
    #dirname = '/scratch1/share/fils'
    dirname = '/Users/dcollins/RESEARCH2/Paper37_Philaments/2015-06-12-disperse/DATA'
    fname1 ='%s/%s.%s'%(dirname,filename_prefix,'fits_c0.1.up.NDskl.fits')
    fname2 ='%s/%s.%s'%(dirname,filename_prefix,'fits_c0.25.up.NDskl.fits')
    fname3 ='%s/%s.%s'%(dirname,filename_prefix,'fits_c0.5.up.NDskl.fits')
    fname4 ='%s/%s.%s'%(dirname,filename_prefix,'fits_c1.up.NDskl.fits')
    fname5 ='%s/%s.%s'%(dirname,filename_prefix,'fits.up.NDskl.fits')
    fname6 ='%s/%s.%s'%(dirname,filename_prefix,'fits_c0.25.up.NDskl.fits')
    fname7 ='%s/%s.%s'%(dirname,filename_prefix,'fits_c0.5.up.NDskl.fits')
    fname8 ='%s/%s.%s'%(dirname,filename_prefix,'fits_c1.up.NDskl.fits')
    fname9 ='%s/%s.%s'%(dirname,filename_prefix,'fits.up.NDskl.fits')

    set1_py = pyfits.open(fname4)
    set1_data_nobuf = set1_py[0].data
    set1_data = np.zeros([256,256])
    set1_data[2:-2,2:-2] = set1_data_nobuf


    fname_full = '%s/%s_512_%04d_smoothed_%04d_density.fits'%(dirname,setname,frame,resolution)
    """these sets lost two zones to smoothing.  Restore."""
    fullset_nobuf = pyfits.open(fname_full)[0].data
    fullset = np.zeros([256,256])
    fullset[2:-2,2:-2] = fullset_nobuf
    fname_full_high = '%s/%s_512_%04d_%04d_projection_density.fits'%(dirname,setname,frame,8192)
    fullset_high = pyfits.open(fname_full_high)[0].data

actual_resolution = fullset.shape[0]
if 1:
    #"""plot""" major image.  
    if 1:
        """all the filaments"""
        filament_image = newfig()
        filament_ax = filament_image.add_subplot(111)
        d1=filament_ax.imshow(set1_data)
        filament_image.colorbar(d1)
        outname = 'filament_image_%s_n%04d_r%04d.png'%(setname, frame, resolution)
        filament_image.savefig(outname)
        print 'filaments', outname

    if 0:
        """the major image"""
        high_res_image = newfig()
        high_res_ax = high_res_image.add_subplot(111)
        image=high_res_ax.imshow( np.log10(fullset),cmap='gray',interpolation='nearest')
        high_res_image.colorbar(image)
        outname = 'high_res_no_filaments_%s_n%04d_r%04d.png'%(setname, frame, resolution)
        high_res_image.savefig(outname)
        print "full image (no filaments)", outname

    """for plotting the profiles."""
    profile_image = newfig()
    profile_ax = profile_image.add_subplot(111)

    profile_high_image = newfig()
    profile_high_ax = profile_high_image.add_subplot(111)

    low_res_image = newfig()
    low_res_ax = low_res_image.add_subplot(111)
    image=low_res_ax.imshow( np.log10(fullset),cmap='gray',interpolation='nearest')

    if use_high_res:
        high_res_image = newfig()
        high_res_ax = high_res_image.add_subplot(111)
        high_res_ax.imshow( np.log10(fullset_high),cmap='gray',interpolation='nearest')
        outname = 'high_res_high_image_%s_n%04d_r%04d.png'%(setname, frame, resolution)
        high_res_image.savefig(outname)
        print "full image (high) ", outname

    dx = 1 #1./256
    left = actual_resolution
    Npoints = actual_resolution*1j
    #yup.  cell centered, needs to start and end at half steps, use the exact number.
    #x, y = np.ogrid[0.5*dx:1-0.5*dx:256j, 0.5*dx:1-0.5*dx:256j]
    #x, y = np.mgrid[0.5*dx:1-0.5*dx:256j, 0.5*dx:1-0.5*dx:256j]
    y, x = np.mgrid[0.5*dx:left-0.5*dx:Npoints, 0.5*dx:left-0.5*dx:Npoints]

    if use_high_res:
        left=alternate_set_res; Npoints=alternate_set_res*1j
        yL, xL = np.mgrid[0.5*dx:left-0.5*dx:Npoints, 0.5*dx:left-0.5*dx:Npoints]

    filament_cmap=rainbow_map(set1_data.max()+1)
    for nfil in np.unique(set1_data)[1:]:
        if nfil not in [48]: #, 55]: # range(48,60): #range(105,170):
            continue
#       if nfil not in [115]: #  [48]: # range(45,60): #[18]: #in skip_list[setname][frame]:
#           continue

        print "filament", nfil
        """for the profiles"""
        profile_ax.clear() # = profile_image.add_subplot(111)
        profile_high_ax.clear() # = profile_image.add_subplot(111)


        """Masks and positions"""
        this_fil = set1_data == nfil
        x_fil = x[ this_fil]
        y_fil = y[ this_fil]
        if nfil in [48]:
            keep = x_fil < 15
            x_fil = x_fil[keep]
            y_fil = y_fil[keep]
        x_centroid = x_fil.sum()/x_fil.size
        y_centroid = y_fil.sum()/y_fil.size

        """Check for periodic jumps."""
        if (x_fil < 0.51).any() or (y_fil < 0.51).any():
            print "Filament", nfil, "probably has periodic wrap."



        """ Overplot the filaments """
        if 0:
            if use_high_res:
                high_res_ax.scatter(x_fil.astype('int')*conversion_to_alternate, y_fil.astype('int')*conversion_to_alternate, 
                                    marker='o', linewidths=0,s=0.5, c=filament_cmap(nfil))
                high_res_ax.text(x_centroid*conversin_to_alternate,y_centroid*conversin_to_alternate, 
                                 "%d"%nfil, fontsize=5, color=filament_cmap(nfil))
            low_res_ax.scatter(x_fil.astype('int'),y_fil.astype('int'), marker='o',c=filament_cmap(nfil), linewidths=0, s=0.5)
            low_res_ax.text(x_centroid,y_centroid, "%d"%nfil, fontsize=5, color=filament_cmap(nfil))



        
        if 0:
            print "Skipping the profiles"
            continue
        """Fit the line"""
        if x_fil.size < 2:
            print "Filament", nfil, "only has one point!"
            continue

        fit = scatter_fit.scatter_fit(None,x_fil,y_fil, plot_points=False)#x_fil,y_fil)
        slope = fit['fit'][0]
        offset = fit['fit'][1]
        """For actual profiles along the filament"""
        n_points = x_fil.size
        rmap = rainbow_map(n_points+1)
        for ind in  range(n_points): 
            if nfil in [48] and x_fil[ind] > 14.8:
                continue
        
            dx=1.
            dy=1.

            #ind = x_fil.size/2 
            spine_point_x = x_fil[ind]; spine_point_y = y_fil[ind]
            width = extraction_width_pc/box_size_pc*filament_set_res
            """the transverse line"""
            r = (spine_point_x-x_fil)**2 + (spine_point_y - y_fil)**2
            r_args = np.argsort(r)
            x_closest = x_fil[r_args][0:n_points_to_fit]
            y_closest = y_fil[r_args][0:n_points_to_fit]
            #poly fit has a problem with vertical lines.
            vertical_problem = np.abs((x_closest-x_closest.mean())).max() < 0.5*dx
            horizontal_problem = np.abs((y_closest-y_closest.mean())).max() < 0.5*dy
            if vertical_problem:
                Deltax = width
                Deltay = 0
                #y0 = max([spine_point_y - 0.5*Deltay,0.5]); y1 = min([spine_point_y + 0.5*Deltay, y.max()])
                y0=y1=spine_point_y
                x0 = max([spine_point_x - 0.5*Deltax,0.5]); x1 = min([spine_point_x + 0.5*Deltax, x.max()])
            else:
                #fit = scatter_fit.scatter_fit(None,x_closest,y_closest, plot_points=False)#x_fil,y_fil)
                to_plot = None; #low_res_ax
                fit = scatter_fit.scatter_fit(to_plot,x_closest,y_closest, plot_points=False)#x_fil,y_fil)
                slope = fit['fit'][0]
                offset = fit['fit'][1]
                perp_slope = -1./slope
                Theta = np.arctan(perp_slope)
                if slope < 1e-5:
                    horizontal_problem = True
                #Spatial extent in line, targeting 0.6 pc for the filament
                Deltax = width * np.cos(Theta) 
                Deltay = width * np.sin(Theta) 
                x0 = max([spine_point_x - 0.5*Deltax,0.5]); x1 = min([spine_point_x + 0.5*Deltax, x.max()])
                y0 = max([spine_point_y - 0.5*Deltay,0.5]); y1 = min([spine_point_y + 0.5*Deltay, y.max()])
            yLeft = min([y0,y1]); yRight = max([y0,y1])
            xLeft = min([x0,x1]); xRight = max([x0,x1])
                

            #plt.plot([x0,x1],[y0,y1],c='g') #fit ranges

            slice_x = slice(np.where(x[0,:]<=xLeft)[0][-1], np.where(x[0,:]>=xRight)[0][0]+1)
            slice_y = slice(np.where(y[:,0]<=yLeft)[0][-1], np.where(y[:,0]>=yRight)[0][0]+1)
            x_sub=x[slice_y,slice_x]
            y_sub=y[slice_y,slice_x]
            keep = point_selector(x_sub,y_sub,x0,x1,y0,y1,dx,dy)
            x_a = x_sub[keep]
            y_a = y_sub[keep]
            sign_of_line = np.sign(x_a-spine_point_x)
            sign_of_line_y = np.sign(y_a-spine_point_y) #this funny sign juggle takes care of vertical points
            sign_of_line[sign_of_line==0] = sign_of_line_y[sign_of_line==0]

            coordinates = sign_of_line*np.sqrt((x_a-spine_point_x)**2+(y_a-spine_point_y)**2) #This centers the profile on the Filament.
            coordinates *= box_size_pc/filament_set_res
            sort_coord = np.argsort(coordinates)
            density = fullset[slice_y,slice_x][keep] # fullsetnar([fullset[ix,iy] for ix,iy in zip(x_i,y_i)])
            profile_ax.plot(coordinates[sort_coord], log_fun(density)[sort_coord], marker='o',c=rmap(ind))
            #profile_ax.scatter((spine_point_x-x_fil[0])*4.6/256,0.3,c=rmap(ind))
            profile_ax.text((spine_point_x-x_fil[0])*box_size_pc/filament_set_res,1.2,"%d"%ind,color=rmap(ind))
            low_res_ax.scatter(x_a,y_a, marker='o',c=rmap(ind), linewidths=0, s=0.1)
            low_res_ax.scatter(spine_point_x,spine_point_y, marker='*',c='k', linewidths=0, s=1)
            low_res_ax.scatter(x_closest,y_closest, marker='*',c='g', linewidths=0, s=0.5)

            if use_high_res:
                """now repeat for high res"""
                thexL = conversion_to_alternate*spine_point_x; theyL = conversion_to_alternate*spine_point_y
                width = extraction_width_pc/box_size_pc*alternate_set_res
                slice_x = slice(np.where(xL[0,:]<conversion_to_alternate*xLeft)[0][-1], 
                                np.where(xL[0,:]>conversion_to_alternate*xRight)[0][0]+1)
                slice_y = slice(np.where(yL[:,0]<conversion_to_alternate*yLeft)[0][-1]
                                , np.where(yL[:,0]>conversion_to_alternate*yRight)[0][0]+1)
                x_sub=xL[slice_y,slice_x]
                y_sub=yL[slice_y,slice_x]
                dx=1.
                dy=1.
                keep = point_selector(x_sub,y_sub,conversion_to_alternate*x0,conversion_to_alternate*x1,
                                      conversion_to_alternate*y0,conversion_to_alternate*y1,dx,dy)
                x_aL = x_sub[keep]
                y_aL = y_sub[keep]

                new_center_x = thexL; new_center_y = theyL
                density = fullset_high[slice_y,slice_x][keep]
                if 1:
                    """ Shift the curve to match the peak.  To avoid multiple maxima, restrict to a window of 0.1 pc """
                    peak_hunt_subset = (x_aL - thexL)**2+(y_aL-theyL)**2 < (0.05/4.6*8192)**2 #0.05 pc in pixels
                    max_coord=int( np.where(density[peak_hunt_subset]==density[peak_hunt_subset].max())[0].mean() )
                    new_center_x =  x_aL[peak_hunt_subset][max_coord]
                    new_center_y =  y_aL[peak_hunt_subset][max_coord]

                if (new_center_x-thexL)**2+(new_center_y-theyL)**2 > width**2:
                    error = "Error: Peak for filament %d index %d shift by more than 0.25 width"%(nfil, ind)
                    error_list.append(error)
                    #print error
                sign_of_line = np.sign(x_aL-new_center_x)
                sign_of_line_y = np.sign(y_aL-new_center_y)
                sign_of_line[sign_of_line==0] = sign_of_line_y[sign_of_line==0]

                coordinates = sign_of_line*np.sqrt((x_aL-new_center_x)**2+(y_aL-new_center_y)**2) #needs to be centered somehow.
                coordinates *= box_size_pc/alternate_set_res
                sort_coord = np.argsort(coordinates)

                profile_high_ax.plot(coordinates[sort_coord], log_fun(density)[sort_coord], marker='o',c=rmap(ind))
                #high_res_ax.scatter(new_center_x,new_center_y, marker='pngs',c=rmap(ind), linewidths=0,  s=0.01)
                high_res_ax.scatter(x_aL,y_aL, marker='o',c=rmap(ind), linewidths=0, s=0.01)
                high_res_ax.scatter(new_center_x,new_center_y, c='y', marker="*", linewidths=0,s=0.1)



        fname_profile='profile_center_%s_n%04d_r%04d_f%02d.pdf'%(setname, frame, resolution, nfil)
        low_res_name = 'low_res_image_%s_n%04d_r%04d_f%02d.png'%(setname, frame, resolution, nfil)
        #low_res_ax.plot([0,256],[0,256],c='b')
        low_res_ax.set_xlim(0,actual_resolution)
        low_res_ax.set_ylim(actual_resolution,0)
        low_res_image.savefig(low_res_name)
        print "save", low_res_name
        #fname_profile='profile_center_unshift_%02d.png'%nfil
        profile_ax.set_ylabel('log density')
        profile_ax.set_xlabel('r[pc]')
        profile_image.savefig(fname_profile)
        print "save", fname_profile
        if use_high_res:
            fname_profile='profile_high_%s_n%04d_r%04d_f%02d.pdf'%(setname, frame, resolution, nfil)
            profile_high_ax.set_ylabel('log density')
            profile_high_ax.set_xlabel('r[pc]')
            profile_high_image.savefig(fname_profile)
            print "save", fname_profile
            high_res_ax.set_xlim(0,alternate_set_res)
            high_res_ax.set_ylim(alternate_set_res,0)
            outname = 'high_res_high_image_%s_n%04d_r%04d_f%02d.pdf'%(setname, frame, resolution,nfil)
            high_res_image.savefig(outname)
            print "save", outname

this_format = 'pdf'
high_res_ax.set_xlim(0,8192)
high_res_ax.set_ylim(8192,0)
high_res_outname = 'high_res_all_filaments_%s_n%04d_r%04d.%s'%(setname,frame,resolution,this_format)
high_res_image.savefig(high_res_outname)
print "save", high_res_outname

low_res_name = 'low_res_all_filaments_%s_n%04d_r%04d.%s'%(setname, frame, resolution,this_format)
low_res_ax.set_xlim(0,actual_resolution)
low_res_ax.set_ylim(actual_resolution,0)
low_res_image.savefig(low_res_name)
print "save", low_res_name
for error in error_list:
    print error
for fig in fig_list:
    plt.close(fig)
    del fig
