#line of sight is z
if 'ef' not in dir():
    execfile('go')
ef('p37_misc.py')

highfil1=[(40,104),(55,162)]
highfil2=[(146,179),(21,108)]

lowfil2=[(119,254),(333,401)]
lowfil1=[(148,311), (196,308)]
#n=0,0.3,0.7
plt.clf()

if 1:
    #frame=40
    #frame=70
    ds = yt.load('/scratch1/dcollins/Paper08/B02/512/RS%04d/restart%04d'%(frame,frame)); simname = 'high_field'
    subset = 'full'
    if 0:
        L =nar([119./512,333./512, 0.])
        R =nar([254./512,401./512,1.])
    if 0:
        L = nar([40./512,55./512,0])
        R = nar([104./512,162./512,1])
    if 0:
        #x=0:20, y
        #high 1.5, zoom in
        L = nar([(70.)/512,55./512,0])
        R = nar([(100.)/512   ,(55.+25.)/512,1])
    if 0:
        subset = 'highfil_1_5'
        #x=0:20, y
        #high 1.5, zoom in
        L = nar([(70.)/512,55./512,0])
        R = nar([(100.)/512   ,(55.+30.)/512,1])
    if 1:
        #high, early, three filaments
        subset = '3fil'
        C = nar([208.,135.])/256 #Center, pixels
        S = nar([50.,50.])/256   #width
        L, R = np.zeros(3), np.zeros(3)
        L[0] = C[0]-S[0]/2; L[1]=C[1]-S[1]/2; L[2]=0
        R[0] = C[0]+S[0]/2; R[1]=C[1]+S[1]/2; R[2]=1
    if 1:
        #high, early, three filaments
        subset = '3fil_Knot'
        C = nar([208.,135.])/256 #Center, pixels
        S = nar([50.,50.])/256   #width
        L, R = np.zeros(3), np.zeros(3)
        L[0] = C[0]-S[0]/2; L[1]=C[1]-S[1]/2; L[2]=0
        R[0] = C[0]+S[0]/2; R[1]=C[1]+S[1]/2; R[2]=1
        L[0] += 25./256
    if 1:
        subset='phil_trace_1'
        L, R = np.zeros(3), np.zeros(3)
        off = [-0.025, 0.14]
        #This is a complex way to do things, due to way things happened with Phil.
        anchor={80:[0.5,1.5],70:[0.5,1.5],60:[0.5,1.0],50:[0.5,1.0],40:[0.5,1.0],30:[0.5,1.0]}
        trace={80:[-571.,-492.],70:[-632.,-584.],60:[-616.,0],50:[-616.,-149.],40:[-616-60,-140],
               30:[-670,-145]}
        L[0] = anchor[frame][0]/4.6+trace[frame][0]/8192. + off[0]
        L[1] = anchor[frame][1]/4.6 +trace[frame][1]/8192. + off[1]
        R[0] = L[0]+700./8192
        R[1] = L[1]+700./8192
        L[2]=0.0; R[2]=1.0
        C = 0.5*(L+R)
        z_extents={80:[0.175,0.255]}
   
    if 0:
        L[2]=z_extents[frame][0]
        R[2]=z_extents[frame][1]
    L=bu(L,1./512)
    R=bu(R,1./512)
    C = 0.5*(L+R)
    reg = ds.region(C,L,R)
    if 0:
        subset='FULL'
        field = 'density'
        proj = ds.proj(field,"z")
        pw=proj.to_pw()
        axis = "z"
        ProjL = nar([0.,0.])
        ProjR = nar([1.,1.])
        if 0:
            pw=proj.to_pw()
            pw.set_cmap('density','gray')
            #pw.annotate_line(L,[R[0],R[1],1], coord_system = 'data', plot_args={'color':'r'})
            pw.annotate_line(L,[R[0],L[1],1], coord_system = 'data', plot_args={'color':'r'})
            pw.annotate_line(L,[L[0],R[1],1], coord_system = 'data', plot_args={'color':'r'})
            pw.annotate_line([L[0],R[1],0],R, coord_system = 'data', plot_args={'color':'r'})
            pw.annotate_line([R[0],L[1],0],R, coord_system = 'data', plot_args={'color':'r'})
            radius = 0.1/4.6
            #pw.annotate_sphere(L-1.5*radius, radius, coord_system='data', circle_args={'color':'r'})
            #pw.annotate_line(L,R, coord_system = 'data', plot_args={'color':'r'})
            #pw.annotate_line(L,R, coord_system = 'data', plot_args={'color':'r'})
            #pw.annotate_line(L,R, coord_system = 'data', plot_args={'color':'r'})
            print pw.save('philaments_high_1.5_t%04d'%frame)
    if 1:
        axis = "z"
        ProjL = nar([L[0], L[1]])
        ProjR = nar([R[0], R[1]])
        proj = ds.proj("density","z",data_source=reg,center=C); field = 'density'
        #proj = ds.proj("grid_level","z",data_source=reg,center=C,method='mip'); field = 'grid_level'
        if 0:
            pw=proj.to_pw(center=C)
            pw.set_origin('domain')
            D = R-L
            pw.set_width(max([D[0],D[1]]))
            pw.set_cmap('density','gray')
            try:
                pw.set_log('grid_level',False)
                pw.set_zlim('grid_level',0,4)
            except:
                pass
            print pw.save('p37_trace1_full_t%04d'%frame)
            pw.set_width(max(ProjR-ProjL))
            pw.annotate_streamlines('Bx','By',factor=32,plot_args={'color':'b'})
            pw.annotate_streamlines('x-velocity','y-velocity',factor=32,plot_args={'color':'r'})
            print pw.save('p37_trace1_vel_t%04d'%frame)
    if 0:
        axis = "y"
        #ds.units_override={'length_unit':(4.6,'pc')}
        proj = ds.proj("density","y",data_source=reg,center=C)
        if 1:
            pw=proj.to_pw(center=C,origin='domain')
            pw.set_cmap('density','gray')
            pw.set_width(max(R-L))
            ProjL = nar([L[2], L[0]])
            ProjR = nar([R[2], R[1]])
            #pw.set_unit('length_unit',(4.6,'pc'))
            #pw.annotate_magnetic_field()
            #pw.annotate_streamlines('Bz','Bx',factor=32,plot_args={'color':'b'})
            #pw.annotate_streamlines('z-velocity','x-velocity',factor=32,plot_args={'color':'r'})
            #pw.annotate_velocity() #plot_args={'color':'r'})
            print pw.save('p37_trace1_vel_t%04d'%frame)
    #proj = projo.to_frb(1,[512,512])

#failsafe as I moved to transverse projections
#del R
#del L
print "on to the next"
if 1:
    mask = proj['density']>0
    nmask = proj['density']<=0
    total_zones = 8192
    Nzones = max(ProjR[0]-ProjL[0], ProjR[1]-ProjL[1])*total_zones
    resolution_name = "%d"%Nzones
    width = max(ProjR-ProjL)
    res = 1./Nzones
    frb = proj.to_frb( width, int(Nzones))
    #mag = copy.copy(frb['magnetic_pressure'].v)
    if 0:
        mask2  = den >0
        mask2o = den <=0
        den[mask2o] = den[mask2].min()
    #make a subset
    if 0:
        den_damned_square = copy.copy(frb[field].v)
        shape_x, shape_y =  den_damned_square.shape
        half_x = shape_x/2; half_y = shape_y/2
        ok_y = np.where( den_damned_square[half_x,:] >0 )[0]
        ok_x = np.where( den_damned_square[:,half_y] >0 )[0]
        den = den_damned_square[ok_x[0]:ok_x[-1]+1, ok_y[0]:ok_y[-1]+1]
    if 0:
        """simple plot."""
        plt.clf()
        plt.imshow( np.log10( den), origin='lower', interpolation='nearest',cmap='gray')
        outname = 'filament_image_%s_%s_t%04d_%s_%s.png'%(simname, subset, frame,resolution_name, field)
        plt.savefig(outname)
        print outname
    if 1:
        den = copy.copy(frb[field])
    density_factor = 1000*4.6*3.08e18 #1000 cm^-3 * 4.6 pc * cm/pc
    if 1:
        if field == 'density':
            to_plot = np.log10( density_factor* den)
            to_write = copy.copy(to_plot)
            cb_label= r'$\log_{10}( N[\rm{cm}^{-2}]=N[\rm{code}]*%0.2e )$'%density_factor
            norm = None #matplotlib.colors.Normalize(21,24)
        elif field == 'grid_level':
            to_plot = den
            cb_label = r'$level$'
            norm = None

        else:
            raise
        if 1:
            distance = 140. #pc
            resolution = 80.#arcsec
            resolution_name = '0.05pc'
            smooth_pc = distance*resolution/206264. #in pc
            smooth_pc = 0.05
            smooth_px = max(int(smooth_pc/(4.6/total_zones)), 2)
            #smooth_px = 10
            ef('p37_blur.py')
            smoothed_buffer = np.log10(blur_image(density_factor* den.v,smooth_px))
            to_plot = np.zeros_like(den.v)
            to_plot[smooth_px:-smooth_px, smooth_px:-smooth_px] = smoothed_buffer
            to_write = copy.copy(to_plot)
            mask = to_plot == 0
            nonmask = to_plot > 0
            non_zero_min  = to_plot[nonmask].min()
            to_plot[mask] = non_zero_min

            

            norm=None #matplotlib.colors.Normalize(21,24)

        """code-heinous, readout-good plot"""
        plt.clf()
        plt.imshow(to_plot, origin='lower', interpolation='nearest',cmap='gray',  norm=norm)
        cb=plt.colorbar()
        cb.set_label(cb_label)
        dx = 4.6/8193.
        stupid_slice = slice(None)
        old_xticks = plt.xticks()[0][stupid_slice]
        plt.xticks( old_xticks, ["$%0.2f$"%n for n in  old_xticks*dx+ProjL[0] ] )
        plt.xlabel(r'$x[\rm{pc}]\  (\Delta x=$%s$\rm{pc})$'%expform(dx))
        old_yticks = plt.yticks()[0][stupid_slice]
        plt.yticks( old_yticks, ["$%0.2f$"%n for n in  old_yticks*dx+ProjL[1] ] )
        plt.ylabel(r'$y[\rm{pc}]\  (\Delta x=$%s$\rm{pc})$'%expform(dx))
        outname = 'p37_%s_%s_t%04d_%s_%s_%s.png'%(simname, subset, frame,resolution_name, field,axis)
        plt.xlim(0,to_plot.shape[0])
        plt.ylim(0,to_plot.shape[1])
        plt.savefig(outname)
        print outname

if 0:
    plt.clf()
    fig = plt.figure(figsize=[81.92]*2, dpi=100)
    fig.figimage(np.log10(den))
    outname = 'full_res_%s_%s_t%04d_%s'%(simname, subset, frame, resolution_name)
    fig.savefig(outname)
    print outname

if 0:
    """ histogram """
    plt.clf()
    outname = 'histogram_%s_%s_t%04d_%s'%(simname, subset, frame, resolution_name)
    mn=-0.75; mx=1.75; nbins=100; dx=(mx-mn)/nbins; bins=na.arange(mn,mx+dx,dx)
    plt.hist(np.log10(den.flatten()),bins=bins,histtype='step',normed=True)
    plt.xlabel(r'$\log \Sigma$')
    plt.ylabel(r'$N$')
    plt.savefig(outname)
    print outname

# hdu = pyfits.PrimaryHDU(den)
# hdulist = pyfits.HDUList([hdu])
# hdulist.writeto('MOVEME.fits')
# hdulist.close()

if 1:
    """ write text"""
    field = 'density'
    fname = 'filament_data_%s_%s_t%04d_%s_%s.txt'%(simname, subset, frame,resolution_name, field)
    outfile = open(fname,'w')
#    den = copy.copy(frb['density'].v)
    thisone = to_write # copy.copy(frb[field].v)
    for n, full_row in enumerate(thisone):
#       row = full_row[ full_row > 0]
#       N = len(row)
#       if N == 0:
#           continue
        outfile.write("%0.16e "*N%tuple(row))
        outfile.write("\n")
    outfile.close()
    print fname

if 0:
    """ write text transposed"""
    field = 'density'
    fname = 'filament_data_transposed_%s_%s_t%04d_%s_%s.txt'%(simname, subset, frame,resolution_name, field)
    outfile = open(fname,'w')
    thisone2 = copy.copy(frb[field].v).transpose()
    plt.clf()
    plt.imshow( np.log10( thisone2), origin='lower', interpolation='nearest',cmap='gray')
    oname2 = 'filament_image_transposed_%s_%s_t%04d_%s_%s.png'%(simname, subset, frame,resolution_name, field)
    plt.savefig(oname2)
    print oname2
    for n, full_row in enumerate(thisone2):
        row = full_row[ full_row > 0]
        N = len(row)
        if N == 0:
            continue
        outfile.write("%0.16e "*N%tuple(row))
        outfile.write("\n")
    outfile.close()
    print fname
   
if 0:
    fname = 'philaments_high_1.5_t0070_full_res_density.txt'
    thing = np.loadtxt(fname)
    plt.clf()
    for n in range(100,420,20):
        plt.clf()
        plt.plot(np.log10(thing[:,n]))
        plt.ylim((-0.2,2.))
        plt.savefig('philaments_high_1.5_yslice_%04d.png'%n)
     
     
if 0:
    den = copy.copy(frb['density'].v)
    thisone = copy.copy(frb[field].v)
    for n in range(100,420,20):
        plt.clf()
        fname = 'tmp_%04d.png'%n
        plt.plot(np.log10(den[:,n]))
        plt.ylim(-0.2,2)
        plt.savefig(fname)
     
   
if 0:
    fname  = 'philaments_high1.5_t%04d.txt'%frame
    outfile = open(fname,'w')
    mask = proj['Density']>0
    d = proj['Density'][mask]
    x = proj['x'][mask]
    y = proj['y'][mask]
    dx = proj['dx'][mask]
    dy = proj['dy'][mask]
    outfile.write("%16s %16s %16s %16s %16s\n"%("n","x","y","dx","dy"))
    for n in range(d.size):
        outfile.write( "%16.4e %16.4e %16.4e %16.4e %16.4e\n"%(d[n],x[n],y[n],dx[n],dy[n]))
    outfile.close()
   

