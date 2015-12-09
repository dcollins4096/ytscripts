#line of sight is z
if 'ef' not in dir():
    execfile('go')
ef('p20_misc.py')

highfil1=[(40,104),(55,162)]
highfil2=[(146,179),(21,108)]

lowfil2=[(119,254),(333,401)]
lowfil1=[(148,311), (196,308)]
#n=0,0.3,0.7
plt.clf()

if 1:
    frame=50
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
        anchor={80:[0.5,1.5],70:[0.5,1.5],60:[0.5,1.0],50:[0.5,1.0],40:[0.5,1.0],30:[0.5,1.0]}
        trace={80:[-571.,-492.],70:[-632.,-584.],60:[-616.,0],50:[-616.,-149.],40:[-616-60,-140],
               30:[-670,-145]}
        L[0] = anchor[frame][0]/4.6+trace[frame][0]/8192. + off[0]
        L[1] = anchor[frame][1]/4.6 +trace[frame][1]/8192. + off[1]
        R[0] = L[0]+700./8192
        R[1] = L[1]+700./8192
        L[2]=0.0; R[2]=1.0
        C = 0.5*(L+R)
        print "L",L
        print "R",R
   
    L=bu(L,1./512)
    R=bu(R,1./512)
    C = 0.5*(L+R)
    reg = ds.region(C,L,R)
    if 0:
        subset='FULL'
        proj = ds.proj("density","z")
        pw=proj.to_pw()
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
        proj = ds.proj("density","z",data_source=reg,center=C)
        pw=proj.to_pw(center=C)
        print pw.save('p37_trace1_t%04d'%frame)
    if 0:
        proj = ds.proj("density","y",data_source=reg,center=C)
        pw=proj.to_pw(center=C)
        #pw.annotate_magnetic_field()
        #pw.annotate_streamlines('By','Bz')
        #pw.annotate_streamlines('y-velocity','z-velocity',plot_args={'color':'k'})
        print pw.save('philaments_high_1.5_t%04d'%frame)
    #proj = projo.to_frb(1,[512,512])

if 1:
    mask = proj['density']>0
    nmask = proj['density']<=0
    total_zones = 8192
    Nzones = max(R[0]-L[0], R[1]-L[1])*total_zones
    resolution_name = "%d"%Nzones
    Delta = np.zeros(3)
    Delta[0] = np.floor((proj['x'][mask].max()-proj['x'][mask].min())*total_zones)/total_zones
    Delta[1] = np.floor((proj['y'][mask].max()-proj['y'][mask].min())*total_zones)/total_zones
    Delta[2] = np.floor((proj['z'][mask].max()-proj['z'][mask].min())*total_zones)/total_zones
    if 0:
       print "shift"
       Delta=[1.0,1.0,0]
    Delta_zones = [(a,'code_length') for a in Delta]
    width = max(Delta)
    res = 1./Nzones
    frb = proj.to_frb( width, int(Nzones))
    field = 'density' #'level'
    den_damned_square = copy.copy(frb[field].v)
    #mag = copy.copy(frb['magnetic_pressure'].v)
    if 0:
        mask2  = den >0
        mask2o = den <=0
        den[mask2o] = den[mask2].min()
    #make a subset
    if 1:
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
        if field == 'density':
            density_factor = 1000*4.6*3.08e18 #1000 cm^-3 * 4.6 pc * cm/pc
            to_plot = density_factor* np.log10( den)
            cb_label= r'$N[\rm{cm}^{-2}]=N[\rm{code}]*%0.2e$'%density_factor
        elif field == 'level':
            to_plot = den
            cb_label = r'$level$'

        else:
            raise
        """code-heinous, readout-good plot"""
        plt.clf()
        plt.imshow(to_plot, origin='lower', interpolation='nearest',cmap='gray')
        cb=plt.colorbar()
        cb.set_label(cb_label)
        dx = 4.6/8193.
        old_xticks = plt.xticks()[0][1:-1]
        plt.xticks( old_xticks, ["$%0.2f$"%n for n in  old_xticks*dx+L[0] ] )
        plt.xlabel(r'$x[\rm{pc}]\  (\Delta x=$%s$\rm{pc})$'%expform(dx))
        old_yticks = plt.yticks()[0][1:-1]
        plt.yticks( old_yticks, ["$%0.2f$"%n for n in  old_yticks*dx+L[1] ] )
        plt.ylabel(r'$y[\rm{pc}]\  (\Delta x=$%s$\rm{pc})$'%expform(dx))
        outname = 'filament_image_%s_%s_t%04d_%s_%s.png'%(simname, subset, frame,resolution_name, field)
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

if 0:
    distance = 140. #pc
    resolution = 40.#arcsec
    resolution_name = '40arcsec'
    smooth_pc = distance*resolution/206264. #in pc

    #Nzones = nar((R[0]-L[0], R[1]-L[1]))*total_zones
    smooth_px = max(int(smooth_pc/(4.6/total_zones)), 2)
    ef('p37_tools.py')
    blur = blur_image(den,smooth_px)
    plt.clf()
    if 1:
        plt.imshow( np.log10( den), origin='lower', interpolation='nearest')
        plt.savefig('full_a.png')
        plt.clf()
        plt.imshow( np.log10( blur), origin='lower', interpolation='nearest')
        plt.savefig('full_b.png')
        print "full_b.png"
#   hdu = pyfits.PrimaryHDU(blur)
#   hdulist = pyfits.HDUList([hdu])
#   fitsname = 'b02_512_%04d_smoothed_%s_density.fits'%(frame,resolution_name)
#   print fitsname
#   hdulist.writeto(fitsname)
#   Nzones
#   hdulist.close()
    field = 'density'
    fname = 'filament_data_%s_%s_t%04d_%s_%s.txt'%(simname, subset, frame,resolution_name, field)
    outfile = open(fname,'w')
    for n, row in enumerate(blur):
        N = len(row)
        outfile.write("%0.16e "*N%tuple(row))
        outfile.write("\n")
    outfile.close()
    print fname


if 0:
    """ write text"""
    field = 'density'
    fname = 'filament_data_%s_%s_t%04d_%s_%s.txt'%(simname, subset, frame,resolution_name, field)
    outfile = open(fname,'w')
    den = copy.copy(frb['density'].v)
    thisone = copy.copy(frb[field].v)
    for n, full_row in enumerate(thisone):
        row = full_row[ full_row > 0]
        N = len(row)
        if N == 0:
            continue
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
   

