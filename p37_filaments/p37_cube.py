if 'ef' not in dir():
    execfile('go')
ef('p37_misc.py')
#frame = 70
if 1:       
    #frame=40
    #frame=70
    ds = yt.load('/scratch1/dcollins/Paper08/B02/512/RS%04d/restart%04d'%(frame,frame)); simname = 'high_field'
    subset = 'full'
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
    zlim = {60:[0.5,10],
            70:[1e-1,1e3]}.get(frame,None)
    zlim=None
    zoom = 1
    if zoom == 0:
        z_extents = [0.0,1.0]
    elif zoom == 1:
        """final zoom"""
        z_extents={80:[0.175,0.255],
                   70:[0.175,0.255],
                   50:[0.175,0.255],
                   60:[0.175,0.255]}.get(frame,[0.0,1.0])
    L[2]=z_extents[0]; R[2]=z_extents[1]
    C = 0.5*(L+R)
    L=bu(L,1./512)
    R=bu(R,1./512)
    C = 0.5*(L+R)
    reg = ds.region(C,L,R)
    if 1:
        """full projection"""
        for axis in ['z']: # ['y']: # 'xyz'[2:]:
            #horizontal = 'yzx'.index(axis) #this worked, but now I don't see why
            #vertical = 'zxy'.index(axis)
            vertical = 'yzx'.index(axis) 
            horizontal = 'zxy'.index(axis)
            ProjL = nar([L[horizontal], L[vertical]])
            ProjR = nar([R[horizontal], R[vertical]])
            print ProjL
            print ProjR
            field = 'density'
            proj = ds.proj(field,axis,data_source=reg,center=C); 
            pw=proj.to_pw(center=C)
            pw.set_origin('domain')
            width = max(ProjR-ProjL)
            pw.set_width(width)
            pw.set_cmap('density','gray')
            if zlim is not None:
                pw.set_zlim('density',zlim[0],zlim[1])
            print pw.save('p37_%s_t%04d_3d_zoom%02d'%(subset,frame, zoom))
            v_horiz, v_vert = '%s-velocity'%'xyz'[horizontal],'%s-velocity'%'xyz'[vertical]
            if 0:
                nzones = int(width*512.)
                frb = proj.to_frb(width, nzones)
                Y,X = np.mgrid[0:nzones:1, 0:nzones:1]
                pw.plots['density'].axes.hold(True)
                pw.plots['density'].axes.streamplot(X,Y,frb[v_horiz].v,frb[v_vert].v)
                pw.plots['density'].axes.hold(False)
                print pw.save('testx.png')

            if 1:
                proj_mass = ds.proj(field,axis,data_source=reg,center=C, weight_field='cell_mass')
                pw_mass=proj_mass.to_pw(center=C)
                #pw.annotate_streamlines_dave(v_horiz, v_vert,
                #                        factor=33,plot_args={'color':'r'},other_proj=pw_mass)
                pw.annotate_streamlines_dave('B%s'%'xyz'[horizontal],'B%s'%'xyz'[vertical],factor=8,plot_args={'color':'b'},
                                             other_proj=pw_mass)
                print pw.save('davetest_mass_ax%s_%04d.png'%(axis,frame)) #'p37_%s_t%04d_3d_zoom%02d_vec'%(subset,frame, zoom))

    if 0:
        level = 1
        dims = [44*2**level]*3
        print dims
        cg = ds.covering_grid(level,L,dims)
        dataset = cg['density'].v
        if 1:
            setname = 'p37_%s_t%04d_cube_level%02d_%s.fits'%(subset,frame,level,'xyz'[axis])
            hdu = pyfits.PrimaryHDU(dataset)
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto(setname)
            hdulist.close()
            print setname

        if 0:
            plt.clf()
            axis=2
            norm = matplotlib.colors.Normalize(1.5,5)
            to_plot = np.log10(np.sum(cg['density'].v,axis=axis))
            plt.imshow(to_plot,origin='lower',interpolation='nearest',cmap='gray',norm=norm)
            plt.colorbar()
            outname = 'p37_%s_t%04d_cube_level%02d_%s'%(subset,frame,level,'xyz'[axis])
            plt.savefig(outname)
            print outname

                 
