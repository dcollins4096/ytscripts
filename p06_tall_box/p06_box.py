

if 0:
    car=taxi.taxi('fg02')
    car.frames=[0,8]
    #car.plot()
    ds = car.load(0)
    ds.print_stats()
    ad=ds.all_data()
    grid=ds.index.grids[0]
    print "density", ad.quantities['Extrema']('density')
    print "Tempera", ad.quantities['Extrema']('Temperature')
    print "TE", ad.quantities['Extrema']('TotalEnergy')

if 1:
    car=taxi.taxi('fg05')
    car.frames=range(0,21,5) #21)
    car.frames=[0,5,20]
    rm = rainbow_map(21)
    plt.clf()
    nfinal = car.return_frames()[-1]
    #fig, (ax_rho, ax_P, ax_gz) = plt.subplots(3, 1) #, sharex=True)
    fig, ((ax_rho, ax_P), (ax_gz,ax_v)) = plt.subplots(2, 2) #, sharex=True)
    for frame in car.return_frames():
        ds = car.load(frame)
        grid=ds.index.grids[0]
        sl=[slice(8,9),slice(8,9),slice(None)]
        z=grid['z'][sl].flatten()
        c = {0:'r',1:'g',nfinal:'b'}.get(frame, 'k')
        c=rm(frame)
        density = grid['density'][sl].flatten()
        dumb_plt(ax_rho, z.flatten(),density,'','density','p06_%s_oned_several.pdf'%(car.outname), c=c)
        dumb_plt(ax_P, z.flatten(),grid['pressure'][sl].flatten(),'z','pressure','p06_%s_oned_several.pdf'%(car.outname), c=c)
        P = grid['pressure'][sl].flatten()
        P = P.in_units('code_mass/(code_length*code_time**2)')
        gradp = ( (P[1:]-P[:-1])/((z[1:]-z[:-1])) ).in_units('code_mass/(code_length**2*code_time**2)')
        zbar = 0.5*(z[1:]+z[:-1])

        dumb_plt(ax_gz, zbar,gradp,'','','p06_%s_oned_several.pdf'%(car.outname), c=c,linestyle="--",label={nfinal:r'$dP/dz$'}.get(frame,None))
        az = grid['External_Acceleration_z'][sl].flatten()
        fz = (density.in_units('code_density')*az).v
        dumb_plt(ax_gz, z.flatten(),fz,'','Fz, gradP','p06_%s_oned_several.pdf'%(car.outname), c=c, label={nfinal:r'$F_g$'}.get(frame,None))
        ax_gz.legend(loc=1)


        vz = grid['velocity_z'][sl].flatten()
        dumb_plt(ax_v, z.flatten(),vz,'','Vz','p06_%s_oned_several.pdf'%(car.outname), c=c,label={nfinal:r'$v_z$'}.get(frame,None))
        cs = grid['sound_speed'][sl].flatten()
        dumb_plt(ax_v, z.flatten(),cs,'','','p06_%s_oned_several.pdf'%(car.outname), c=c,label={nfinal:r'$c_s$'}.get(frame,None), linestyle="--")
        ax_v.legend(loc=1)
    plt.close(fig)

if 0:
    #car=taxi.taxi('fg04')
    car = taxi.taxi('fg05')
    dumpname = car.directory+"/dump"
    stuff = []
    k=[]
    g=[]
    int_g=[]
    max_g=[]
    rho=[]
    rho_c=[]
    shift=[]
    
    import re
    m1 = re.compile(r'CLOWNc \[(.*)\]')
    for line in open(dumpname,'r'):
        if line.startswith("CLOWNc"):
            match= m1.match(line)
            if match is  None:
                print "Shoot, must have broken something."
                raise
            else:
                vals = match.group(1).split(",")
                k1,g1,G1,Gm1,r0,s,r1=vals
                k.append(k1)
                g.append(g1)
                int_g.append(G1)
                max_g.append(Gm1)
                rho_c.append(r0)
                shift.append(s)
                rho.append(r1)
    k=nar(map(int,k))
    g=nar(map(float,g))
    int_g=nar(map(float,int_g))
    max_g=nar(map(float,max_g))
    rho=nar(map(float,rho))
    rho_c=nar(map(float,rho_c))
    density=nar(map(float,density))
    shift=nar(map(float,shift))

    fig, (ax_rho, ax_gz) = plt.subplots(2, 1) #, sharex=True)
    ax_gz.scatter(k, g,c='r',label='g')
    ax_gz.scatter(k, int_g,c='g',label='int(g)')
    ax_gz.set_ylim(min([int_g.min(),int_g.max()]))
    ax_rho.scatter(k, rho,c='b')
    ax_rho.set_ylim(min(rho),max(rho))
    outname = 'p06_dumb1.pdf'
    fig.savefig(outname)
    print outname




"""
fprintf(stderr,"CLOWNc [%d, %"GOUTSYM", %"GOUTSYM", %"GOUTSYM", %"GOUTSYM", %"GOUTSYM", %"GOUTSYM" %"GOUTSYM"],\n",
 k, AccelerationField[2][index], integrated_acceleration[k], max_integrated_accel, CentralDensity, 1.0*shift,
 density);"""
                


#   b=nar(
#       [ [13, 8.009136e-05, 8.009136e-05, 0.00096109632, 0.016579514556299, 0.016570730947099],
#      [14, -8.009136e-05, 0, 0.00096109632, 0.016570730947099, 0.016570730947099],
#      [15, -8.009136e-05, -8.009136e-05, 0.00096109632, 0.0165619473379, 0.016570730947099],
#      [16, -8.009136e-05, -0.00016018272, 0.00096109632, 0.016545861278657, 0.016570730947099],
#      [17, -8.009136e-05, -0.00024027408, 0.00096109632, 0.01652501462486, 0.016570730947099],
#      [18, -8.009136e-05, -0.00032036544, 0.00096109632, 0.016500315654888, 0.016570730947099],
#      [19, -8.009136e-05, -0.0004004568, 0.00096109632, 0.016472289607292, 0.016570730947099],
#      [20, -8.009136e-05, -0.00048054816, 0.00096109632, 0.016441291015211, 0.016570730947099],
#      [21, -8.009136e-05, -0.00056063952, 0.00096109632, 0.016407580232291, 0.016570730947099],
#      [22, -8.009136e-05, -0.00064073088, 0.00096109632, 0.016371359033937, 0.016570730947099],
#      [23, -8.009136e-05, -0.00072082224, 0.00096109632, 0.016332789791971, 0.016570730947099],
#      [24, -8.009136e-05, -0.0008009136, 0.00096109632, 0.016292006856983, 0.016570730947099],
#      [25, -8.009136e-05, -0.00088100496, 0.00096109632, 0.016249123804613, 0.016570730947099],
#      [12, 8.009136e-05, 0, 0.00096109632, 0.016570730947099, 0.016570730947099],
#      [11, 8.009136e-05, -8.009136e-05, 0.00096109632, 0.0165619473379, 0.016570730947099],
#      [10, 8.009136e-05, -0.00016018272, 0.00096109632, 0.016545861278657, 0.016570730947099],
#      [9, 8.009136e-05, -0.00024027408, 0.00096109632, 0.01652501462486, 0.016570730947099],
#      [8, 8.009136e-05, -0.00032036544, 0.00096109632, 0.016500315654888, 0.016570730947099],
#      [7, 8.009136e-05, -0.0004004568, 0.00096109632, 0.016472289607292, 0.016570730947099],
#      [6, 8.009136e-05, -0.00048054816, 0.00096109632, 0.016441291015211, 0.016570730947099],
#      [5, 8.009136e-05, -0.00056063952, 0.00096109632, 0.016407580232291, 0.016570730947099],
#      [4, 8.009136e-05, -0.00064073088, 0.00096109632, 0.016371359033937, 0.016570730947099],
#      [3, 8.009136e-05, -0.00072082224, 0.00096109632, 0.016332789791971, 0.016570730947099],
#      [2, 8.009136e-05, -0.0008009136, 0.00096109632, 0.016292006856983, 0.016570730947099],
#      [1, 8.009136e-05, -0.00088100496, 0.00096109632, 0.016249123804613, 0.016570730947099],
#      [0, 8.009136e-05, -0.00096109632, 0.00096109632, 0.016204238299692, 0.016570730947099] ])


if 0:
    plt.clf()
    fig, (ax_rho, ax_gz) = plt.subplots(2, 1) #, sharex=True)
    ax_gz.scatter(b[:,0], b[:,1],c='r',label='g')
    ax_gz.scatter(b[:,0], b[:,2],c='g',label='int(g)')
    ax_gz.set_ylim(min(b[:,2]),max(b[:,2]))
    ax_rho.scatter(b[:,0], b[:,4],c='b')
    ax_rho.set_ylim(min(b[:,4]),max(b[:,4]))
    fig.savefig('p06_dumb1.pdf')

