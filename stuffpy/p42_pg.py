from numpy.fft import fftn, ifftn, fftfreq
dx = 0
dy = 1
dz = 2

def deriv2(var,dim):
    """Returns first derivative (2-point central finite difference)"""
         
    sl_m1 = slice(None,-2,None)
    sl_p1 = slice(2,None,None)
    sl_c = slice(1,-1,None)
    ds = 2.0/float(var.shape[0])        
            
    # extend array due to periodic boundary conditions
    tmp = np.concatenate((var[-1:,:,:],var,var[:1,:,:]),axis=0)
    tmp = np.concatenate((tmp[:,-1:,:],tmp,tmp[:,:1,:]),axis=1)
    tmp = np.concatenate((tmp[:,:,-1:],tmp,tmp[:,:,:1]),axis=2)    
            
    if (dim == dx):
        p1 = tmp[sl_p1,sl_c,sl_c]
        m1 = tmp[sl_m1,sl_c,sl_c]
    elif (dim == dy):
        p1 = tmp[sl_c,sl_p1,sl_c]
        m1 = tmp[sl_c,sl_m1,sl_c]
    elif (dim == dz):
        p1 = tmp[sl_c,sl_c,sl_p1]
        m1 = tmp[sl_c,sl_c,sl_m1]
    else:
        print "watch out for dimension!"
                
    del tmp 
            
    return np.array((p1 - m1)/ds)

def div2(varX,varY=None,varZ=None):        
    """Returns divergence (2-point central finite difference)"""
    
    if varY is None:
        varY = varX[1]
        varZ = varX[2]
        varX = varX[0]
        
    xdx = deriv2(varX,dx)
    ydy = deriv2(varY,dy)
    zdz = deriv2(varZ,dz)
    
    return np.array(xdx + ydy + zdz)

def curl2(varX,varY=None,varZ=None):
    """Returns curl (2-point central finite difference)"""
    
    if varY is None:
        varY = varX[1]
        varZ = varX[2]
        varX = varX[0]
    
    xdy = deriv2(varX,dy)
    xdz = deriv2(varX,dz)
    ydx = deriv2(varY,dx)
    ydz = deriv2(varY,dz)
    zdx = deriv2(varZ,dx)
    zdy = deriv2(varZ,dy)
    
    return np.array([zdy - ydz,xdz - zdx,ydx - xdy])

def grad2(var):
    """Returns gradient (2-point central finite difference)"""
    return np.array([deriv2(var,0), deriv2(var,1),deriv2(var,2)])


def getRotFreeField(V):   
    """
    returns the rotation free component of a 3D 3 component vector field
    based on 2nd order finite central differences by solving
    discrete La Place eqn div V = - div (grad phi)
    """
    
    # set up left side in Fourier space
    divV = div2(V)
    FTdivV = fftn(divV)
    
    
    N = FTdivV.shape[1]
    # assumes 3d periodic domain with L = 1
    kx = fftfreq(N)
    # is actually k/N
    kx, ky, kz = np.meshgrid(kx,kx,kx)
    
    # discrete fourier representation of -div grad based on consecutive 2nd order first derivatives 
    denom = -1/2. * N**2. * (np.cos(4.*np.pi*kx) + np.cos(4.*np.pi*ky) + np.cos(4.*np.pi*kz) - 3.)
    # these are 0 in the nominator anyway, so set this to 1 to avoid division by zero
    denom[denom == 0.] = 1.    

    FTdivV /= denom
    
    phi = ifftn(FTdivV).real
    
    return -grad2(phi)


def getSolWeight(Dil,Rot):
    """
    Calculates the solenoidal weight, i.e. the amount of power in the purely 
    solenoidal (source free) part of a vector field """
    
    RotSqr = 0.
    DilSqr = 0.
    
    for i in range(3):
        RotSqr += np.sum(np.abs(fftn(Rot[i]))**2.)
        DilSqr += np.sum(np.abs(fftn(Dil[i]))**2.)   
    
    
    return (RotSqr/(DilSqr + RotSqr))



SolWeights = {
    0:{'values':[], 'dir':'/scratch1/dcollins/Paper42_new_turb/aq30_m2.9_drive0_32',},
  0.5:{'values':[], 'dir':'/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',},
    1:{'values':[], 'dir':'/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',},
      }

plt.clf()
if 0:
    aq32 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',name='aq32',frames=range(0,42,2),fields=['density'])
    aq31 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',name='aq31',frames=range(0,42,2),fields=['density'])
    aq16 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq16_ppm_m9_drive0_noamr_128',name='aq15',frames=range(17),fields=['density'])
    car = aq32
    frame = 40
    car.fill(frame)
    ds=car.ds
    reg = car.get_region(frame)
    vx = reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume')
    vy = reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume')
    vz = reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume')
    px = reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume')
    py = reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume')
    pz = reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume')
    ex = reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume')
    ey = reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume')
    ez = reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume')
    print "%s vx (%0.2e %0.2e %0.2e) vy/vx %0.2f vz/vx %0.2f"%(car.name, vx,vy,vz,vy/vx,vz/vx)
    print "%s px (%0.2e %0.2e %0.2e) py/px %0.2f pz/px %0.2f"%(car.name, px,py,pz,py/px,pz/px)
    print "%s ex (%0.2e %0.2e %0.2e) ey/ex %0.2f ez/ex %0.2f"%(car.name, ex,ey,ez,ey/ex,ez/ex)


if 0:
    aq32 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',name='aq32',frames=range(0,42,2),fields=['density'])
    aq31 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',name='aq31',frames=range(0,42,2),fields=['density'])
    aq16 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq16_m2.9_drive0_noamr_128',name='aq15',frames=range(17),fields=['density'])
    car = aq31
    frame = 40
    car.fill(frame)
    ds=car.ds

    fields = ['%s-acceleration'%s for s in 'xyz']
    #fields = ['DrivingField%s'%s for s in '123']
    cube = ds.covering_grid(0, left_edge=ds.domain_left_edge,
                            dims=ds.domain_dimensions,
                            fields=[("enzo", fields[0]), 
                                    ("enzo", fields[1]),
                                    ("enzo", fields[2])])
    Acc = [cube[("enzo",fields[0])].d,
          cube[("enzo", fields[1])].d,
          cube[("enzo", fields[2])].d]

    Dil = getRotFreeField(Acc)

    Rot = np.array([Acc[0] - Dil[0],
                Acc[1] - Dil[1],
                Acc[2] - Dil[2]])
    
    # safety check
    if (np.abs(curl2(Dil)).max() > 1e-10):
        print "whoopsie, sth went wrong in calculation the dilatational part"
        
    if (np.abs(div2(Rot)).max() > 1e-10):
        print "whoopsie, sth went wrong in calculation the rotational part"            

    RotSqr = 0.
    DilSqr = 0.
    
    for i in range(3):
        RotSqr += np.sum(np.abs(fftn(Rot[i]))**2.)
        DilSqr += np.sum(np.abs(fftn(Dil[i]))**2.)   
    
    
    rat= (RotSqr/(DilSqr + RotSqr))



    import fourier_tools.fourier_filter as Filter
    dpower = 0
    rpower = 0
    for dim in range(3):
        vdilhat = fftn(Dil[dim])
        vrothat = fftn(Rot[dim])
        dpower = vdilhat*np.conj(vdilhat) + dpower
        rpower = vrothat*np.conj(vrothat) + rpower
    ffd = Filter.FourierFilter(dpower)
    ffr = Filter.FourierFilter(rpower)
    power_1d = np.array([dpower[ffd.get_shell(bin)].sum() for bin in range(ffd.nx)])
    power_1r = np.array([rpower[ffr.get_shell(bin)].sum() for bin in range(ffr.nx)])
    kspace=ffd.get_shell_k()
    mk = kspace[ kspace != 0].min()
    plt.clf()
    plt.plot((kspace/mk)[1:], power_1d[1:],label='dil')
    plt.plot((kspace/mk)[1:], power_1r[1:],label='rot')
    plt.xscale('log')
    plt.yscale('log')
    outname = 'p42_power_%s_%04d.pdf'%(car.name,frame)
    plt.legend(loc=0)
    plt.savefig(outname)
    plt.clf()
    plt.plot((kspace/mk)[1:], power_1r[1:]/(power_1r[1:]+power_1d[1:]),label='rot/(rot+dil)')
    plt.xscale('log')
    outname = 'p42_power_ratio_%s_%04d.pdf'%(car.name,frame)
    plt.legend(loc=0)
    plt.savefig(outname)
    print outname




if 0:
    for SolWeight in SolWeights.keys():
        t = []
        ts = yt.load(["%s/DD%04d/data%04d"%(SolWeights[SolWeight]['dir'],frame,frame) for frame in range(0,44,4) ])


        for ds in ts:
            if ds.current_time == 0.:
                continue
                
            cube = ds.covering_grid(0, left_edge=ds.domain_left_edge,
                                    dims=ds.domain_dimensions,
                                    fields=[("enzo", "x-acceleration"), 
                                            ("enzo", "y-acceleration"),
                                            ("enzo", "z-acceleration")])
            Acc = [cube[("enzo", "x-acceleration")].d,
                  cube[("enzo", "y-acceleration")].d,
                  cube[("enzo", "z-acceleration")].d]

            Dil = getRotFreeField(Acc)

            Rot = np.array([Acc[0] - Dil[0],
                        Acc[1] - Dil[1],
                        Acc[2] - Dil[2]])
            
            # safety check
            if (np.abs(curl2(Dil)).max() > 1e-10):
                print "whoopsie, sth went wrong in calculation the dilatational part"
                
            if (np.abs(div2(Rot)).max() > 1e-10):
                print "whoopsie, sth went wrong in calculation the rotational part"            


            t.append(ds.current_time/0.172)
            SolWeights[SolWeight]['values'].append(getSolWeight(Dil,Rot))

        plt.plot(t,SolWeights[SolWeight]['values'],"x-",
             label="$\zeta = %.1f$; data mean = %.2g" % (SolWeight,np.mean(SolWeights[SolWeight]['values'])))

        plt.ylim(-0.05,1.05)
        plt.ylabel("Solenoidal weight")
        plt.xlabel("time t/T")
        plt.legend(loc="lower left", bbox_to_anchor=(0.05,0.1))
    outname = 'p42_pgtest.pdf'
    plt.savefig(outname)
    print outname
