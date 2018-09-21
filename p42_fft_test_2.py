import taxi
aq24 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq24_sto_b',name='p42_aq24',frames=[20])
aq22 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq22_solenoidal_test',name='p42_aq22',frames=[0])
oober = aq24
field = 'acceleration'
frame = 20
ef('p42_helmholtz.py')
debug = 0
dtype = 'float32'
ifft_func = np.fft.ifft
fft_func = np.fft.fft
dx = 1./32**3

if 0:
    nx = [64,64,64]
    x_coord = np.ogrid[0:nx[0],0:nx[1],0:nx[2]]
    x = np.mgrid[0:nx[0],0:nx[1],0:nx[2]]
    vx = np.sin(x[1])

if 0:
    #simple 1d test form http://math.stackexchange.com/questions/740840/derivative-of-function-using-discrete-fourier-transform-matlab
    a = 0
    b = 4*np.pi
    N = 4096.
    dx = (b-a)/N
    t = a + dx*np.arange(N)
    if 0:
        f = np.sin(t)*np.cos(t)
        dfa=np.cos(2*t)
    if 1:
        #wants an 8.  Hm.
        f = np.sin(t)**2
        dfa=np.sin(2*t)
    fftx = fft_func(f)
    k = 2*np.pi/(b-a)*np.arange(-N/2,N/2)/(N)*8
    df = ifft_func(-1j*k*fftx)
    plt.clf()
    plt.plot(t,f,label='f')
    plt.plot(t,df, label='df')
    plt.plot(t,dfa,label='dfa')
    plt.legend(loc=0)
    plt.savefig('p42_fft_2.pdf')



if 1:
    for oober in [aq22]:
        #oober.frames = [0]
        oober.fill()
        #field = 'acceleration'
        #field = 'Driving'
        field = 'velocity'
        all_fields = ['%s-%s'%(s,field) for s in 'xyz']
        #all_fields = ['DrivingField%s'%s for s in '123']
        for tmpfield in all_fields:
            fff = oober.fft(field = tmpfield)
        output = MakeHelmholz_2(oober,oober.frames[0], field,debug=3)
        g=oober.ds.index.grids[0]
        #v1L= g['x-acceleration']
        #v2L= g['y-acceleration']
        #v3L= g['z-acceleration']

        #v1hat_b = fft_func(v1)
        v1hat = output['Vhat'][1]
        v1 = output['V'][1] 
        #v1_recon= ifft_func( output['VhatPhi'][0] + output['VhatA'][0])
        #stat(v1-v1_recon, 'v1 - v1 recon')
        HatDiff = output['Vhat'][0] - (output['VhatPhi'][0] + output['VhatA'][0])
        stat(HatDiff, 'hat diff')

        #kdotv = output['KdotV']
        #diva = g['acceleration_divergence']
        #print 'div vs kdotv', (diva*diva*dx).sum()/(kdotv*np.conj(kdotv)).sum()

if 1:
    power_phi = output['VhatPhi'][0]*np.conj(output['VhatPhi'][0])
    power_phi+= output['VhatPhi'][1]*np.conj(output['VhatPhi'][1])
    power_phi+= output['VhatPhi'][2]*np.conj(output['VhatPhi'][2])
    ff_phi = Filter.FourierFilter(power_phi)
    power_1d_phi = np.array([power_phi[ff_phi.get_shell(bin)].sum() for bin in range(ff_phi.nx)])
    kspace_phi=ff_phi.get_shell_k()
    power_A = output['VhatA'][0]*np.conj(output['VhatA'][0])
    power_A+= output['VhatA'][1]*np.conj(output['VhatA'][1])
    power_A+= output['VhatA'][2]*np.conj(output['VhatA'][2])
    ff_A = Filter.FourierFilter(power_A)
    power_1d_A = np.array([power_A[ff_A.get_shell(bin)].sum() for bin in range(ff_A.nx)])
    kspace_A=ff_A.get_shell_k()
    plt.clf()
    plt.plot(MinK(kspace_phi), power_1d_phi,label='phi')
    plt.plot(MinK(kspace_A), power_1d_A,label='A')
    #plt.plot(MinK(kspace_A), power_1d_phi/(power_1d_phi+power_1d_A),label='A')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('test.pdf')

if 0:
    # this works.
    print ((v1*v1*dx)).sum()/(v1hat*np.conj(v1hat)).sum()
    v1hathat = np.fft.ifftn(v1hat)/dx
    stat(v1-v1hathat, 'inv-fft')



if 0:
    if 1:
        vhat=[]
        oober.fill(frame)
        if field == 'velocity':
            fieldlist=['%s-velocity'%s for s in 'xyz']
        if field == 'acceleration':
            fieldlist=['%s-acceleration'%s for s in 'xyz']
        if field == 'Driving':
            fieldlist=['DrivingField%s'%s for s  in '123']
        for fld in fieldlist:
            oober.fft(field=fld)
        vhat = oober.fft(frame,fieldlist[0],debug=debug)
    if 0:
        nx = vhat.shape
        kvec = np.ogrid[0:nx[0],0:nx[1],0:nx[2]]
        kdotv = vhat*kvec[0]
                
        for cmpt, cmpt_name in enumerate(fieldlist[1:]):
            vhat = oober.fft(frame,cmpt_name,debug=debug)
            kdotv += vhat*kvec[cmpt+1]

        del vhat
        NormK = kvec[0]**2+kvec[1]**2+kvec[2]**2
        NormK[0,0,0]=1.0
        if dtype == 'float32':
            fft_dtype = 'complex64'
        elif dtype == 'float64':
            fft_dtype = 'complex128'
        for cmpt, cmpt_name in enumerate(fieldlist):
            this_set = kvec[cmpt]*kdotv/(NormK)
            filename = '%s/fft_converging-%s.%s'%(oober.product_dir(frame),cmpt_name,dtype)
            fptr = h5py.File(filename,'w')
            fptr.create_dataset('converging-%s'%(cmpt_name),this_set.shape,data=this_set)
            fptr.close()
            print "Created",filename
            vhat = oober.fft(frame,cmpt_name,debug=debug)
            this_set = vhat - this_set
            filename = '%s/fft_solenoidal-%s.%s'%(oober.product_dir(frame),cmpt_name,dtype)
            fptr = h5py.File(filename,'w')
            fptr.create_dataset('solenoidal-%s'%(cmpt_name),this_set.shape,data=this_set)
            fptr.close()
            print "Created",filename

    if 0:
#ConvergingPower(oobername,frame,field,debug,dtype)

        #def HelmholzPower(oober,frame,field,debug=1,dtype='float32'):
#mark_time=time_marker()
        if field == 'velocity':
            fieldlist=['%s-velocity'%s for s in 'xyz']
        if field == 'acceleration':
            fieldlist=['%s-acceleration'%s for s in 'xyz']
        if field == 'Driving':
            fieldlist=['DrivingField%s'%s for s  in '123']
#mark_time = time_marker()
        if dtype == 'float32':
            fft_dtype = 'complex64'
        elif dtype == 'float64':
            fft_dtype = 'complex128'
        else:
            fft_type=dtype
        power = 0
        for i,x in enumerate('xyz'):
            if mark_time is not None:
                mark_time('Start loop %s'%x)
            debug = 200
            Vhat = oober.fft(frame,'converging-%s-%s'%(x,field),num_ghost_zones=-1,debug=debug,dtype=dtype)
            power += (Vhat.conjugate()*Vhat)
            if mark_time is not None:
                mark_time('power addition')
        shell_average(power,oober,frame,'converging-%s'%field,debug,mark_time)
        div_power = power
        power = 0
        for i,x in enumerate('xyz'):
            if mark_time is not None:
                mark_time('Start loop %s'%x)
            debug = 200
            Vhat = oober.fft(frame,'solenoidal-%s-%s'%(x,field),num_ghost_zones=0,debug=debug,dtype=dtype)
            power += (Vhat.conjugate()*Vhat)
            if mark_time is not None:
                mark_time('power addition')
        shell_average(power,oober,frame,'solenoidal-%s'%field,debug,mark_time)
        sol_power = power
