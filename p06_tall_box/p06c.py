
if 0:
    fname = '../Enzo/grackle_3/input/CloudyData_UVB=HM2012.h5'
    fptr = h5py.File(fname,'r')
    Lambda=fptr['CoolingRates']['Metals']['Cooling'][:,0,:]
    Gamma=fptr['CoolingRates']['Metals']['Heating'][:,0,:]
    nG=10**fptr['CoolingRates']['Metals']['Heating'].attrs['Parameter1']
    TG=fptr['CoolingRates']['Metals']['Heating'].attrs['Temperature']
    logn=fptr['CoolingRates']['Metals']['Cooling'].attrs['Parameter1']
    n=10**logn
    T=fptr['CoolingRates']['Metals']['Cooling'].attrs['Temperature']
    T_array=np.repeat([T],29,axis=0)
    n_array=np.repeat([n],161,axis=0)
    P_array = n_array.T*T_array
    nL = n_array.T*Lambda #yes this is what you want.

    def closest(density,target_density):
        #n=10**fptr['CoolingRates']['Metals']['Cooling'].attrs['Parameter1']
        dist = np.abs(density - target_density)
        ind = np.where( dist == dist.min() )[0][0]
        return ind


    i15em2 = closest(n,1.5e-2)
    i15em2l = closest(logn,np.log10(1.5e-2))

    plt.clf()
    rm = rainbow_map(len(n))
    for n in range(0,len(n),2):
        a=plt.plot(T, Lambda[n,:],c=rm(n))
    #plt.colorbar()
    plt.xscale('log');plt.yscale('log')
    plt.xlabel(r'$T[K]$')
    plt.ylabel(r'$\Lambda$')
    outname = 'p06_cloudy1.pdf'
    plt.savefig(outname)
    print outname
    plt.clf()
    plt.imshow(np.log(Lambda),interpolation='nearest')
    outname = 'p06_cloudy2.pdf'
    plt.savefig(outname)
    print outname


    #b = Fn(Lambda, 1.5e-2)
    #print b



if 0:
    car=taxi.taxi('fa07')
    car.frames='every 10'
    mins={}
    maxs={}
    times={}
    axes={}
    fig, ax90210 = plt.subplots(2, 1, sharex=True)
    units={'Cooling_Time':lambda x: x.in_units('Myr')}
    nothing = lambda x:x
    for nf,field  in enumerate(['Temperature','Cooling_Time']):
        mins[field]=[]
        maxs[field]=[]
        times[field]=[]
        axes[field] = ax90210[nf]
        for frame in car.return_frames():
            both = car.stat(field,frame=frame)
            this_fun = units.get(field,nothing)
            mins[field].append(this_fun(both['min'][0]))
            maxs[field].append(this_fun(both['max'][0]))
            times[field].append( car.ds.current_time.in_units('Myr'))
        axes[field].plot(times[field],mins[field],label=field)
        axes[field].set_ylabel("%s [%s]"%(field, mins[field][0].units))
        axes[field].set_xlabel('time [Myr]')
    outname = 'p06_cloudytest_%s.pdf'%car.outname
    fig.savefig(outname)
    print outname




