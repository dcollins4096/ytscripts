from go import *
import cooling_rate

if 1:
    from importlib import reload
    reload(cooling_rate)
    plt.clf()
    density=1.6737352238051868e-24
    Cloudy1 = "/Users/dcollins/Enzo/grackle_3_py3/input/"+"CloudyData_UVB=HM2012.h5"
    my_chem_1 = cooling_rate.make_chem(Cloudy1)
    print(my_chem_1.primordial_chemistry)
    print(my_chem_1.grackle_data_file)
    data_1=cooling_rate.plot_stuff(my_chem_1, plot_object=plt,color="k")
    my_chem_2 = cooling_rate.make_chem(Cloudy1,primordial_chemistry=0)
    data_2=cooling_rate.plot_stuff(my_chem_2,plot_object=plt,color='r')
    plt.savefig('cooling_rate.png')

density=1.6737352238051868e-24
if 1:
    fname = '../Enzo/grackle_3/input/CloudyData_UVB=HM2012.h5'
    fptr = h5py.File(fname,'r')
    Lambda_all=fptr['CoolingRates']['Metals']['Cooling'][:]
    Gamma_all=fptr['CoolingRates']['Metals']['Heating'][:]
    Lambda = Lambda_all[:,0,:]
    Gamma = Gamma_all[:,0,:]
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


    i15em2 = closest(n,density)
    ni = i15em2
#   grackle_data="/home/dcollins4096/Enzo/grackle_3/input/"+"CloudyData_UVB=HM2012.h5"
#   my_chem_2 = cooling_rate.make_chem(grackle_data,primordial_chemistry=0)
#   cooling_rate.plot_stuff(my_chem_2, density = n[ni], plot_object=plt)#,color=rm(b))
    #i15em2l = closest(logn,np.log10(1.5e-2))
    rm = rainbow_map(26*29)
#    b = 0
    #for ni in range(0,29,5):
    for z in [0,-1]: #range(0,26,5):
        print("DENSITY", n[ni])
        b=z+29*ni
#           #a=plt.plot(T, Lambda_all[ni,z,:],c=rm(b))
#           #a=plt.plot(T, Gamma_all[ni,z,:],c=rm(b),linestyle='--')
        annet2 =  Lambda_all[ni,z,:]- Gamma_all[ni,z,:]
        heat = annet2>0
        cool = annet2<0
        plt.plot(T,Lambda_all[ni,z,:], c='y') #- Gamma_all[ni,z,:]
        #a=plt.plot(T[heat],annet2[heat],c=rm(b))
        #a=plt.plot(T[heat],annet2[heat],c=rm(b))
        a=plt.plot(T[cool],np.abs(annet2[cool]),c=rm(b),linestyle='--')
        a=plt.plot(T[cool],np.abs(annet2[cool]),c=rm(b),linestyle='--')
#           #anette= Lambda_all[ni,z,:]* n[ni] - Gamma_all[ni,z,:]
#           #a=plt.plot(T, anette,c=rm(b),linestyle='--')
    plt.xscale('log'); plt.yscale('log')
    outname = 'cooling_rate.png'
    plt.savefig(outname)
    print(outname)
if 0:
    rm = rainbow_map(26)
    for z in range(0,26,2):
        a=plt.plot(T, Lambda_all[i15em2,z,:],c=rm(z))
    #a=plt.plot(T, Lambda_all[i15em2,0,:],c='g')
    plt.savefig('b.png')

if 0:
    #plt.clf()
    rm = rainbow_map(len(n))
    for n in range(0,len(n),2):
        a=plt.plot(T, Lambda[n,:],c=rm(n))
    #plt.colorbar()
    plt.xscale('log');plt.yscale('log')
    plt.xlabel(r'$T[K]$')
    plt.ylabel(r'$\Lambda$')
    outname = 'p06b_cloudy1.pdf'
    plt.savefig(outname)
    print(outname)
    plt.clf()
    p1 = plt.imshow(np.log10(Lambda),interpolation='nearest')
    outname = 'p06b_cloudy2.pdf'
    cb=plt.colorbar(p1,orientation='horizontal')
    cb.set_label('Cooling')
    plt.savefig(outname)
    print(outname)
    plt.clf()
    p2=plt.imshow(np.log10(Gamma),interpolation='nearest')
    outname = 'p06b_cloudy3.pdf'
    cb2=plt.colorbar(p2,orientation='horizontal')
    cb2.set_label('Heating')
    plt.savefig(outname)
    print(outname)


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
    print(outname)




