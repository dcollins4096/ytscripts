def g(n):
    aj01.fill(n)
    aj02.fill(n)
    aj01.ds.print_stats()
    aj02.ds.print_stats()
    for car in [aj01,aj02]:
        print "%s %d"%(car.name, car.ds['NumberOfParticles'])

#taxi_list = [aj01,aj02]
if 0:
    aj15=taxi.taxi('aj15')
    aj16=taxi.taxi('aj16')
    aj17=taxi.taxi('aj17')
    fleet  = taxi.fleet([aj15,aj16,aj17])
taxi_list = [aj15,aj16,aj17]
if 0:
    t1 = {}
    npart = {}
    ncycle={}
    mass={}
    for car in taxi_list:
        car.frames=range(0,600)
        t1[car.name] = []
        npart[car.name] = []
        ncycle[car.name]=[]
        mass[car.name]=[]
        car.fill(0)
        for n in car.frames:
            if not car.frame_dict.has_key(n):
                continue
            car.fill(n)
            ncycle[car.name].append(car.ds['InitialCycleNumber'])
            t1[car.name].append(car.ds['InitialTime'])
            npart[car.name].append(car.ds['NumberOfParticles'])
            if car.ds['NumberOfParticles'] > 0:
                mass[car.name].append( (car.ds.all_data()['particle_mass'].in_units('msun')).sum() )
            else:
                mass[car.name].append(0)



if 0:
    plt.clf()
    for car in taxi_list:
        plt.plot(ncycle[car.name],npart[car.name],label=car.name,marker='x')
    plt.legend(loc=0)
    plt.xlabel('cycle'); plt.ylabel('nparticles')
    plt.xlim(0,1000)
    plt.ylim(0,500)
    carnames = "_%s"*len(taxi_list)%tuple([car.name for car in taxi_list])
    outname = 'p33%s_cycle_particles.png'%carnames
    plt.savefig(outname)
    print outname

    plt.clf()
    for car in taxi_list:
        plt.plot(t1[car.name],npart[car.name],label=car.name,marker='x')
    plt.legend(loc=0)
    plt.xlabel('t'); plt.ylabel('nparticles')
    carnames = "_%s"*len(taxi_list)%tuple([car.name for car in taxi_list])
    outname = 'p33%s_time_particles.png'%carnames
    plt.savefig(outname)
    print outname

    plt.clf()
    for car in taxi_list:
        plt.plot(t1[car.name],mass[car.name],label=car.name,marker='x')
    plt.legend(loc=0)
    plt.xlabel('t'); plt.ylabel('total particle mass [msun]')
    carnames = "_%s"*len(taxi_list)%tuple([car.name for car in taxi_list])
    outname = 'p33%s_time_mass.png'%carnames
    plt.savefig(outname)
    print outname

    plt.clf()
    for car in taxi_list:
        plt.plot(ncycle[car.name],mass[car.name],label=car.name,marker='x')
    plt.legend(loc=0)
    plt.xlabel('cycle'); plt.ylabel('total particle mass [msun]')
    carnames = "_%s"*len(taxi_list)%tuple([car.name for car in taxi_list])
    outname = 'p33%s_cycle_mass.png'%carnames
    plt.savefig(outname)
    print outname

if 0:
    x_bounds = [1e12,-1e12]
    y_bounds = [1e12,-1e12]
    z_bounds = [1e12,-1e12]
    for car in taxi_list:
        for frame in [150]:
            car.fill(frame)
            xp = car.ds.all_data()['particle_position_x']
            yp = car.ds.all_data()['particle_position_y']
            zp = car.ds.all_data()['particle_position_z']
            x_bounds[0] = min([x_bounds[0], xp.min()])
            y_bounds[0] = min([y_bounds[0], yp.min()])
            z_bounds[0] = min([z_bounds[0], zp.min()])
            x_bounds[1] = max([x_bounds[1], xp.max()])
            y_bounds[1] = max([y_bounds[1], yp.max()])
            z_bounds[1] = max([z_bounds[1], zp.max()])



if 0:
    for car in taxi_list:
        car.region_type='sphere'
        car.center=[0.5]*3
        car.radius=(150,'kpc')
        car.operation = 'RegionProjection'
        car.callbacks = ['particles']
        car.cmap['density']='gray'
        car.frames=[150]
        car.plot()

if 1:
    for field in aj15.fields:
        min1 = aj15.proj_zlim[field][0]
        max1 = aj15.proj_zlim[field][1]
        for car in fleet.taxi_list:
            min1 = min( [min1, car.proj_zlim[field][0]])
            max1 = max( [max1, car.proj_zlim[field][1]])
        for car in fleet.taxi_list:
            car.proj_zlim[field] = [min1,max1]


