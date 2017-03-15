if 'ef' not in dir():
    execfile('go')
def g(n):
    aj01.fill(n)
    aj02.fill(n)
    aj01.ds.print_stats()
    aj02.ds.print_stats()
    for car in [aj01,aj02]:
        print "%s %d"%(car.name, car.ds['NumberOfParticles'])

#taxi_list = [aj01,aj02]
if 0:
    pass
    #aj15=taxi.taxi('aj15')
    #aj16=taxi.taxi('aj16')
    #aj17=taxi.taxi('aj17')
    #fleet  = taxi.fleet([aj15,aj16,aj17])
    #fleet = taxi.fleet(['aj19_sphere','aj20_sphere'])
    #fleet = taxi.fleet(['aj23','aj24'])
    #fleet = taxi.fleet(['aj25','aj26'])
    #fleet = taxi.fleet(['aj27','aj28'])
    #fleet = taxi.fleet(['aj25','aj26','aj27','aj28','aj29','aj30','aj31','aj32'])
    #fleet = taxi.fleet(['aj37','aj38'])
    #fleet = taxi.fleet(['aj43','aj44'])
    #fleet['frames'] = range(15)
    #fleet[4].frames=range(0,500,15) + [500]
    #fleet[5].frames=range(0,500,15) + [500]
    #fleet[6].frames=range(0,500,15) + [500]
    #fleet[7].frames=range(0,500,15) + [500]
    #fleet = taxi.fleet(['aj31','aj41','aj43','aj46'])
    #fleet = taxi.fleet(['aj50','aj51','aj52'])
    #fleet['frames']=range(0,41,2)

#taxi_list = [aj15,aj16,aj17]

if 1:
    t1 = {}
    npart = {}
    ncycle={}
    mass={}
    frames={}
    for car in fleet.taxi_list:

        t1[car.name] = t1.get(car.name, [])
        npart[car.name] = npart.get(car.name,[])
        ncycle[car.name]= ncycle.get(car.name,[])
        mass[car.name]= mass.get(car.name,[])
        frames[car.name] = frames.get(car.name,[])
        car.fill(0)
        for n in car.frames:
            if n in ncycle[car.name]:
                continue
            car.fill(n)
            ncycle[car.name].append(car.ds['InitialCycleNumber'])
            t1[car.name].append(car.ds['InitialTime'])
            this_npart=car.count_particles() 
            npart[car.name].append(this_npart)
            if this_npart > 0:
                mass[car.name].append( (car.ds.all_data()['particle_mass'].in_units('msun')).sum() )
            else:
                mass[car.name].append(0)



if 1:
    #fleet[0].outname = 'aj25_ppm_def'
    #fleet[1].outname = 'aj25_mhd_def'
    #fleet[2].outname = 'aj25_ppm_no_def'
    #fleet[3].outname = 'aj25_mhd_no_def'
    carnames = "all_0_100_%s"%fleet.allnames()
    plt.clf()
    for car in fleet.taxi_list:
        plt.plot(ncycle[car.name],npart[car.name],label=car.outname,marker='x')

    plt.legend(loc=0)
    plt.xlabel('cycle'); plt.ylabel('nparticles')
    outname = 'p33%s_cycle_particles.pdf'%carnames
    plt.savefig(outname)
    print outname

    plt.clf()
    for car in fleet.taxi_list:
        plt.plot(t1[car.name],npart[car.name],label=car.outname,marker='x')
    plt.legend(loc=0)
    plt.xlabel('t'); plt.ylabel('nparticles')
    outname = 'p33%s_time_particles.pdf'%carnames
    plt.savefig(outname)
    print outname

    plt.clf()
    for car in fleet.taxi_list:
        plt.plot(t1[car.name],nar(mass[car.name])/nar(npart[car.name]),label=car.outname,marker='x')
    plt.legend(loc=0)
    plt.xlabel('t'); plt.ylabel('average particle mass [msun]')
    outname = 'p33%s_time_avgmass.pdf'%carnames
    plt.savefig(outname)
    print outname

    plt.clf()
    for car in fleet.taxi_list:
        plt.plot(t1[car.name],mass[car.name],label=car.outname,marker='x')
    plt.legend(loc=0)
    plt.xlabel('t'); plt.ylabel('total particle mass [msun]')
    outname = 'p33%s_time_mass.pdf'%carnames
    plt.savefig(outname)
    print outname

    plt.clf()
    for car in fleet.taxi_list:
        plt.plot(ncycle[car.name],mass[car.name],label=car.outname,marker='x')
    plt.legend(loc=0)
    plt.xlabel('cycle'); plt.ylabel('total particle mass [msun]')
    outname = 'p33%s_cycle_mass.pdf'%carnames
    plt.savefig(outname)
    print outname

if 0:
    x_bounds = [1e12,-1e12]
    y_bounds = [1e12,-1e12]
    z_bounds = [1e12,-1e12]
    for car in fleet.taxi_list:
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
    for car in fleet.taxi_list:
        car.region_type='sphere'
        car.center=[0.5]*3
        car.radius=(150,'kpc')
        car.operation = 'RegionProjection'
        car.callbacks = ['particles']
        car.cmap['density']='gray'
        car.frames=[150]
        car.plot()

if 0:
    for field in aj15.fields:
        min1 = aj15.proj_zlim[field][0]
        max1 = aj15.proj_zlim[field][1]
        for car in fleet.taxi_list:
            min1 = min( [min1, car.proj_zlim[field][0]])
            max1 = max( [max1, car.proj_zlim[field][1]])
        for car in fleet.taxi_list:
            car.proj_zlim[field] = [min1,max1]


