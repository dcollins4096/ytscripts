if 'ef' not in dir():
    execfile('go')

#taxi_list = [aj01,aj02]
if 'fleet' not in dir():
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
    #fleet=taxi.fleet(['g15','g16'])
    #fleet = taxi.fleet(['g17','g18','g17c','g18c'])
    fleet = taxi.fleet(['g17e','g18e','g19e','g20e'])
    #fleet['frames']=range(0,250,10)
    fleet['frames'] = 'every 10'
    particle_pickle_name = "gXXe_particles.pickle"
#taxi_list = [aj15,aj16,aj17]

if 1:
    if 1:
        if 'all_the_stuff' not in dir():
            all_the_stuff=fPickle.load(particle_pickle_name)
            t1    =all_the_stuff['t1']
            npart =all_the_stuff['npart']
            ncycle=all_the_stuff['ncycle']
            mass  =all_the_stuff['mass']
            frames=all_the_stuff['frames']
    else:
        t1 = {}
        npart = {}
        ncycle={}
        mass={}
        frames={}

    for car in fleet:
        t1[car.name] = t1.get(car.name, [])
        npart[car.name] = npart.get(car.name,[])
        ncycle[car.name]= ncycle.get(car.name,[])
        mass[car.name]= mass.get(car.name,[])
        frames[car.name] = frames.get(car.name,[])
        #frames[car.name] = range(0,250,10)
        car.fill(0)
        for n in car.return_frames():
            if n in frames[car.name]:
                continue
            car.fill(n)
            frames[car.name].append(n)
            ncycle[car.name].append(car.ds['InitialCycleNumber'])
            t1[car.name].append(car.ds['InitialTime'])
            this_npart=car.count_particles()
            npart[car.name].append(this_npart)
            if this_npart > 0:
                mass[car.name].append( (car.ds.all_data()['particle_mass'].in_units('msun')).sum() )
            else:
                mass[car.name].append(0)

if 1:
	all_the_stuff={'t1':t1,'npart':npart,'ncycle':ncycle,'mass':mass,'frames':frames}
	all_the_stuff=fPickle.dump(all_the_stuff,particle_pickle_name)


if 1:
    #fleet[0].outname = 'aj25_ppm_def'
    #fleet[1].outname = 'aj25_mhd_def'
    #fleet[2].outname = 'aj25_ppm_no_def'
    #fleet[3].outname = 'aj25_mhd_no_def'
    format = 'png'
    max_frame = max([max(car.return_frames()) for car in fleet])
    min_frame = min([min(car.return_frames()) for car in fleet])
    carnames = ("all_%04d_%04d_%s"%(min_frame,max_frame,fleet.allnames()), format)

    if 0:
        plt.clf()
        for car in fleet.taxi_list:
            plt.plot(ncycle[car.name],npart[car.name],label=car.outname,marker='x')

        plt.legend(loc=0)
        plt.xlabel('cycle'); plt.ylabel('nparticles')
        outname = 'p33%s_cycle_particles.%s'%carnames
        plt.savefig(outname)
        print outname

    if 0:
        plt.clf()
        for car in fleet.taxi_list:
            plt.plot(t1[car.name],npart[car.name],label=car.outname,marker='x')
        plt.legend(loc=0)
        plt.xlabel('t'); plt.ylabel('nparticles')
        outname = 'p33%s_time_particles.%s'%carnames
        plt.savefig(outname)
        print outname

    if 0:
        plt.clf()
        for car in fleet.taxi_list:
            plt.plot(t1[car.name],nar(mass[car.name])/nar(npart[car.name]),label=car.outname,marker='x')
        plt.legend(loc=0)
        plt.xlabel('t'); plt.ylabel('average particle mass [msun]')
        outname = 'p33%s_time_avgmass.%s'%carnames
        plt.savefig(outname)
        print outname

    if 0:
        plt.clf()
        for car in fleet.taxi_list:
            plt.plot(t1[car.name],mass[car.name],label=car.outname,marker='x')
        plt.legend(loc=0)
        plt.xlabel('t'); plt.ylabel('total particle mass [msun]')
        outname = 'p33%s_time_mass.%s'%carnames
        plt.savefig(outname)
        print outname

    if 0:
        plt.clf()
        for car in fleet.taxi_list:
            plt.plot(ncycle[car.name],mass[car.name],label=car.outname,marker='x')
        plt.legend(loc=0)
        plt.xlabel('cycle'); plt.ylabel('total particle mass [msun]')
        outname = 'p33%s_cycle_mass.%s'%carnames
        plt.savefig(outname)
        print outname
        
    if 1:
        dt = {}
        tmid = {}
        dmass = {}
        dmdt = {}
        for car in fleet.taxi_list:
            t_temp = nar(t1[car.name])
            m_temp = nar(mass[car.name])
            dt[car.name] = t_temp[1:] - t_temp[:-1]
            tmid[car.name] = 0.5*(t_temp[1:] + t_temp[:-1])
            dmass[car.name] = m_temp[1:]-m_temp[:-1]
            dmdt[car.name] = dmass[car.name]/(dt[car.name]*1e6)

        plt.clf()
        captions={'g17e':'PPM stock', 'g18e':r'CT, $B_0=10^{-15}\rm{G}$','g19e':r'CT, $B_0=0$','g20e':r'PPM, RK2 timestep'}
        for car in fleet.taxi_list:
            plt.plot(tmid[car.name],dmdt[car.name],label=captions[car.name],marker='x')
        plt.legend(loc=0)
        plt.xlabel('t [Myr]'); plt.ylabel(r'$\dot{M}_*\ [\rm{M}_\odot/\rm{Myr}]$')
        outname = 'p33%s_SFR_time.%s'%carnames
        plt.savefig(outname)
        print outname

        plt.xlim(0,3000)
        plt.ylim(-1,10)
        outname = 'p33%s_restricted_SFR_time.%s'%carnames
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


