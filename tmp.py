if 'ef' not in dir():
    execfile('go')

reload(taxi)
if 'flt' not in dir() and 1:
    flt=taxi.fleet(['h03','h02','h01'])

for frame in [9,10]: #range(11):
    reg=flt('output.append(car.get_region(%d))'%frame)
    colors = ['r','g','b']
    for field in ['density']: # ['Metal_Density']: #,'Electron_Density']:
        plt.clf()
        for ncar,car in enumerate(flt):
            reg = car.get_region(frame)
            regfield =  reg[field].flatten()
            plt.scatter(reg['radius'].flatten(),regfield,c=colors[ncar],s=5,label=car.name)
            print "MIN", reg[field].min()
        outname='p33_m7_%04d_%s_scatter.png'%(frame,field)
        plt.legend(loc=1)
        plt.yscale('log')
        plt.ylim([5e-28,9e-27])
        plt.savefig(outname)
        print outname


if 0:
    plt.clf()
    fix,ax1=plt.subplots()
    ax2=ax1.twinx()
    c=['r','b','m','c']
    storage = {'h01':{}, 'h02':{}}
    for ncar,car in enumerate(flt.taxi_list[0:]):
        metal = []
        time = []
        starmass = []

        for n in range(0,12,2):
            reg = car.get_region(n)
            metal.append((reg['cell_volume']*reg['Metal_Density']).sum().in_units('msun'))
            time.append(car.ds['InitialTime'])
            if car.count_particles() > 0:
                starmass.append(reg['particle_mass'].sum().in_units('msun'))
            else:
                starmass.append(0)

        storage[car.name]['metal']=nar(metal)
        storage[car.name]['time']=nar(time)
        storage[car.name]['starmass']=nar(starmass)
        #storage[car.name]['starmass']=starmass
        ax1.plot(time,metal,label=car.name + " metal",c=c[ncar], linestyle="-")
        ax2.plot(time,starmass,c=c[ncar+2], linestyle="--",label=car.name + " starmass")

    ax1.set_xlabel('time [Myr]')
    ax1.set_ylabel('metal [Msun]')
    ax2.set_ylabel('Stellar [Msun]')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.legend(loc=0)
    ax2.legend(loc=0)
    plt.savefig('h01h02_Msun.png')
    print metal
#h01.plot(frames=[1],fields=['density'])
