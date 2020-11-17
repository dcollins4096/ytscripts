from go import *
import turb_quan 
reload(turb_quan)

def add_nickle_mass(obj):
    def nickle_mass(field,data):
        Ni = data.ds.arr(data['Density_56Ni'],'g/cm**3')
        return Ni*data['cell_volume']
    obj.add_field('NickleMass',function=nickle_mass,units='g',sampling_type='cell')

if 'flt' not in dir():
    flt = taxi.fleet(['p52_441','p52_432','p52_433','p52_434'])
labelmap = {'p52_441':'D22','p52_432':'D26','p52_433':'D27','p52_434':'D28'}
flt['frames']=[0,10,50,100]
flt['derived_fields']={'nickle_mass':add_nickle_mass}

if 0:
    """make"""
    for car in flt.taxi_list:
        car.frames=range(0,110,10)
        qb=turb_quan.quan_box(car)
        qb.load()
        qb()
        qb.dump()

if 1:
    qbset={}
    for car in flt.taxi_list:
        qbset[car.name]=turb_quan.quan_box(car)
        qbset[car.name].load()

    c=rm(ncar)
    fig,axes=plt.subplots(2,1,sharex=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    ax = axes[0]
    ax1 = axes[1]
    for n,car in enumerate(flt.taxi_list):
        this=qbset[car.name].stuff
        argsort = np.argsort(this['t'])
        ax.plot( nar(this['t'])[argsort], nar(this['56Ni'])[argsort],label=labelmap[car.name],c=rm(n))
        if n > 0:
            delta = nar(this['56Ni'])[argsort] - nar(qbset['p52_441'].stuff['56Ni'])[argsort]
            delta /= nar(this['56Ni'])[argsort]
            ax1.plot( nar(this['t'])[argsort], delta,label=car.name,c=rm(n))

    axbonk(ax,xlabel=r'$t[s]$',ylabel=r'$M_{^{56}\rm{Ni}}\ [M_{\odot}]$')
    axbonk(ax1,xlabel=r'$t[s]$',ylabel=r'$\Delta M_{^{56}\rm{Ni}}/M\ \ [M_{\odot}]$')
    ax.legend(loc=0)
    fig.savefig('p52_nickle_time.pdf')
    plt.close(fig)

