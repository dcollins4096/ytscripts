
from go import *
import performance_tools as pt
reload(pt)

if 1:
    car = taxi.load('p49d_a04')
    filename = "%s/performance.out"%car.directory
    p = pt.perform(filename)
#print(p.data['Stochastic']['Mean Time']/p.data['MHDRK2']['Mean Time'])

    flt1 = taxi.fleet(['g40_i01','g40_i02'])
    flt1.name = 'iso'
    flt2 = taxi.fleet(['g40_w01','g40_w02'])
    flt2.name = 'wd'
    armada = [flt1,flt2]
    pdict = {}
    for flt in armada:
        pdict[flt.name]={}
    
    for flt in armada:
        for car in flt.taxi_list:
            pdict[flt.name][car.name]= pt.perform("%s/performance.out"%car.directory)


if 0:
    ntask=[]
    mean_time = {}
    for key in pdict['g40_i02'].data.keys():
        mean_time[key]=[]
    zone_up_su = {}
    for key in ['Total','Level 00']:
        zone_up_su[key]=[]
    for sim in pdict:
        zone_up_su['Total'].append( pdict[sim].data['Total']['Updates/processor/sec'].mean())
        zone_up_su['Level 00'].append( pdict[sim].data['Level 00']['Updates/processor/sec'].mean())


        for key in pdict[sim].data.keys():
            if key == 'mpi_tasks':
                mean_time[key].append( pdict[sim].data[key])
            else:
                mean_time[key].append( pdict[sim].data[key]['Mean Time'].mean())


if 0:
    lab_mpi=r'$\rm{mpi\ tasks}$'
    fig,ax=plt.subplots(1,1)
    for key in zone_up_su:
        plt.plot(mean_time['mpi_tasks'],zone_up_su[key],label=key)
    ax.legend(loc=0)
    axbonk(ax,xlabel=lab_mpi, ylabel=r'$\rm{updates/proc/sec}$')
    fig.savefig('g40_zoneup.png')

if 0:
    fig,ax=plt.subplots(1,1)
    for key in mean_time:
        if key is not 'mpi_tasks':
            ax.plot(mean_time['mpi_tasks'], mean_time[key],label=key)
    axbonk(ax,xlabel=lab_mpi,ylabel='time')
    ax.legend(loc=0)
    fig.savefig('g40_timing.png')
    plt.close(fig)



if 0:
    smooth_len = 2
    plt.clf()
    y1 = pt.smooth( p1.data['Total']['Mean Time'],smooth_len)
    y2 = pt.smooth( p2.data['Total']['Mean Time'],smooth_len)
    y3 = pt.smooth( p2.data['driving']['Mean Time'],smooth_len)
    plt.plot(p1.data['Total']['Cycle'],y1,label=c1.name)
    plt.plot(p2.data['Total']['Cycle'],y2,label=c2.name)
    plt.plot(p2.data['Total']['Cycle'],y1[:len(y3)]+y3,label=car_base+"_quan est")
    plt.legend(loc=0)
    plt.ylim(0,11)
    outname = "g40_perf_total_%s.png"%car_base
    plt.savefig(outname)
    print("plot "+outname)
