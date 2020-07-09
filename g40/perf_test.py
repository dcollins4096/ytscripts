
from go import *
import performance_tools as pt
reload(pt)
plt.close('all')
if 1:

    flt1 = taxi.fleet(['g40_i01','g40_i02'])
    flt1.name = 'iso'
    flt2 = taxi.fleet(['g40_w01s','g40_w01','g40_w02'])
    flt2.name = 'wd'
    flt3 = taxi.fleet(['g40_sa','g40_sb','g40_sc','g40_sd'])
    flt3.name = 'spheres'
    armada = [flt1,flt2,flt3]
    armada=[flt3]
    pdict = {}
    for flt in armada:
        pdict[flt.name]={}
    
    for flt in armada:
        for car in flt.taxi_list:
            pdict[flt.name][car.name]= pt.perform("%s/performance.out"%car.directory)

lab_mpi=r'$\rm{mpi\ tasks}$'
fig_zu,ax_zu=plt.subplots(1,1)
fig_time, ax_time=plt.subplots(1,1)
if 1:
    for suite_name in pdict:
        print("suite "+suite_name)
        suite = pdict[suite_name]
        mean_time = {}
        for simname in suite:
            for key in suite[simname].data:
                mean_time[key]=[]
        zone_up_su = {}
        for key in ['Total','Level 00']:
            zone_up_su[key]=[]
        for simname in suite:
            sim=suite[simname]
            zone_up_su['Total'].append( sim.data['Total']['Updates/processor/sec'].mean())
            zone_up_su['Level 00'].append( sim.data['Level 00']['Updates/processor/sec'].mean())


            ax_time.plot( sim.data['Total']['Cycle'], sim.data['Total']['Updates/processor/sec'],label=simname)
            fig_i, axi = plt.subplots(1,1)
            for key in sim.data:
                if key == 'mpi_tasks':
                    mean_time[key].append( sim.data[key])
                else:
                    mean_time[key].append( sim.data[key]['Mean Time'].mean())
                    axi.plot( sim.data[key]['Cycle'], sim.data[key]['Mean Time'],label="%s %s"%(simname,key))
            
            axi.legend(loc=0)
            axbonk(axi,xlabel='cycle',ylabel='Mean Time')
            outname = 'g40_%s_timing.png'%simname
            print(outname)
            fig_i.savefig(outname)

        for key in zone_up_su:
            ax_zu.plot(mean_time['mpi_tasks'],zone_up_su[key],label="%s %s"%(key,suite_name))

        fig,ax=plt.subplots(1,1)
        for key in mean_time:
            if key is not 'mpi_tasks':
                ax.plot(mean_time['mpi_tasks'], mean_time[key],label=key)
        axbonk(ax,xlabel=lab_mpi,ylabel='time')
        ax.legend(loc=0)
        outname = 'g40_timing_%s.png'%suite_name
        print(outname)
        fig.savefig(outname)
        plt.close(fig)


ax_time.legend(loc=0)
axbonk(ax_time,xlabel='cycle',ylabel='up/proc/sec')
fig_time.savefig('g40_timing.png')
ax_zu.legend(loc=0)
axbonk(ax_zu,xlabel=lab_mpi, ylabel=r'$\rm{updates/proc/sec}$')
fig_zu.savefig('g40_zoneup.png')
plt.close(fig_time)
plt.close(fig_zu)
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
