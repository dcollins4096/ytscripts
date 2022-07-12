
from go import *
import performance_tools as pt
reload(pt)
plt.close('all')
if 1:

    #flt1 = taxi.fleet(['g40_i01','g40_i02', 'g40_i02_8', 'g40_i03', 'g40_i04','g40_i05'])
    flt1 = taxi.fleet(['g40_t2_i01','g40_t2_i02', 'g40_t2_i03', 'g40_t2_i04','g40_t2_i05'])
    flt1.name = 'turbulence'

    flt2 = taxi.fleet(['g40_w01', 'g40_w02','g40_w03', 'g40_w04','g40_w05'])
    flt2.name = 'supernovae'
#    flt3 = taxi.fleet(['g40_sa','g40_sb','g40_sc','g40_sd'])
#    flt3 = taxi.fleet(['g40_sc'])
#    flt3.name = 'spheres'
    flt5 = taxi.fleet(['g40_g01','g40_g02','g40_g03', 'g40_g04', 'g40_g05'])
    flt5.name = 'cores'
    flt6 = taxi.fleet(['g40_c01','g40_c02', 'g40_c03'])#,'g40_c04']) #g40_c05
    flt6.name = 'cooling'
#    armada = [flt1,flt2,flt3]
    armada=[flt1,flt5]
#    armada=[flt1,flt2, flt5, flt6]#, flt2,flt4,flt5,flt6]
    pdict = {}
    for flt in armada:
        pdict[flt.name]={}
    
    for flt in armada:
        for car in flt.taxi_list:
            pdict[flt.name][car.name]= pt.perform("%s/performance.out"%car.directory)

lab_mpi=r'$\rm{mpi\ tasks}$'
fig_zu,ax_zu=plt.subplots(1,1)
fig_time, ax_time=plt.subplots(1,1)
areas = ['CommunicationTranspose', 'ComputePotentialFieldLevelZero', 'EL 00', 'EL 01', 'EL 02', 'EL 03', 'EL 04', 'Group WriteAllData', 'Level 00', 'RebuildHierarchy', 'SetBoundaryConditions', 'SolveHydroEquations', 'Total', 'driving', 'mpi_tasks']
skip_plot = ['CommunicationTranspose', 'ComputePotentialFieldLevelZero', 'EL 04', 'driving', 'mpi_tasks']
do_plots=None
#do_plots = ['Total','EH 00','EH 01','EH 02','EH EL', 'mpi_tasks']
#do_plots = ['Total', 'Level 00', 'EL 00','EL 01','EL 02','EL 03','EL small', 'mpi_tasks']
#do_plots += ['SolveHydroEquations', 'SetBoundaryConditions']
if 1:
    for suite_name in pdict:
        print("suite "+suite_name)
        suite = pdict[suite_name]
        mean_time = {}
        for simname in suite:
            if do_plots is None:
                do_plots = list(suite[simname].data.keys())
            for key in do_plots:
                mean_time[key]=[]
        zone_up_su = {}
        for key in ['Total','Level 00']:
            zone_up_su[key]=[]
        for simname in suite:
            sim=suite[simname]
            zone_up_su['Total'].append( sim.data['Total']['Updates/processor/sec'].mean())
            zone_up_su['Level 00'].append( sim.data['Level 00']['Updates/processor/sec'].mean())


            ax_time.plot( sim.data['Total']['Cycle'], sim.data['Total']['Updates/processor/sec'],label=simname,marker='*')
            fig_i, axi = plt.subplots(1,1)

            total_usage = 0
            for key in do_plots:
                if key not in sim.data:
                    continue
                if key == 'mpi_tasks':
                    mean_time[key].append( sim.data[key])
                else:
                    mean_time[key].append( sim.data[key]['Mean Time'].mean())
                    if key not in skip_plot:
                        axi.plot( sim.data[key]['Cycle'], sim.data[key]['Mean Time'],label="%s %s"%(simname,key),marker='*')

                #if key.startswith("EL "): # not in ['Total','Level 00','mpi_tasks']:
                #    total_usage += sim.data[key]['Mean Time'];
                if key not in ['Total','Level 00','mpi_tasks']:
                    total_usage += sim.data[key]['Mean Time'];
                else:
                    print(" skip key ", key)
            
            axi.plot(sim.data['Total']['Cycle'], total_usage, label='Sum',c='k',marker='*')
            axi.legend(loc=0)
            axbonk(axi,xlabel='cycle',ylabel='Mean Time')
            outname = 'plots_to_sort/g40_%s_timing.png'%simname
            print(outname)
            fig_i.savefig(outname)

        for key in ['Level 00']: #zone_up_su:
            #ax_zu.plot(mean_time['mpi_tasks'],zone_up_su[key],label="%s %s"%(key,suite_name), marker='*')
            ax_zu.plot(mean_time['mpi_tasks'],zone_up_su[key],label=suite_name, marker='*')

        fig,ax=plt.subplots(1,1)
        for key in mean_time:
            if len(mean_time[key]) == 0:
                continue
            if key is not 'mpi_tasks':
                ax.plot(mean_time['mpi_tasks'], mean_time[key],label=key,marker='*')
        axbonk(ax,xlabel=lab_mpi,ylabel='time',xscale='log')
        ax.legend(loc=0)
        outname = 'plots_to_sort/g40_timing_%s.png'%suite_name
        print(outname)
        fig.savefig(outname)
        plt.close(fig)


ax_time.legend(loc=0)
axbonk(ax_time,xlabel='cycle',ylabel='up/proc/sec')
fig_time.savefig('plots_to_sort/g40_timing.png')
ax_zu.legend(loc=0)
axbonk(ax_zu,xlabel=lab_mpi, ylabel=r'$\rm{updates/proc/sec}$',xscale='log',yscale='log')
outname = 'plots_to_sort/g49_zoneup.pdf'
fig_zu.savefig(outname)
print(outname)
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
    outname = "plots_to_sort/g40_perf_total_%s.png"%car_base
    plt.savefig(outname)
    print("plot "+outname)
