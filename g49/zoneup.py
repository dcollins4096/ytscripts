
from go import *
import performance_tools as pt
reload(pt)
plt.close('all')
if 1:

    flt1 = taxi.fleet(['g40_t2_i01','g40_t2_i02', 'g40_t2_i03', 'g40_t2_i04','g40_t2_i05'])
    flt1.name = 'cmb, turbulence'

    flt5 = taxi.fleet(['g40_g01','g40_g02','g40_g03', 'g40_g04', 'g40_g05'])
    flt5.name = 'cores, galaxies'

    armada=[flt1,flt5]

    pdict = {}
    for flt in armada:
        pdict[flt.name]={}
    
    for flt in armada:
        for car in flt.taxi_list:
            pdict[flt.name][car.name]= pt.perform("%s/performance.out"%car.directory)

lab_mpi=r'$\rm{mpi\ tasks}$'
fig_zu,ax_zu=plt.subplots(1,1)

if 1:
    for suite_name in pdict:
        print("suite "+suite_name)
        suite = pdict[suite_name]
        zone_up_su = {}
        mpi_tasks=[]
        for key in ['Total','Level 00']:
            zone_up_su[key]=[]
        for simname in suite:
            sim=suite[simname]
            zone_up_su['Total'].append( sim.data['Total']['Updates/processor/sec'].mean())
            zone_up_su['Level 00'].append( sim.data['Level 00']['Updates/processor/sec'].mean())
            mpi_tasks.append( sim.data['mpi_tasks'])


        for key in ['Level 00']: #zone_up_su:
            #ax_zu.plot(mean_time['mpi_tasks'],zone_up_su[key],label="%s %s"%(key,suite_name), marker='*')
            print( suite_name, key)
            print(zone_up_su[key])
            ax_zu.plot(mpi_tasks,zone_up_su[key],label=suite_name, marker='*')



ax_zu.legend(loc=0)
axbonk(ax_zu,xlabel=lab_mpi, ylabel=r'$\rm{(zone\ updates)/(core\ second)}$',xscale='log',yscale='log', ylim=[1e4,2e5])
outname = 'plots_to_sort/g49_zoneup.pdf'
fig_zu.savefig(outname)
print(outname)
plt.close(fig_zu)
