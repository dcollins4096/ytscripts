
from go import *
import performance_tools as pt
reload(pt)
plt.close('all')
if 1:

    flt1 = taxi.fleet([  'g47_turb_N4', 'g47_turb_N8', 'g47_turb_N16', 'g47_turb_N32',])
    flt1.name = 'turbulence'
    #flt2 = taxi.fleet([  'g47_cmb_N8', 'g47_cmb_N16', 'g47_cmb_N32'])
    flt2 = taxi.fleet([   'g47_cmb_N16', 'g47_cmb_N32'])
    flt2.name = 'cmb'
    #flt3 = taxi.fleet([ 'g47_cores_r4','g47_cores_r8','g47_cores_r16','g47_cores_r32'])
    flt3 = taxi.fleet([ 'g47_cores_r4','g47_cores_r8','g47_cores_r16'])
    flt3.name = 'cores'
    flt4 = taxi.fleet([ 'g47_galaxies_u1', 'g47_galaxies_u2', 'g47_galaxies_u3', 'g47_galaxies_u4', 'g47_galaxies_u5'])
    flt4.name = 'galaxies'
    flt7 = taxi.fleet([ 'g47_blank_ppm','g47_blank_mhd','g47_blank_ppm_grav','g47_blank_mhd_grav'])
    flt7.name = 'blanks'

    #armada = [ flt1, flt3, flt2, flt4]
    armada = [ flt1, flt2, flt3, flt4]
    slices = [slice(7,9), slice(1,None), slice(11,None), slice(1,None)]
    #slices = [slice(11,13), slice(1,None), slice(11,None), slice(1,None)]
    #armada = [flt1]

    pdict = {}
    for flt in armada:
        pdict[flt.name]={}
    
    for flt in armada:
        for car in flt.taxi_list:
            pdict[flt.name][car.name]= pt.perform("%s/performance.out"%car.directory)

lab_mpi=r'$\rm{mpi\ tasks}$'
areas = ['CommunicationTranspose', 'ComputePotentialFieldLevelZero', 'EL 00', 'EL 01', 'EL 02', 'EL 03', 'EL 04', 'Group WriteAllData', 'Level 00', 'RebuildHierarchy', 'SetBoundaryConditions', 'SolveHydroEquations', 'Total', 'driving', 'mpi_tasks']
skip_plot = ['CommunicationTranspose', 'ComputePotentialFieldLevelZero', 'EL 04', 'driving', 'mpi_tasks']
do_plots=None
#do_plots = ['Total','Level 00' 'mpi_tasks']
#do_plots = ['Total','EH 00','EH 01','EH 02','EH EL', 'mpi_tasks']
#do_plots = ['Total', 'Level 00', 'EL 00','EL 01','EL 02','EL 03','EL small', 'mpi_tasks']
#do_plots += ['SolveHydroEquations', 'SetBoundaryConditions']

CorePerNode=64
if 0:
    fig_zu,ax_zu=plt.subplots(1,1, figsize=(8,4))
    for suite_name in pdict:
        print("suite "+suite_name)
        suite = pdict[suite_name]
        mean_time = {}
        core_hr_zone_up = {}
        mpi_tasks = []

        for key in ['Total','Level 00']:
            core_hr_zone_up[key]=[]
        for simname in suite:
            sim=suite[simname]
            core_hr_zone_up['Total'].append(    1./sim.data['Total']['Updates/processor/sec'].mean()/3600)
            core_hr_zone_up['Level 00'].append( 1./sim.data['Level 00']['Updates/processor/sec'].mean()/3600)
            mpi_tasks.append( sim.data['mpi_tasks'])
        for key in ['Total']: #core_hr_zone_up:
            print(nar(core_hr_zone_up[key])/CorePerNode)
            ax_zu.plot(nar(mpi_tasks)/CorePerNode,nar(core_hr_zone_up[key])/CorePerNode,label=suite_name, marker='*')
    ax_zu.legend(loc=0)
    #axbonk(ax_zu,xlabel=lab_mpi, ylabel=r'$SU/zone-up$',xscale='log',yscale='log')
    axbonk(ax_zu,xlabel='Nodes', ylabel=r'$SU/zone-up$',xscale='log',yscale='log', ylim=[1e-11,4e-10])
    outname = 'plots_to_sort/g47_zoneup.pdf'
    fig_zu.savefig(outname)
    print(outname)
    plt.close(fig_zu)

if 1:
    fig_tot,ax_tot=plt.subplots(1,1, figsize=(8,4))
    for nsuite,suite_name in enumerate(pdict):
        print("suite "+suite_name)
        suite = pdict[suite_name]
        mean_time = {}
        core_hr_zone_up = {}
        total_wall = {}
        mpi_tasks = []
        routines = ['Total']
        #sl = slice(7,9)
        sl = slices[nsuite]
        if 0:
            routines=['Total','Level 00', 'SolveHydroEquations']
            routines=[ 'Group WriteAllData', 'Level 00', 'SetBoundaryConditions', 'SolveHydroEquations', 'Total']
            el_routines=['EL%s'%s for s in '123456789']
            routines=['Total', 'EL1','EL5', 'EL2','EL3', 'EL4', 'EL6', 'EL7', 'EL9', 'EL8']
            routines+=['EH1','EH2','EH3','EH4']
            skip_plotting = ['EL2','EL6','EL7', 'EL4'] #these are basically zero
            skip_plotting += ['EL8','EL9'] #boundary routines
            skip_plotting += ['EL1','EL5','EL3'] #the major EL routines

        for key in routines:
            core_hr_zone_up[key]=[]
            total_wall[key]=[]
        for simname in suite:
            sim=suite[simname]
            #core_hr_zone_up['Total'].append(    1./sim.data['Total']['Updates/processor/sec'].mean()/3600)
            #core_hr_zone_up['Level 00'].append( 1./sim.data['Level 00']['Updates/processor/sec'].mean()/3600)
            mpi_tasks.append( sim.data['mpi_tasks'])
            for RRR in routines:
                #sl=slice(3,6)
                #sl=slice(6,9)
                subset =  sim.data[RRR]['Mean Time'][sl]
                print(subset)
                total_wall[RRR].append( subset.mean())
        n_node = nar(mpi_tasks)/CorePerNode
        for key in routines: #core_hr_zone_up:
            print(nar(core_hr_zone_up[key])/CorePerNode)
            total_wall[key]=nar(total_wall[key])
            if key in skip_plotting:
                continue
            line=total_wall[key]
            fit = np.polyfit( np.log10(n_node), np.log10(line),1)
            ax_tot.plot( n_node, 10**(-1*np.log10(n_node)+fit[1]),'k--')
            #label=suite_name + " " + key
            #label = '%s %0.2f'%(key, fit[0])
            label = '%s %0.2f'%(suite_name, fit[0])
            ax_tot.plot(n_node,line,label=label, marker='*')
            #ax_tot.plot(mpi_tasks, mpi_tasks, marker='*')


    if 0:
        #
        # Fitting sub-parts
        #
        EL_all = np.stack([total_wall[key] for key in el_routines])
        EL_tot=EL_all.sum(axis=0)
        ax_tot.plot( n_node, EL_tot, label='EL_all %0.2f'%fit[0], marker='*')
        fit = np.polyfit( np.log10(n_node), np.log10(EL_tot),1)
        #fit = np.polyfit( np.log10(n_node), np.log10(total_wall['SolveHydroEquations']),1)
        ax_tot.plot( n_node, 10**(fit[0]*np.log10(n_node)+fit[1]),c='k')
        #ax_tot.plot( n_node, 10**(-1*np.log10(n_node)+fit[1]),c=[0.5]*3)
        #a = np.stack([total_wall[key] for key in [ 'Group WriteAllData', 'Level 00', 'SetBoundaryConditions', 'SolveHydroEquations']])
        #a = np.stack([total_wall[key] for key in [ 'Group WriteAllData',  'SetBoundaryConditions', 'SolveHydroEquations']])
        #ax_tot.plot( n_node, nar(total_wall['Total']) - nar(total_wall['Group WriteAllData']),c='g')
        #ax_tot.plot( n_node, a.sum(axis=0))
        print(fit)
    ax_tot.legend(loc=0)
    axbonk(ax_tot, xscale='log',yscale='log')
    #axbonk(ax_zu,xlabel=lab_mpi, ylabel=r'$SU/zone-up$',xscale='log',yscale='log')
    #axbonk(ax_tot,xlabel='Nodes', ylabel=r'$SU/zone-up$',xscale='log',yscale='log', ylim=[1e-11,4e-10])
    outname = 'plots_to_sort/g47_wall.pdf'
    fig_tot.savefig(outname)
    print(outname)
    plt.close(fig_tot)
