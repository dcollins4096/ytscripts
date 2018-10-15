if 'ef' not in dir():
    execfile('go')
import taxi
a20 = taxi.taxi('taxi_stops/p42_aq20.taxi')
a23 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq23_stochastic_stock',name='p42_aq23')
a24 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq24_stochastic_test_1',name='p42_aq24')
a20.fill(0)
a23.fill(0)
a24.fill(0)
a20.ds.index.grids[0]['acceleration_divergence']
if 1:
    def all_ratios(oober):
        times = []
        div = []
        curl = []
        rel = []
        for n in sorted(oober.frame_dict.keys())[::2]:
            oober.fill(n)
            cube = oober.ds.all_data()
            #divv = np.sum(cube['velocity_divergence']**2)
            #vort = np.sum(cube['vorticity_magnitude']**2)
            divv = np.sum(cube['acceleration_divergence']**2)
            vort = np.sum(cube['Avorticity_magnitude']**2)
            div.append(divv)
            curl.append(vort)
            rel.append( divv/(divv+vort) )
            times.append( oober.ds['InitialTime']/oober.ds['DrivenFlowAutoCorrl'].max() )
        return times, rel, div, vort

    plt.clf()
    #t20, r20, d20, c20 = all_ratios(a20)
    plt.plot(t20, r20, label = 'a20 DFW = %f'%a20.ds['DrivenFlowWeight'])
    #t23, r23, d23, c23 = all_ratios(a23)
    plt.plot(t23, r23, label = 'a23 DFW = %f'%a23.ds['DrivenFlowWeight'])
    #t24, r24, d24, c24 = all_ratios(a24)
    plt.plot(t24, r24, label = 'a24 DFW = %f'%a24.ds['DrivenFlowWeight'])
    ntests = 3
    sub_out = "_%s"*ntests%tuple([a20.outname, a23.outname, a24.outname])
    plt.legend(loc=0)
    plt.xlabel('t')
    plt.ylabel(r'$|\nabla\cdot v|^2/( |\nabla\cdot v|^2+|\nabla\times b|^2) $')
    plt.ylim(-0.1,1.1)
    fname = 'p42_ZetaTest%s.pdf'%sub_out
    plt.savefig(fname)
