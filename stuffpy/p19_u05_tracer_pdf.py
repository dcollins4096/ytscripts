
if 1:
    for frame in [0,10,20,30,40,50,60,70,80,90,100,110,120,125]:
        scratchdir = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
        scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
        fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)
        ds = yt.load(fname)
        prof_tracer = yt.create_profile(ds.all_data(),('deposit','all_cic'),fields='cell_volume',weight_field=None)
        prof_tracer.set_x_unit('code_density')
        plot_specs = [dict(color='r')]
        prof_gas = yt.create_profile(ds.all_data(),'density',fields='cell_volume',weight_field=None)
        prof_gas.set_x_unit('code_density')
        plot_specs.append(dict(color='g'))
        #plot = yt.ProfilePlot.from_profiles([prof_tracer,prof_gas],plot_specs=plot_specs)
        plt.clf()
        scaler = 1e10 #*(128.**-3)
        plt.plot(0.5*(prof_tracer.x_bins[1:]+prof_tracer.x_bins[0:-1])*scaler,prof_tracer['cell_volume'],c='r',
                 label = 'tracer')
        plt.plot(0.5*(prof_gas.x_bins[1:]+prof_gas.x_bins[0:-1]),prof_gas['cell_volume'],c='g',
                 label = 'gas')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc=1)
        plt.ylabel(r'$V(\rho)$')
        plt.xlabel(r'$\rho$')
        fname   = "prof_u05_%04d.pdf"%frame
        print fname
        print plt.savefig(fname)

